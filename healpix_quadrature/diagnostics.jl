# Diagnostics for the quadrature stability experiments, see
# docs/dev/2026-08/healpix-quadrature-exactness.md, "Verifying `PerOrderQuadrature` on a GPU node".
#
# Blow-up is the least informative signal, so this records, on a schedule, the quantities that
# discriminate earlier: the vorticity power spectrum (a shallowing or upturning tail is the first
# sign of an ill-behaved transform), the quantities that should be conserved, and the extrema.
#
# Included from `run_case.jl`; not part of the package.

using SpeedyWeather
using SpeedyWeather: AbstractCallback, Schedule, isscheduled
import SpeedyWeather: initialize!, callback!, finalize!
using SpeedyWeather.RingGrids
using SpeedyWeather.SpeedyTransforms
using SpeedyWeather.LowerTriangularArrays
using Dates

"""Per-gridpoint solid angle `ΔΩ_ij`, as a `Field` on the grid, summing to 4π."""
function solid_angle_field(grid, NF)
    ΔΩ = RingGrids.get_solid_angles(grid)            # one entry per ring of one hemisphere
    nlat = RingGrids.get_nlat(grid)
    field = zeros(NF, grid)
    for (j, ring) in enumerate(RingGrids.eachring(grid))
        j_hemisphere = min(j, nlat - j + 1)          # rings are symmetric about the equator
        for ij in ring
            field[ij] = ΔΩ[j_hemisphere]
        end
    end
    return field
end

"""Per-gridpoint `cos(latitude)`, as a `Field` on the grid."""
function coslat_field(grid, NF)
    coslat = cosd.(RingGrids.get_latd(grid))
    field = zeros(NF, grid)
    for (j, ring) in enumerate(RingGrids.eachring(grid))
        for ij in ring
            field[ij] = coslat[j]
        end
    end
    return field
end

"""
Records, on `schedule`, the diagnostics that discriminate between quadrature schemes long before
a run diverges. All reductions run on the architecture the model runs on; only the vorticity
coefficients are copied to the host, for the power spectrum.
"""
Base.@kwdef mutable struct QuadratureDiagnostics{NF, FieldType} <: AbstractCallback
    "[OPTION] how often to record"
    schedule::Schedule = Schedule(every = Day(1))

    "per-gridpoint solid angle ΔΩ_ij, Σ = 4π"
    area::FieldType

    "per-gridpoint ΔΩ_ij cos(φ_ij), for the angular momentum integral"
    area_coslat::FieldType

    "layer thicknesses Δσ_k, summing to 1"
    Δσ::Vector{NF} = NF[]

    time::Vector{DateTime} = DateTime[]

    "the l = m = 0 coefficient of ln(pₛ), the mode the quadrature constrains exactly"
    mean_lnpressure::Vector{Float64} = Float64[]

    "∫ pₛ dΩ / 4π [Pa], the global mean surface pressure, i.e. the total mass"
    mean_pressure::Vector{Float64} = Float64[]

    "Σ_lm Δσ (|ζ_lm|² + |D_lm|²) / (l(l+1)) · R²/2, the rotational + divergent kinetic energy"
    energy::Vector{Float64} = Float64[]

    "Σ_lm Δσ |ζ_lm|² / 2, the potential enstrophy of the barotropic-equivalent flow"
    enstrophy::Vector{Float64} = Float64[]

    "∫ u a cos(φ) dσ dΩ, the relative angular momentum"
    angular_momentum::Vector{Float64} = Float64[]

    max_vorticity::Vector{Float64} = Float64[]
    max_vertical_velocity::Vector{Float64} = Float64[]
    min_pressure::Vector{Float64} = Float64[]
    max_humidity::Vector{Float64} = Float64[]

    "vorticity power spectrum, one column per record, summed over layers"
    spectrum::Vector{Vector{Float64}} = Vector{Float64}[]

    "[OPTION] also keep the surface-layer vorticity coefficients every `snapshot_every` records,
    so two runs on the same spectrum can be differenced directly. 0 disables."
    snapshot_every::Int = 0

    "surface-layer vorticity coefficients, one column per snapshot"
    snapshots::Vector{Vector{ComplexF64}} = Vector{ComplexF64}[]

    "times the snapshots were taken at"
    snapshot_time::Vector{DateTime} = DateTime[]
end

function QuadratureDiagnostics(spectral_grid::SpectralGrid; kwargs...)
    (; NF, grid, architecture) = spectral_grid
    # the weight fields are filled by scalar indexing, so build them on a host copy of the grid
    host_grid = RingGrids.nonparametric_type(typeof(grid))(grid.nlat_half, SpeedyWeather.CPU())
    area = on_architecture(architecture, solid_angle_field(host_grid, NF))
    area_coslat = on_architecture(
        architecture, solid_angle_field(host_grid, NF) .* coslat_field(host_grid, NF)
    )
    return QuadratureDiagnostics{NF, typeof(area)}(; area, area_coslat, kwargs...)
end

function initialize!(callback::QuadratureDiagnostics, vars, model)
    initialize!(callback.schedule, vars.prognostic.clock)
    callback.Δσ = Vector{eltype(callback.Δσ)}(model.geometry.σ_levels_thick)
    n = callback.schedule.steps + 1                 # + 1 for the initial conditions
    for field in (
            :mean_lnpressure, :mean_pressure, :energy, :enstrophy, :angular_momentum,
            :max_vorticity, :max_vertical_velocity, :min_pressure, :max_humidity,
        )
        empty!(getfield(callback, field))
        sizehint!(getfield(callback, field), n)
    end
    empty!(callback.time); sizehint!(callback.time, n)
    empty!(callback.spectrum); sizehint!(callback.spectrum, n)
    empty!(callback.snapshots); empty!(callback.snapshot_time)
    record!(callback, vars, model)                  # record the initial conditions
    return nothing
end

function callback!(callback::QuadratureDiagnostics, vars, model)
    isscheduled(callback.schedule, vars.prognostic.clock) || return nothing
    record!(callback, vars, model)
    return nothing
end

finalize!(::QuadratureDiagnostics, args...) = nothing

"""Take one record. Vorticity and divergence are radius-scaled inside `run!`, so both are
unscaled here to keep the recorded numbers comparable across resolutions and runs."""
function record!(callback::QuadratureDiagnostics, vars, model)
    (; time_stepping) = model
    radius = Float64(model.geometry.radius[])
    scale = vars.prognostic.scale[]                 # vorticity/divergence are scaled by radius
    (; Δσ, area, area_coslat) = callback

    vorticity = SpeedyWeather.get_prognostic_step(vars.prognostic.vorticity, time_stepping, callback, model)
    divergence = SpeedyWeather.get_prognostic_step(vars.prognostic.divergence, time_stepping, callback, model)
    lnpressure = SpeedyWeather.get_prognostic_step(vars.prognostic.pressure, time_stepping, callback, model)

    push!(callback.time, vars.prognostic.clock.time)

    # THE EXACTLY CONSTRAINED MODE. Only l = m = 0 is pinned by Σ_j g_j = 4π, so any drift here
    # beyond roundoff is a bug rather than a property of the flow. Y₀₀ = 1/(2√π).
    lnps00 = Complex{Float64}(on_architecture(SpeedyWeather.CPU(), lnpressure)[1])
    push!(callback.mean_lnpressure, real(lnps00) / (2sqrt(π)))

    # SPECTRAL DIAGNOSTICS, on the host: small arrays, and `power_spectrum` indexes scalars
    vorticity_host = on_architecture(SpeedyWeather.CPU(), vorticity)
    vorticity_host .*= inv(scale)
    divergence_host = on_architecture(SpeedyWeather.CPU(), divergence)
    divergence_host .*= inv(scale)

    power = SpeedyTransforms.power_spectrum(vorticity_host, normalize = false)
    push!(callback.spectrum, Float64.(vec(sum(power .* Δσ', dims = 2))))

    if callback.snapshot_every > 0 &&
            (length(callback.time) - 1) % callback.snapshot_every == 0
        push!(callback.snapshots, ComplexF64.(vec(vorticity_host.data[:, end])))
        push!(callback.snapshot_time, vars.prognostic.clock.time)
    end

    power_divergence = SpeedyTransforms.power_spectrum(divergence_host, normalize = false)
    enstrophy = 0.5 * sum(Float64.(power) .* Δσ')
    # (ζ, D) -> kinetic energy: the inverse Laplacian brings in R²/(l(l+1)), l 0-based
    l = 0:(size(power, 1) - 1)
    inverse_eigenvalue = [li == 0 ? 0.0 : radius^2 / (li * (li + 1)) for li in l]
    energy = 0.5 * sum((Float64.(power) .+ Float64.(power_divergence)) .* inverse_eigenvalue .* Δσ')
    push!(callback.enstrophy, enstrophy)
    push!(callback.energy, energy)

    # GRID DIAGNOSTICS, reduced on the device
    u = SpeedyWeather.get_prognostic_step(vars.grid.u, time_stepping, callback, model)
    lnpressure_grid = SpeedyWeather.get_prognostic_step(vars.grid.pressure, time_stepping, callback, model)
    humidity = SpeedyWeather.get_prognostic_step(vars.grid.humidity, time_stepping, callback, model)
    vorticity_grid = SpeedyWeather.get_prognostic_step(vars.grid.vorticity, time_stepping, callback, model)
    w = vars.dynamics.w

    # ∫ u a cos(φ) dσ dΩ; the layer sum is a matrix-vector product over the flattened grid
    angular_momentum = 0.0
    for k in eachindex(Δσ)
        angular_momentum += Float64(Δσ[k]) * Float64(sum(view(u.data, :, k) .* area_coslat.data))
    end
    push!(callback.angular_momentum, radius * angular_momentum)

    # `vars.grid.pressure` holds ln(pₛ); the mass integral needs pₛ itself
    push!(callback.mean_pressure, Float64(sum(exp.(lnpressure_grid.data) .* area.data)) / 4π)
    push!(callback.min_pressure, exp(Float64(minimum(lnpressure_grid.data))))
    push!(callback.max_vorticity, Float64(maximum(abs, vorticity_grid.data)) / scale)
    push!(callback.max_vertical_velocity, Float64(maximum(abs, w.data)))
    push!(callback.max_humidity, Float64(maximum(humidity.data)))
    return nothing
end

"""Everything recorded, as a plain `NamedTuple` of arrays for JLD2."""
function as_namedtuple(callback::QuadratureDiagnostics)
    spectrum = isempty(callback.spectrum) ? zeros(0, 0) : reduce(hcat, callback.spectrum)
    return (;
        time = callback.time,
        mean_lnpressure = callback.mean_lnpressure,
        mean_pressure = callback.mean_pressure,
        energy = callback.energy,
        enstrophy = callback.enstrophy,
        angular_momentum = callback.angular_momentum,
        max_vorticity = callback.max_vorticity,
        max_vertical_velocity = callback.max_vertical_velocity,
        min_pressure = callback.min_pressure,
        max_humidity = callback.max_humidity,
        spectrum,
        snapshot_time = callback.snapshot_time,
        snapshots = isempty(callback.snapshots) ? zeros(ComplexF64, 0, 0) : reduce(hcat, callback.snapshots),
    )
end
