abstract type AbstractVerticalDiffusion <: AbstractParameterization end

export BulkRichardsonDiffusion
@parameterized @kwdef struct BulkRichardsonDiffusion{NF, VectorType} <: AbstractVerticalDiffusion
    "[OPTION] von Kármán constant [1]"
    von_Karman::NF = 0.4

    "[OPTION] roughness length [m]"
    @param roughness_length::NF = 3.21e-5 (bounds = Positive,)

    "[OPTION] Critical Richardson number for stable mixing cutoff [1]"
    @param critical_Richardson::NF = 1 (bounds = Positive,)

    "[OPTION] Fraction of surface boundary layer"
    @param surface_layer_fraction::NF = 0.1 (bounds = 0 .. 1,)

    "[OPTION] diffuse static energy?"
    diffuse_static_energy::Bool = true

    "[OPTION] diffuse momentum?"
    diffuse_momentum::Bool = true

    "[OPTION] diffuse humidity? Ignored for PrimitiveDryModels"
    diffuse_humidity::Bool = true

    "[DERIVED] Vertical Laplace operator, operator for cells above"
    ∇²_above::VectorType

    "[DERIVED] Vertical Laplace operator, operator for cells below"
    ∇²_below::VectorType
end

Adapt.@adapt_structure BulkRichardsonDiffusion

# generator function
function BulkRichardsonDiffusion(SG::SpectralGrid; kwargs...)
    arch = SG.architecture
    ∇²_above = on_architecture(arch, zeros(SG.NF, SG.nlayers))
    ∇²_below = on_architecture(arch, zeros(SG.NF, SG.nlayers))
    return BulkRichardsonDiffusion{SG.NF, SG.VectorType}(; ∇²_above, ∇²_below, kwargs...)
end

variables(::BulkRichardsonDiffusion) = (
    ParameterizationVariable(:boundary_layer_height, Grid2D(), desc = "Boundary layer height index", units = "1"),
    ScratchVariable(:vertical_diffusion_c, GridXYZ(), desc = "Scratch array for the tridiagonal solver", units = "?", namespace = :grid),
    ScratchVariable(:vertical_diffusion_d, GridXYZ(), desc = "Scratch array for the tridiagonal solver", units = "?", namespace = :grid),
)

function initialize!(diffusion::BulkRichardsonDiffusion, model::PrimitiveEquation)

    (; nlayers) = model.geometry
    nlayers == 1 && return nothing     # no diffusion for 1-layer model

    # ∇² operator on σ levels like 1/Δσ² but for variable Δσ
    # also includes a 1/2 so that the diffusion coefficients on full levels can be added
    # which is equivalent to interpolating them on half levels for a ∂σ (K ∂σ) formulation
    # with σ-dependent diffusion coefficient K
    # TODO: vertical diffusion coefficients are computed in sigma-coordinate space using
    # σ[k] - σ[k±1] and σ_half spacings. Generalising to hybrid coordinates would require
    # reformulating the ∂σ(K ∂σ) operator in terms of pressure.
    σ = on_architecture(CPU(), model.geometry.σ_levels_full)
    σ_half = on_architecture(CPU(), model.geometry.σ_levels_half)
    ∇²_above = on_architecture(CPU(), diffusion.∇²_above)
    ∇²_below = on_architecture(CPU(), diffusion.∇²_below)

    for k in 1:nlayers
        σ₋ = k <= 1 ? -Inf : σ[k - 1]   # sets the gradient across surface and top to 0
        σ₊ = k >= nlayers ? Inf : σ[k + 1]   # = no flux boundary conditions
        ∇²_above[k] = inv(2 * (σ[k] - σ₋) * (σ_half[k + 1] - σ_half[k]))
        ∇²_below[k] = inv(2 * (σ₊ - σ[k]) * (σ_half[k + 1] - σ_half[k]))
    end

    arch = model.architecture
    diffusion.∇²_above .= on_architecture(arch, ∇²_above)
    diffusion.∇²_below .= on_architecture(arch, ∇²_below)
    return nothing
end

# function barrier
@propagate_inbounds parameterization!(ij, vars, diffusion::BulkRichardsonDiffusion, model) =
    vertical_diffusion!(ij, vars, diffusion, model.time_stepping, model.atmosphere, model.planet, model.orography, model.land_sea_mask, model.geopotential, model.geometry.σ_levels_full)

@propagate_inbounds function vertical_diffusion!(
        ij,
        vars,
        diffusion::BulkRichardsonDiffusion,
        time_stepping,
        atmosphere,
        planet,
        orography,
        land_sea_mask,
        geopot,
        σ_levels_full,
    )

    (; diffuse_momentum, diffuse_static_energy, diffuse_humidity) = diffusion

    # escape immediately if all diffusions disabled
    any((diffuse_momentum, diffuse_static_energy, diffuse_humidity)) || return nothing

    K, kₕ = get_diffusion_coefficients!(ij, vars, diffusion, time_stepping, atmosphere, planet, orography, land_sea_mask, geopot, σ_levels_full)

    u_tend = get_tendency_step(vars.tendencies.grid.u, time_stepping, diffusion)
    v_tend = get_tendency_step(vars.tendencies.grid.v, time_stepping, diffusion)
    temp_tend = get_tendency_step(vars.tendencies.grid.temperature, time_stepping, diffusion)

    u = get_prognostic_step(vars.grid.u, time_stepping, diffusion)
    v = get_prognostic_step(vars.grid.v, time_stepping, diffusion)

    # implicit (backward Euler) time step, unconditionally stable
    Δt = implicit_vertical_diffusion_time_step(time_stepping)
    diffuse_momentum && _vertical_diffusion!(ij, vars, u_tend, u, K, kₕ, diffusion, Δt)
    diffuse_momentum && _vertical_diffusion!(ij, vars, v_tend, v, K, kₕ, diffusion, Δt)

    if atmosphere isa AbstractWetAtmosphere && diffuse_humidity
        humid_tend = get_tendency_step(vars.tendencies.grid.humidity, time_stepping, diffusion)
        humid = get_prognostic_step(vars.grid.humidity, time_stepping, diffusion)
        _vertical_diffusion!(ij, vars, humid_tend, humid, K, kₕ, diffusion, Δt)
    end

    if diffuse_static_energy
        # compute dry static energy on the fly
        dry_static_energy = vars.scratch.grid.a
        cₚ = atmosphere.heat_capacity
        T = get_prognostic_step(vars.grid.temperature, time_stepping, diffusion)
        Φ = vars.dynamics.geopotential

        for k in 1:size(T, 2)
            dry_static_energy[ij, k] = cₚ * T[ij, k] + Φ[ij, k]
        end

        # diffuse dry static energy but convert its tendency back to temperature with 1/cₚ
        _vertical_diffusion!(ij, vars, temp_tend, dry_static_energy, K, kₕ, diffusion, Δt, inv(cₚ))
    end
    return nothing
end

@propagate_inbounds function get_diffusion_coefficients!(
        ij,
        vars,
        diffusion::BulkRichardsonDiffusion,
        time_stepping::AbstractTimeStepper,
        atmosphere::AbstractAtmosphere,
        planet::AbstractPlanet,
        orog,
        land_sea_mask,
        geopot::AbstractGeopotential,
        σ_levels_full,
    )
    nlayers = length(diffusion.∇²_above)

    # parameters
    Ri_c = diffusion.critical_Richardson
    fb = diffusion.surface_layer_fraction
    κ = diffusion.von_Karman
    z₀ = diffusion.roughness_length
    gravity⁻¹ = inv(planet.gravity)

    # Typical height Z of lowermost layer from geopotential of reference surface temperature
    # minus surface geopotential (orography * gravity), simplification compared to
    # Frierson to reduce the number of expensive log calls given that z ≈ Z for most
    # surface temperature variations
    T₀ = atmosphere.reference_temperature
    gravity = planet.gravity
    Δp_geopot_full = geopot.Δp_geopot_full
    Z = T₀ * Δp_geopot_full[nlayers] / gravity
    logZ_z₀ = log(Z / z₀)

    u = get_prognostic_step(vars.grid.u, time_stepping, diffusion)
    v = get_prognostic_step(vars.grid.v, time_stepping, diffusion)
    T = get_prognostic_step(vars.grid.temperature, time_stepping, diffusion)
    geopotential = vars.dynamics.geopotential
    (; orography) = orog

    # Boundary layer depth is highest layer for which Ri < Ri_c (the "critical" threshold)
    # as well as all layers below
    Ri = bulk_richardson!(ij, vars, diffusion, time_stepping, atmosphere, planet, orog, land_sea_mask)
    kₕ::Int = nlayers
    while kₕ > 0 && Ri[ij, kₕ] < Ri_c
        kₕ -= 1
    end
    kₕ += 1  # uppermost layer where Ri < Ri_c

    # for output, TODO as layer index or height?
    vars.parameterizations.boundary_layer_height[ij] = kₕ

    # reuse scratch array for diffusion coefficients
    K = vars.scratch.grid.b

    # diffusion above boundary layer is 0, reset for all layers
    for k in 1:nlayers
        K[ij, k] = 0
    end

    if kₕ <= nlayers    # boundary layer depth is at least 1 layer thick (calculate diffusion)

        # Calculate diffusion coefficients following Frierson 2006, eq. 16-20
        # h always non-negative to avoid error in log
        h = max(geopotential[ij, kₕ] * gravity⁻¹ - orography[ij], 0)
        Ri_N = Ri[ij, nlayers]                      # surface bulk Richardson number
        Ri_N = clamp(Ri_N, 0, Ri_c)                 # cases of eq. 12-14
        sqrtC = (κ / logZ_z₀) * (1 - Ri_N / Ri_c)           # sqrt of eq. 12-14
        surface_speed = sqrt(u[ij, nlayers]^2 + v[ij, nlayers]^2)
        K0 = κ * surface_speed * sqrtC              # height-independent K eq. 19, 20

        for k in kₕ:nlayers
            # height [m] above surface
            z = max(geopotential[ij, k] * gravity⁻¹ - orography[ij], z₀)
            zmin = min(z, fb * h)         # height [m] to evaluate Kb(z) at
            K_k = K0 * zmin             # = κ*u_N*√Cz in eq. (19, 20)

            # multiply with z-dependent factor in eq. (18) ?
            K_k *= z < fb * h ? one(K0) : zfac(z, h, fb)

            # multiply with Ri-dependent factor in eq. (20) for stable surface layer, using
            # the surface bulk Richardson number Ri_N (Ri at kₕ is ≈ Ri_c by construction)
            K_k *= Ri_N <= 0 ? one(K0) : Rifac(Ri_N, Ri_c, logZ_z₀)
            # convert K [m²/s] from z to σ coordinates, ∂z = -gρ/pₛ ∂σ = -gσ/(RT) ∂σ, into [1/s]
            dσdz = gravity * σ_levels_full[k] / (atmosphere.R_dry * T[ij, k])
            K[ij, k] = K_k * dσdz^2     # write diffusion coefficient into array
        end
    end

    # return diffusion coefficients and height index of boundary layer
    return K, kₕ
end

# z-dependent factor in Frierson, 2006 eq (18)
@inline zfac(z, h, fb) = z / (fb * h) * (1 - (z - fb * h) / ((1 - fb) * h))^2

# Ri-dependent factor in Frierson, 2006 eq (20)
@inline function Rifac(Ri, Ri_c, z, z₀)
    Ri_Ri_c = Ri / Ri_c
    return inv(1 + Ri_Ri_c * log(z / z₀) / (1 - Ri_Ri_c))
end

# Approximate: Ri-dependent factor in Frierson, 2006 eq (20)
# because 1 / (1 + log(z/z₀)) is so weakly dependent on z for 10-10000m
@inline function Rifac(Ri, Ri_c, logz_z₀)
    Ri_Ri_c = Ri / Ri_c
    return inv(1 + Ri_Ri_c * logz_z₀ / (1 - Ri_Ri_c))
end

"""$(TYPEDSIGNATURES)
Time step over which the implicit vertical diffusion is solved. Leapfrog evaluates the
parameterizations at the previous step and applies their tendencies over 2Δt."""
implicit_vertical_diffusion_time_step(time_stepping::AbstractLeapfrog) = 2 * time_stepping.Δt
implicit_vertical_diffusion_time_step(time_stepping::AbstractTimeStepper) = time_stepping.Δt

"""
$(TYPEDSIGNATURES)
Vertical diffusion of `var` within the boundary layer (layers `kₕ` to the surface) with
diffusion coefficients `K` [1/s] in σ coordinates, no flux through the top of the boundary
layer and the surface. Implicit (backward Euler) in time over `Δt` for unconditional stability,
solved with the tridiagonal Thomas algorithm, the change is accumulated as
`scale * (var_new - var) / Δt` into `tend`."""
@propagate_inbounds function _vertical_diffusion!(
        ij,         # horizontal grid point ij
        vars,       # for scratch arrays
        tend,       # tendency to accumulate diffusion into
        var,        # variable calculate diffusion from
        K,          # diffusion coefficients [1/s] in σ coordinates
        kₕ,         # uppermost layer that's still within the boundary layer
        diffusion::BulkRichardsonDiffusion,
        Δt,         # time step [s]
        scale = 1,  # scale the tendency, e.g. 1/cₚ to convert dry static energy to temperature
    )
    (; ∇²_above, ∇²_below) = diffusion
    nlayers = size(tend, 2)
    c′ = vars.scratch.grid.vertical_diffusion_c    # modified upper diagonal
    d′ = vars.scratch.grid.vertical_diffusion_d    # modified right-hand side, then solution

    # tridiagonal system -a[k]*x[k-1] + (1 + a[k] + c[k])*x[k] - c[k]*x[k+1] = var[k]
    # with a, c the diffusive exchange (×Δt) with the layer above, below
    # 1/2 of the average diffusion coefficient K is already baked into the ∇² operators
    # forward sweep, no flux through the top of the boundary layer (a = 0 at k = kₕ)
    for k in kₕ:nlayers
        a = k > kₕ ? Δt * ∇²_above[k] * (K[ij, k] + K[ij, k - 1]) : zero(Δt)
        c = k < nlayers ? Δt * ∇²_below[k] * (K[ij, k] + K[ij, k + 1]) : zero(Δt)
        c′ₖ₋₁ = k > kₕ ? c′[ij, k - 1] : zero(Δt)
        d′ₖ₋₁ = k > kₕ ? d′[ij, k - 1] : zero(Δt)
        denominator = 1 + a + c - a * c′ₖ₋₁
        c′[ij, k] = c / denominator
        d′[ij, k] = (var[ij, k] + a * d′ₖ₋₁) / denominator
    end

    # back substitution, d′ becomes the solution var_new
    for k in (nlayers - 1):-1:kₕ
        d′[ij, k] += c′[ij, k] * d′[ij, k + 1]
    end

    for k in kₕ:nlayers
        tend[ij, k] += scale * (d′[ij, k] - var[ij, k]) / Δt
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
Calculate the bulk Richardson number following Frierson, 2007.
For vertical stability in the boundary layer."""
@propagate_inbounds function bulk_richardson!(
        ij,
        vars,
        diffusion::BulkRichardsonDiffusion,
        time_stepping::AbstractTimeStepper,
        atmosphere::AbstractAtmosphere,
        planet::AbstractPlanet,
        orog,
        land_sea_mask,
    )
    # reuse work array
    Ri = vars.scratch.grid.a
    nlayers = size(Ri, 2)
    surface = nlayers       # surface index
    cₚ = atmosphere.heat_capacity

    u = get_prognostic_step(vars.grid.u, time_stepping, diffusion)
    v = get_prognostic_step(vars.grid.v, time_stepping, diffusion)
    Φ = vars.dynamics.geopotential
    T = get_prognostic_step(vars.grid.temperature, time_stepping, diffusion)

    # for dry models, use scratch array to bypass access to non-existing humidity variable
    for k in 1:nlayers
        vars.scratch.grid.b[ij, k] = 0
    end   # reset to zero humidity
    q = haskey(vars.grid, :humidity) ?
        get_prognostic_step(vars.grid.humidity, time_stepping, diffusion) :
        vars.scratch.grid.b

    # surface layer, between the surface (skin temperature, orography) and the lowermost layer
    V² = u[ij, surface]^2 + v[ij, surface]^2
    Tₛ = surface_skin_temperature(ij, vars, land_sea_mask, time_stepping, diffusion, T[ij, surface])
    Tᵥₛ = virtual_temperature(Tₛ, q[ij, surface], atmosphere)
    Φₛ = planet.gravity * orog.orography[ij]                                    # surface geopotential
    Θ₀ = cₚ * Tᵥₛ + Φₛ                                                        # virtual dry static energy at surface
    Θ₁ = cₚ * virtual_temperature(T[ij, surface], q[ij, surface], atmosphere) + Φ[ij, surface]  # and at lowermost layer
    Ri[ij, surface] = (Φ[ij, surface] - Φₛ) * (Θ₁ - Θ₀) / (cₚ * Tᵥₛ * V²)

    for k in 1:(nlayers - 1)
        V² = u[ij, k]^2 + v[ij, k]^2
        Tᵥ = virtual_temperature(T[ij, k], q[ij, k], atmosphere)
        virtual_dry_static_energy = cₚ * Tᵥ + Φ[ij, k]
        Ri[ij, k] = Φ[ij, k] * (virtual_dry_static_energy - Θ₁) / (Θ₁ * V²)
    end

    return Ri
end
