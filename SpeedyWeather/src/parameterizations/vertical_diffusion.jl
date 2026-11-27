abstract type AbstractVerticalDiffusion <: AbstractParameterization end

export BulkRichardsonDiffusion

"""Vertical diffusion of momentum, dry static energy and humidity with a diffusivity K at the layer
interfaces that is the maximum of three parts: a boundary-layer profile following Frierson et al. 2006
with the boundary layer height from a bulk Richardson number, a local Richardson-number-dependent
diffusivity following Louis et al. 1982 and Holtslag and Boville 1993, and a background diffusivity
decaying with height. The diffusion is formulated as a flux divergence on the pressure thickness of
the layers, conserving the mass-weighted column integral, and integrated implicitly so that it is stable
for any diffusivity and time step, by default over-implicitly (`implicitness`) as the diffusivity is taken
from the previous time step. The boundary-layer profile uses the surface drag
and wind of the boundary layer (`model.boundary_layer`), without one only the local and the background
diffusivity act. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct BulkRichardsonDiffusion{NF} <: AbstractVerticalDiffusion
    "[OPTION] von Kármán constant [1]"
    @param von_Karman::NF = 0.4 (bounds = 0 .. 1,)

    "[OPTION] Critical bulk Richardson number defining the boundary layer height [1]"
    @param critical_Richardson::NF = 1 (bounds = Positive,)

    "[OPTION] Surface bulk Richardson number at and above which the boundary-layer profile is zero [1]"
    @param Richardson_cutoff::NF = 1 (bounds = Positive,)

    "[OPTION] Fraction of the boundary layer height that is the surface layer [1]"
    @param surface_layer_fraction::NF = 0.1 (bounds = 0 .. 1,)

    "[OPTION] Roughness length [m] for the stability correction, only used if the model has no surface roughness"
    @param roughness_length::NF = 3.21e-5 (bounds = Positive,)

    "[OPTION] Asymptotic mixing length [m] of the local diffusivity below `mixing_length_height`"
    @param mixing_length::NF = 150 (bounds = Nonnegative,)

    "[OPTION] Asymptotic mixing length [m] of the local diffusivity in the free atmosphere"
    @param free_mixing_length::NF = 30 (bounds = Nonnegative,)

    "[OPTION] Height [m] above which the asymptotic mixing length decays to `free_mixing_length`"
    @param mixing_length_height::NF = 1000 (bounds = Positive,)

    "[OPTION] Minimum asymptotic mixing length [m] of the buoyancy-driven diffusivity in unstable stratification"
    @param convective_mixing_length::NF = 150 (bounds = Nonnegative,)

    "[OPTION] Background diffusivity [m²/s] at the surface"
    @param background_diffusivity::NF = 0.1 (bounds = Nonnegative,)

    "[OPTION] e-folding height [m] of the background diffusivity"
    @param background_height::NF = 1000 (bounds = Positive,)

    "[OPTION] Maximum diffusivity [m²/s]"
    @param max_diffusivity::NF = 1000 (bounds = Positive,)

    "[OPTION] Implicitness factor of the time integration [1]: 1 is backward Euler, larger is over-implicit and damps the oscillations that a large diffusivity from the previous time step can cause (1.5 as in the IFS)"
    @param implicitness::NF = 1.5 (bounds = Positive,)

    "[OPTION] Diffuse dry static energy (i.e. temperature)?"
    diffuse_static_energy::Bool = true

    "[OPTION] Diffuse momentum?"
    diffuse_momentum::Bool = true

    "[OPTION] Diffuse humidity? Ignored for models without humidity"
    diffuse_humidity::Bool = true
end

Adapt.@adapt_structure BulkRichardsonDiffusion

# generator function
BulkRichardsonDiffusion(SG::SpectralGrid; kwargs...) = BulkRichardsonDiffusion{SG.NF}(; kwargs...)

# The conductance a = ρK/Δz [kg/m²/s] of the interface k+½ below layer k is stored, not the diffusivity K [m²/s]:
# the implicit solver only needs the product, which saves the array for ρ/Δz. No flux through the surface
# interface (k = nlayers), where the conductance is zero. K = a Δz/ρ can be recovered from the state
variables(::BulkRichardsonDiffusion) = (
    ParameterizationVariable(:boundary_layer_height, Grid2D(), desc = "Boundary layer height above the surface", units = "m"),
    ParameterizationVariable(:vertical_conductance_momentum, GridXYZ(), desc = "Conductance (density times diffusivity over layer distance) of the vertical diffusion of momentum at the interface below each layer", units = "kg/m^2/s"),
    ParameterizationVariable(:vertical_conductance_heat, GridXYZ(), desc = "Conductance (density times diffusivity over layer distance) of the vertical diffusion of heat and humidity at the interface below each layer", units = "kg/m^2/s"),

    #TODO: reusing the scratch memory for a parameterization has caused problems with Reactant before
    # so far this seems fine but it might need be revisited
    ScratchVariable(:a, GridXYZ(), desc = "Scratch array", namespace = :grid),
    ScratchVariable(:b, GridXYZ(), desc = "Scratch array", namespace = :grid),
    ScratchVariable(:a_2D, Grid2D(), desc = "Scratch array", namespace = :grid),
    ScratchVariable(:b_2D, Grid2D(), desc = "Scratch array", namespace = :grid),
)

initialize!(::BulkRichardsonDiffusion, ::PrimitiveEquation) = nothing

# The scheme is written for efficient performance on both CPU and GPU.
# On CPU the outer loop is over the layers, the inner loop over the columns, so that the inner loop reads
# contiguous memory and vectorises. On GPU the outer loop is over the columns (as are all GPU kernels in
# SpeedyWeather), the inner loop over the layers, so that memory access is coalesced (hopefully).
#
# To enable both the parameterization is written for ranges of columns with the GPU just always passing a range
# of one column.

"""Helper tool for GPU/Reactant. The single column `ij` as a range of columns whose length of one is known at compile
time, so that every loop over it compiles to straight-line code without run-time bounds."""
struct SingleColumn{T <: Integer} <: AbstractUnitRange{T}
    ij::T
end

Base.first(column::SingleColumn) = column.ij
Base.last(column::SingleColumn) = column.ij
Base.length(::SingleColumn) = 1
Base.iterate(column::SingleColumn) = (column.ij, nothing)
Base.iterate(::SingleColumn, ::Nothing) = nothing

# One column of the fused kernel (GPU, Reactant), for which the inner loops have one iteration
@propagate_inbounds parameterization!(ij, vars, diffusion::BulkRichardsonDiffusion, model) =
    vertical_diffusion!(SingleColumn(ij), vars, diffusion, model)

# On CPU all columns at once
@inline function column_parameterization_cpu!(vars, diffusion::BulkRichardsonDiffusion, model)
    @inbounds vertical_diffusion!(1:model.geometry.npoints, vars, diffusion, model)
    return nothing
end

"""$(TYPEDSIGNATURES)
Vertical diffusion tendencies of u, v, temperature (via dry static energy) and humidity for the range
`columns` of horizontal grid points, solved implicitly with the conductances from `vertical_conductance!`.
The variables that share a conductance (u and v, dry static energy and humidity) are solved together.

The conductances are those of the state before the time step. With such explicit coefficients backward
Euler is stable for a given diffusivity, but a diffusivity that grows with the diffused state itself
(K ∝ u⋆ ∝ |V|, the boundary layer height) makes columns with a large diffusivity oscillate between two
states and grow (Kalnay and Kanamitsu, 1988). The solve is therefore over-implicit as in the IFS (IFS
documentation Cy49r1, part IV, section 3.5): with the implicitness factor γ = `diffusion.implicitness`
the tendency y = (x⁺ - x⁻)/Δt is the flux divergence of the state x̂ = γx⁺ + (1 - γ)x⁻ = x⁻ + γΔt y, which
is the implicit solve with the time step γΔt. γ = 1 is backward Euler."""
@propagate_inbounds function vertical_diffusion!(columns::AbstractUnitRange, vars, diffusion::BulkRichardsonDiffusion, model)
    (; diffuse_momentum, diffuse_static_energy, diffuse_humidity) = diffusion

    # escape immediately if all diffusions disabled
    any((diffuse_momentum, diffuse_static_energy, diffuse_humidity)) || return nothing

    vertical_conductance!(columns, vars, diffusion, model)

    # time step of the implicit solve: the implicitness factor times the time step (2Δt for leapfrog), so
    # that the tendency is that of the over-implicit state x⁻ + γΔt y. The model steps with Δt as always
    Δt = diffusion.implicitness * default_time_step(model.time_stepping)
    conductance_momentum = vars.parameterizations.vertical_conductance_momentum
    conductance_heat = vars.parameterizations.vertical_conductance_heat
    solver_arguments = (
        vars.scratch.grid.b,                        # general scratch array for the upper diagonal
        vars.parameterizations.surface_pressure,
        Δt, model.planet.gravity, model.geometry.vertical_coordinates,
    )

    if diffuse_momentum
        u = diffused_wind(vars.grid.u, vars.tendencies.grid.u, diffusion, model)
        v = diffused_wind(vars.grid.v, vars.tendencies.grid.v, diffusion, model)
        implicit_vertical_diffusion!(columns, (u, v), conductance_momentum, solver_arguments...)
    end

    has_humidity = haskey(vars.grid, :humidity)
    if has_humidity && diffuse_humidity && diffuse_static_energy
        static_energy = diffused_static_energy(vars, diffusion, model)
        humidity = diffused_humidity(vars, diffusion, model)
        implicit_vertical_diffusion!(columns, (static_energy, humidity), conductance_heat, solver_arguments...)
    elseif diffuse_static_energy
        static_energy = diffused_static_energy(vars, diffusion, model)
        implicit_vertical_diffusion!(columns, (static_energy,), conductance_heat, solver_arguments...)
    elseif has_humidity && diffuse_humidity
        humidity = diffused_humidity(vars, diffusion, model)
        implicit_vertical_diffusion!(columns, (humidity,), conductance_heat, solver_arguments...)
    end
    return nothing
end

# A variable to diffuse is described by (tendency, variable, geopotential_term), see implicit_vertical_diffusion!
@propagate_inbounds function diffused_wind(wind, wind_tendency, diffusion, model)
    (; time_stepping) = model
    return (get_tendency_step(wind_tendency, time_stepping, diffusion), get_prognostic_step(wind, time_stepping, diffusion), nothing)
end

# diffuse the dry static energy s = cₚT + Φ (divided by cₚ) with Φ held fixed, so that the
# tendency of s/cₚ is the temperature tendency
@propagate_inbounds function diffused_static_energy(vars, diffusion, model)
    (; time_stepping) = model
    T = get_prognostic_step(vars.grid.temperature, time_stepping, diffusion)
    T_tend = get_tendency_step(vars.tendencies.grid.temperature, time_stepping, diffusion)
    return (T_tend, T, (vars.dynamics.geopotential, inv(model.atmosphere.heat_capacity)))
end

@propagate_inbounds function diffused_humidity(vars, diffusion, model)
    (; time_stepping) = model
    q = get_prognostic_step(vars.grid.humidity, time_stepping, diffusion)
    q_tend = get_tendency_step(vars.tendencies.grid.humidity, time_stepping, diffusion)
    return (q_tend, q, nothing)
end

# humidity of layer k, zero for models without humidity
@propagate_inbounds column_humidity(q, ij, k) = q[ij, k]
@inline column_humidity(::Nothing, ij, k) = false     # strong zero

# roughness length for the stability correction, from the surface roughness scheme if present
@propagate_inbounds column_roughness_length(roughness, default, ij) = ifelse(roughness[ij] > 0, roughness[ij], default)
@inline column_roughness_length(::Nothing, default, ij) = default

"""$(TYPEDSIGNATURES)
Friction velocity u⋆ = √C_M Vₛ [m/s] at grid point `ij` from the momentum drag coefficient `drag` and the
`surface_wind_speed` of the boundary layer. Zero without them (`nothing`), which switches the
boundary-layer profile of the diffusivity off."""
@propagate_inbounds friction_velocity(drag, surface_wind_speed, ij) = sqrt(abs(drag[ij])) * surface_wind_speed[ij]
@inline friction_velocity(::Nothing, ::Nothing, ij) = false     # strong zero

# drag and surface wind speed of a boundary layer, or nothing for both if `parameterizations` lacks one
@inline function boundary_layer_drag_and_wind(parameterizations)
    has_boundary_layer = haskey(parameterizations, :boundary_layer_drag_momentum) &&
        haskey(parameterizations, :surface_wind_speed)
    has_boundary_layer || return nothing, nothing
    return parameterizations.boundary_layer_drag_momentum, parameterizations.surface_wind_speed
end

"""$(TYPEDSIGNATURES)
Conductances a = ρK/Δz [kg/m²/s] of the vertical diffusion of momentum and heat at the interfaces k+½
below every layer k for the range `columns` of horizontal grid points (typically a single column on GPU, or 
all columns on CPU), with the diffusivity K [m²/s], the air density ρ at the interface and the distance Δz 
between the two layers, written into `vars.parameterizations.vertical_conductance_momentum` and `_heat`, and 
the boundary layer height into `vars.parameterizations.boundary_layer_height`. The diffusivity is
K = min(max(K_profile, K_local, K_background), K_max) with

- K_profile the boundary-layer profile of Frierson et al. 2006, eqs. 16–20, from the friction velocity
  of the surface drag, a stability correction from the bulk Richardson number between surface and
  lowermost layer that vanishes at `Richardson_cutoff`, and the boundary layer height h where the bulk Richardson number relative to the
  lowermost layer reaches `critical_Richardson`, interpolated linearly between layers;
- K_local = l² |∂V/∂z| f(Ri) with the gradient Richardson number Ri between adjacent layers, the
  mixing length 1/l = 1/(κz) + 1/λ, long-tail stability functions of Louis et al. 1982 for Ri > 0 and
  f = √(1 - 18Ri) for Ri < 0 (Holtslag and Boville, 1993) with the buoyancy-driven part using a mixing
  length of at least `convective_mixing_length`;
- K_background = `background_diffusivity` ⋅ exp(-z/`background_height`).

Computed in five loops over the layers with the columns inside: virtual temperature, boundary layer
height, the stability of the surface layer (per column), the stability correction of every interface 
and the conductances."""
@propagate_inbounds function vertical_conductance!(columns::AbstractUnitRange, vars, diffusion::BulkRichardsonDiffusion, model)
    (; time_stepping, atmosphere, geometry) = model
    (; nlayers, vertical_coordinates) = geometry
    g = model.planet.gravity
    inv_g = inv(g)
    cₚ = atmosphere.heat_capacity
    inv_cₚ = inv(cₚ)
    inv_R_dry = inv(atmosphere.R_dry)
    κ = diffusion.von_Karman
    Ri_c = diffusion.critical_Richardson
    Ri_cutoff = diffusion.Richardson_cutoff
    fb = diffusion.surface_layer_fraction
    (; mixing_length, free_mixing_length, mixing_length_height, convective_mixing_length) = diffusion
    (; background_diffusivity, background_height, max_diffusivity) = diffusion
    inv_mixing_length_height = inv(mixing_length_height)
    inv_background_height = inv(background_height)

    u = get_prognostic_step(vars.grid.u, time_stepping, diffusion)
    v = get_prognostic_step(vars.grid.v, time_stepping, diffusion)
    T = get_prognostic_step(vars.grid.temperature, time_stepping, diffusion)
    q = haskey(vars.grid, :humidity) ? get_prognostic_step(vars.grid.humidity, time_stepping, diffusion) : nothing
    Φ = vars.dynamics.geopotential      # includes the orography
    (; orography) = model.orography     # surface geopotential Φₛ = g * orography
    (; surface_pressure, boundary_layer_height) = vars.parameterizations
    drag, surface_wind_speed = boundary_layer_drag_and_wind(vars.parameterizations)

    # roughness length for the stability correction, from the surface roughness scheme if present
    roughness = haskey(vars.parameterizations, :momentum_roughness) ? vars.parameterizations.momentum_roughness : nothing
    default_roughness = diffusion.roughness_length

    conductance_momentum = vars.parameterizations.vertical_conductance_momentum
    conductance_heat = vars.parameterizations.vertical_conductance_heat
    NF = eltype(conductance_momentum)

    # 1. VIRTUAL TEMPERATURE of every layer, in a general scratch array that the implicit solver reuses
    # afterwards. The loops below read it instead of temperature and humidity
    Tᵥ = vars.scratch.grid.a
    for k in 1:nlayers
        @simd ivdep for ij in columns
            Tᵥ[ij, k] = virtual_temperature(T[ij, k], column_humidity(q, ij, k), atmosphere)
        end
    end

    # 2. BOUNDARY LAYER HEIGHT h [m] above the surface
    boundary_layer_height!(columns, boundary_layer_height, u, v, Tᵥ, Φ, orography, g, cₚ, Ri_c)

    # 3. SURFACE LAYER, per column. Stable correction of Frierson 2006, eq. 20,
    # [1 + Ri/Ri_cutoff log(z/z₀)/(1 - Ri/Ri_cutoff)]⁻¹ with the surface Richardson number, here written as
    # (1 - r)/(1 - r + r log(z/z₀)), r = Ri/Ri_cutoff ∈ [0, 1], to avoid a division by zero at Ri = Ri_cutoff.
    # It is 1 for unstable conditions (Ri ≤ 0) and 0 for Ri ≥ Ri_cutoff. Frierson uses the critical Richardson
    # number here; a separate cutoff (as for the drag) keeps boundary-layer mixing in stable conditions until
    # buoyancy dominates the shear. Stored are r and the correction at the surface layer height z = fb*h,
    # which holds for all interfaces above the surface layer
    stability_ratio = vars.scratch.grid.a_2D
    stability_surface_layer = vars.scratch.grid.b_2D
    for ij in columns
        ΔΦ₀ = Φ[ij, nlayers] - g * orography[ij]       # geopotential of the lowermost layer relative to the surface
        r = clamp(surface_richardson(ij, ΔΦ₀, vars, diffusion, model) / Ri_cutoff, 0, 1)
        z₀ = column_roughness_length(roughness, default_roughness, ij)
        stability_ratio[ij] = r
        stability_surface_layer[ij] = stable_correction(r, fb * boundary_layer_height[ij], z₀)
    end

    # 4. STABILITY CORRECTION of the boundary-layer profile at every interface k+½, with its `log` only for the
    # interfaces inside the surface layer (few). That branch does not vectorise, so it is evaluated here and
    # not in the loop of the conductances, which then does. The conductance of heat holds the correction
    # until that loop has read it
    profile_stability = conductance_heat
    for k in 1:(nlayers - 1)
        for ij in columns
            z = ((Φ[ij, k] + Φ[ij, k + 1]) / 2 - g * orography[ij]) * inv_g     # height of the interface above the surface [m]

            stability = stability_surface_layer[ij]
            if z < fb * boundary_layer_height[ij]       # inside the surface layer
                stability = stable_correction(stability_ratio[ij], z, column_roughness_length(roughness, default_roughness, ij))
            end
            profile_stability[ij, k] = stability
        end
    end

    # 5. CONDUCTANCES of the interfaces k+½ between layer k (above) and k+1 (below). Multiplications with
    # reciprocals where the same divisor occurs several times, and sqrt(abs(x)) for the square roots: all
    # their arguments are positive, but `sqrt` alone checks that and throws, which prevents vectorisation.
    # `@inline exp`: Base declares `exp` `@noinline`, and a call of it prevents the vectorisation of the loop
    # too, inlined at the call site it vectorises
    for k in 1:(nlayers - 1)
        @simd ivdep for ij in columns
            Tᵥ_above = Tᵥ[ij, k]
            Tᵥ_below = Tᵥ[ij, k + 1]
            Φ_above = Φ[ij, k]
            Φ_below = Φ[ij, k + 1]
            ΔΦ = Φ_above - Φ_below                                      # geopotential difference between layers
            inv_Δz = g / ΔΦ                                             # 1 / distance between layers [1/m]
            z = ((Φ_above + Φ_below) / 2 - g * orography[ij]) * inv_g   # height of the interface above the surface [m]
            inv_Tᵥ = 2 / (Tᵥ_above + Tᵥ_below)                          # 1 / virtual temperature at the interface

            # air density at the interface over the distance between layers [kg/m⁴]
            p = pressure_half(k + 1, surface_pressure[ij], vertical_coordinates)
            density_over_distance = p * inv_R_dry * inv_Tᵥ * inv_Δz

            # 1. BOUNDARY LAYER PROFILE, Frierson 2006, eqs. 16-20, Kb(z) = κ u⋆ z in the surface layer
            # z < fb*h, cubic decay to zero at z = h above, zero above the boundary layer
            h = boundary_layer_height[ij]
            surface_layer_height = fb * h
            above_surface_layer = (z > surface_layer_height) & (z < h)
            ξ = ifelse(above_surface_layer, z - surface_layer_height, zero(z)) /
                ifelse(above_surface_layer, (1 - fb) * h, one(z))
            K_profile = κ * friction_velocity(drag, surface_wind_speed, ij) * z * profile_stability[ij, k] * (1 - ξ)^2
            K_profile = ifelse(z < h, K_profile, zero(K_profile))

            # 2. LOCAL DIFFUSIVITY, l² S f(Ri) with S = |∂V/∂z|, squared buoyancy frequency
            # N² = g/(cₚTᵥ) ∂sᵥ/∂z from the virtual dry static energy sᵥ = cₚTᵥ + Φ
            ΔV² = max((u[ij, k] - u[ij, k + 1])^2 + (v[ij, k] - v[ij, k + 1])^2, NF(MINIMUM_WIND_SHEAR²))
            S² = ΔV² * inv_Δz^2
            S = sqrt(abs(S²))
            N² = g * (Tᵥ_above - Tᵥ_below + ΔΦ * inv_cₚ) * inv_Tᵥ * inv_Δz
            unstable = N² < 0

            # gradient Richardson number, stable part, and N² of the unstable part. With ifelse and not as
            # max(N²/S², 0) and min(N², 0): from a `max` LLVM (Julia 1.12) knows the sign of the argument of the
            # square roots in the stability functions, removes their `abs` and keeps the check of `sqrt`, which
            # costs the vectorisation of this loop
            Ri = ifelse(unstable, zero(N²), N² / S²)
            N²_unstable = ifelse(unstable, N², zero(N²))

            # asymptotic mixing length λ, constant below mixing_length_height, then decaying towards
            # free_mixing_length (Holtslag and Boville, 1993, eq. 3.6), mixing length 1/l = 1/(κz) + 1/λ
            λ = free_mixing_length + (mixing_length - free_mixing_length) * @inline(exp(min(1 - z * inv_mixing_length_height, 0)))
            l = κ * z * λ / (κ * z + λ)

            # unstable: l²S√(1 - 18Ri) = l²√(S² - 18N²) (Holtslag and Boville, 1993, eq. 3.7), split into a shear
            # part l²S and a buoyancy part l_c²(√(S² - 18N²) - S) whose mixing length l_c uses at least the
            # convective_mixing_length, so that superadiabatic layers are mixed within a few time steps also in
            # the free atmosphere where λ is small. Continuous at N² = 0, l_c = l reproduces eq. 3.7
            λ_convective = max(λ, convective_mixing_length)
            l_convective = κ * z * λ_convective / (κ * z + λ_convective)
            K_unstable = l^2 * S + l_convective^2 * (sqrt(abs(S² - 18 * N²_unstable)) - S)
            K_local_momentum = ifelse(unstable, K_unstable, l^2 * S * stability_function_momentum(Ri))
            K_local_heat = ifelse(unstable, K_unstable, l^2 * S * stability_function_heat(Ri))

            # 3. BACKGROUND DIFFUSIVITY
            K_background = background_diffusivity * @inline(exp(-z * inv_background_height))

            K_momentum = min(max(K_profile, K_local_momentum, K_background), max_diffusivity)
            K_heat = min(max(K_profile, K_local_heat, K_background), max_diffusivity)

            # store the conductances ρK/Δz [kg/m²/s], all the implicit solver needs
            conductance_momentum[ij, k] = density_over_distance * K_momentum
            conductance_heat[ij, k] = density_over_distance * K_heat
        end
    end

    # no diffusive flux through the surface, which is handled by the surface fluxes
    @simd ivdep for ij in columns
        conductance_momentum[ij, nlayers] = 0
        conductance_heat[ij, nlayers] = 0
    end
    return nothing
end

"""$(TYPEDSIGNATURES)
Bulk Richardson number between the surface and the lowermost layer at grid point `ij`, the land and ocean
values of `bulk_richardson_surface` weighted with the land fraction, with the geopotential `ΔΦ₀` of the
lowermost layer relative to the surface. Read from `vars.parameterizations.surface_bulk_richardson` if the
drag of the boundary layer has stored it (`BulkRichardsonDrag`, which computes it for the same state),
otherwise computed here. Zero for a model without the drag or the surface wind of a boundary layer."""
@propagate_inbounds function surface_richardson(ij, ΔΦ₀, vars, diffusion::BulkRichardsonDiffusion, model)
    (; parameterizations) = vars
    haskey(parameterizations, :surface_bulk_richardson) && return parameterizations.surface_bulk_richardson[ij]

    drag, _ = boundary_layer_drag_and_wind(parameterizations)
    isnothing(drag) && return zero(eltype(parameterizations.vertical_conductance_momentum))

    Ri_land, Ri_ocean = bulk_richardson_surface(ij, ΔΦ₀, vars, model.atmosphere, diffusion, model.time_stepping)
    land_fraction = model.land_sea_mask.land_fraction[ij]
    return land_fraction * Ri_land + (1 - land_fraction) * Ri_ocean
end

# minimum wind speed [m/s] squared for the bulk Richardson number (Holtslag and Boville, 1993)
const MINIMUM_WIND_SPEED² = 1

# minimum wind speed difference between layers [m/s] squared for the shear in the local diffusivity
const MINIMUM_WIND_SHEAR² = 0.01

"""$(TYPEDSIGNATURES)
Stability correction [1 + Ri/Ri_cutoff log(z/z₀)/(1 - Ri/Ri_cutoff)]⁻¹ of Frierson 2006, eq. 20, written
with r = Ri/Ri_cutoff ∈ [0, 1] as (1 - r)/(1 - r + r log(z/z₀)) to be finite at r = 1. z is capped one order of
magnitude above the roughness length `z₀` as in `neutral_transfer_coefficient`."""
@inline stable_correction(r, z, z₀) = (1 - r) / (1 - r + r * log(max(z, 10z₀) / z₀))

"""$(TYPEDSIGNATURES)
Long-tail stability function for momentum of Louis et al. 1982, 1/(1 + 2bRi/√(1 + dRi)) with b = d = 5,
for the gradient Richardson number `Ri` ≥ 0. `sqrt(abs(x))` as the argument is positive but `sqrt` alone
checks that and throws, which prevents the vectorisation of the loop that calls it."""
@inline stability_function_momentum(Ri) = 1 / (1 + 10Ri / sqrt(abs(1 + 5Ri)))

"""$(TYPEDSIGNATURES)
Long-tail stability function for heat of Louis et al. 1982, 1/(1 + 3bRi√(1 + dRi)) with b = d = 5,
for the gradient Richardson number `Ri` ≥ 0, see `stability_function_momentum` for `sqrt(abs(x))`."""
@inline stability_function_heat(Ri) = 1 / (1 + 15Ri * sqrt(abs(1 + 5Ri)))

"""$(TYPEDSIGNATURES)
Boundary layer height h [m] above the surface for the range `columns` of horizontal grid points, written
into `height`: the height where the bulk Richardson number relative to the lowermost layer N,
Ri(z) = gz (sᵥ(z) - sᵥ_N) / (cₚTᵥ_N |V(z)|²) with the virtual dry static energy sᵥ = cₚTᵥ + Φ (Frierson et
al. 2006) from the virtual temperature `Tᵥ` of the layers, first reaches `Ri_c` from below, interpolated
linearly between the layers that bracket Ri_c with Ri = 0 at the lowermost layer. The height of the
uppermost layer if Ri < Ri_c everywhere. Layers from the surface upwards in the outer loop, with a fixed
trip count and no early exit: a negative height marks the columns without a crossing so far."""
@propagate_inbounds function boundary_layer_height!(columns::AbstractUnitRange, height, u, v, Tᵥ, Φ, orography, g, cₚ, Ri_c)
    nlayers = size(Tᵥ, 2)
    inv_g = inv(g)

    # top layer height if no crossing of Ri_c, negative until a crossing is found
    @simd ivdep for ij in columns
        height[ij] = -(Φ[ij, 1] - g * orography[ij]) * inv_g
    end

    for k in (nlayers - 1):-1:1
        @simd ivdep for ij in columns
            Φₛ = g * orography[ij]                              # surface geopotential
            Tᵥ_N = Tᵥ[ij, nlayers]
            sᵥ_N = cₚ * Tᵥ_N + Φ[ij, nlayers]                   # virtual dry static energy of the lowermost layer

            # Ri = numerator/denominator of layer k and of the layer below, where Ri = 0 for the lowermost
            # layer. With gz = Φ - Φₛ
            ΔΦ = Φ[ij, k] - Φₛ
            ΔΦ_below = Φ[ij, k + 1] - Φₛ
            numerator = ΔΦ * (cₚ * Tᵥ[ij, k] + Φ[ij, k] - sᵥ_N)
            numerator_below = ΔΦ_below * (cₚ * Tᵥ[ij, k + 1] + Φ[ij, k + 1] - sᵥ_N)
            denominator = cₚ * Tᵥ_N * max(u[ij, k]^2 + v[ij, k]^2, MINIMUM_WIND_SPEED²)
            denominator_below = cₚ * Tᵥ_N * max(u[ij, k + 1]^2 + v[ij, k + 1]^2, MINIMUM_WIND_SPEED²)

            # first crossing of Ri_c from below, Ri_below < Ri_c ≤ Ri, found without a division. The weight
            # (Ri_c - Ri_below)/(Ri - Ri_below) of the interpolation with one, guarded for the other layers
            h = height[ij]
            crossing = (h < 0) & (numerator >= Ri_c * denominator)
            weight = ifelse(crossing, (Ri_c * denominator_below - numerator_below) * denominator, zero(h)) /
                ifelse(crossing, numerator * denominator_below - numerator_below * denominator, one(h))
            height[ij] = ifelse(crossing, (ΔΦ_below + weight * (ΔΦ - ΔΦ_below)) * inv_g, h)
        end
    end

    @simd ivdep for ij in columns
        height[ij] = abs(height[ij])
    end
    return nothing
end

# value of the variable in layer k that is diffused, with or without a geopotential term
@propagate_inbounds diffused_value(variable, ::Nothing, ij, k) = variable[ij, k]
@propagate_inbounds diffused_value(variable, (Φ, scale)::Tuple, ij, k) = variable[ij, k] + scale * Φ[ij, k]

# TODO: we may also factor out the triadiagonal solve into one of the subpackages, as it could beused elsewhere
# this is fairly generic code.

"""$(TYPEDSIGNATURES)
Implicit (backward Euler) vertical diffusion for the range `columns` of horizontal grid points of every
variable in `fields`, a tuple with one `(tendency, variable, geopotential_term)` per variable, with the
`conductance` a = ρK/Δz [kg/m²/s] (density times diffusivity over the distance between layers) on the
interfaces below every layer that these variables share. Diffused is `variable` (plus `scale * Φ` with
`geopotential_term = (Φ, scale)`, or `nothing`). With the diffusion operator L as a flux divergence on the
pressure thickness Δp of the layers,

    (L x)ₖ = g/Δpₖ [aₖ₊½ (xₖ₊₁ - xₖ) - aₖ₋½ (xₖ - xₖ₋₁)],     aₖ₊½ = ρₖ₊½ Kₖ₊½ / Δzₖ₊½,

and no flux through the top and the surface, the tendency y = (x⁺ - x⁻)/Δt is the solution of the
tridiagonal system (1 - ΔtL) y = L x⁻ + S, solved with the Thomas algorithm. The source S is the `tendency`
accumulated so far (the surface fluxes into the lowermost layer, the parameterizations called before), so
these tendencies are applied first and the result is diffused, which spreads the surface fluxes over the
column within the time step. `tendency` is replaced by y. Conserves the mass-weighted column integral
Σ xₖ Δpₖ exactly, with that of the source.

`Δt` is the time step of the implicit solve. With the time step of the model this is backward Euler. With
γ times that time step (γ > 1, over-implicit) y is the tendency L x̂ + S of the state x̂ = x⁻ + γΔt y, and a
model step x⁺ = x⁻ + Δt y with it lies between x⁻ and x̂: x⁺ = (1 - 1/γ) x⁻ + x̂/γ.

The matrix is the same for all variables in `fields`: it is eliminated once, and the right-hand side of every
variable is stored in its tendency. Layers in the outer, columns in the inner loops: what the elimination of
layer k needs of the layer above (upper diagonal, right-hand side) is read from the arrays."""
@propagate_inbounds function implicit_vertical_diffusion!(
        columns::AbstractUnitRange,
        fields::Tuple,              # (tendency, variable, geopotential_term) of every variable to diffuse
        conductance,                # ρK/Δz [kg/m²/s] at the interface below every layer
        tridiagonal_upper,          # scratch, upper diagonal after forward elimination
        surface_pressure,           # [Pa] of every column
        Δt,                         # time step [s] of the implicit solve
        g,                          # gravity [m/s²]
        vertical_coordinates,
    )
    nlayers = size(conductance, 2)
    NF = eltype(conductance)

    # FORWARD ELIMINATION, from top to surface
    for k in 1:nlayers
        k_above = max(k - 1, 1)             # clamped, the terms of a layer above the top or below the
        k_below = min(k + 1, nlayers)       # surface are multiplied with a zero conductance
        @simd ivdep for ij in columns
            # conductances [kg/m²/s] at the interfaces above and below layer k, no flux through top and surface
            a_above = ifelse(k > 1, conductance[ij, k_above], zero(NF))
            a_below = ifelse(k < nlayers, conductance[ij, k], zero(NF))
            g_Δp = g / pressure_thickness(k, surface_pressure[ij], vertical_coordinates)

            # row k of (1 - ΔtL): lower, diagonal and upper entry, diagonally dominant
            α = Δt * g_Δp
            lower = -α * a_above
            upper = -α * a_below
            diagonal = 1 + α * (a_above + a_below)

            upper_above = ifelse(k > 1, tridiagonal_upper[ij, k_above], zero(NF))
            inv_pivot = inv(diagonal - lower * upper_above)     # pivot ≥ 1
            tridiagonal_upper[ij, k] = upper * inv_pivot
            forward_elimination_all!(ij, k, k_above, k_below, fields, a_above, a_below, g_Δp, lower, inv_pivot)
        end
    end

    # BACK SUBSTITUTION, from surface to top
    for k in nlayers:-1:1
        k_below = min(k + 1, nlayers)
        @simd ivdep for ij in columns
            back_substitution_all!(ij, k, k_below, fields, tridiagonal_upper[ij, k])
        end
    end
    return nothing
end

# forward elimination of one variable in layer k: the right-hand side, with the tendency accumulated so far as
# the source, replaces the tendency
@propagate_inbounds function forward_elimination!(
        ij, k, k_above, k_below, (tendency, variable, geopotential_term)::Tuple, a_above, a_below, g_Δp, lower, inv_pivot,
    )
    x = diffused_value(variable, geopotential_term, ij, k)
    x_above = diffused_value(variable, geopotential_term, ij, k_above)
    x_below = diffused_value(variable, geopotential_term, ij, k_below)
    explicit_tendency = g_Δp * (a_below * (x_below - x) - a_above * (x - x_above))
    rhs_above = ifelse(k > 1, tendency[ij, k_above], zero(x))
    tendency[ij, k] = (explicit_tendency + tendency[ij, k] - lower * rhs_above) * inv_pivot
    return nothing
end

# back substitution of one variable in layer k: the solution replaces the right-hand side in the tendency
@propagate_inbounds function back_substitution!(ij, k, k_below, (tendency, variable, geopotential_term)::Tuple, upper)
    y_below = ifelse(k < k_below, tendency[ij, k_below], zero(upper))
    tendency[ij, k] -= upper * y_below
    return nothing
end

# the same for every variable of a tuple, by recursion (no closures, so that it inlines and @inbounds propagates)
@inline forward_elimination_all!(ij, k, k_above, k_below, ::Tuple{}, arguments...) = nothing
@propagate_inbounds function forward_elimination_all!(ij, k, k_above, k_below, fields::Tuple, arguments...)
    forward_elimination!(ij, k, k_above, k_below, first(fields), arguments...)
    forward_elimination_all!(ij, k, k_above, k_below, Base.tail(fields), arguments...)
    return nothing
end

@inline back_substitution_all!(ij, k, k_below, ::Tuple{}, arguments...) = nothing
@propagate_inbounds function back_substitution_all!(ij, k, k_below, fields::Tuple, arguments...)
    back_substitution!(ij, k, k_below, first(fields), arguments...)
    back_substitution_all!(ij, k, k_below, Base.tail(fields), arguments...)
    return nothing
end
