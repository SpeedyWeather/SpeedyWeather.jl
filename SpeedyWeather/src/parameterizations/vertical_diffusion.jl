abstract type AbstractVerticalDiffusion <: AbstractParameterization end

export BulkRichardsonDiffusion

"""Vertical diffusion of momentum, dry static energy and humidity with a diffusivity K at the layer
interfaces that is the maximum of three parts: a boundary-layer profile following Frierson et al. 2006
with the boundary layer height from a bulk Richardson number, a local Richardson-number-dependent
diffusivity following Louis et al. 1982 and Holtslag and Boville 1993, and a background diffusivity
decaying with height. The diffusion is formulated as a flux divergence on the pressure thickness of
the layers, conserving the mass-weighted column integral, and integrated implicitly (backward Euler)
so that it is stable for any diffusivity and time step. Requires a boundary layer (`model.boundary_layer`)
for the surface drag and wind. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct BulkRichardsonDiffusion{NF} <: AbstractVerticalDiffusion
    "[OPTION] von Kármán constant [1]"
    @param von_Karman::NF = 0.4 (bounds = 0 .. 1,)

    "[OPTION] Critical bulk Richardson number defining the boundary layer height [1]"
    @param critical_Richardson::NF = 1 (bounds = Positive,)

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

    "[OPTION] Diffuse dry static energy (i.e. temperature)?"
    diffuse_static_energy::Bool = true

    "[OPTION] Diffuse momentum?"
    diffuse_momentum::Bool = true

    "[OPTION] Diffuse humidity? Ignored for models without humidity"
    diffuse_humidity::Bool = true

    "[OPTION] Spread the tendencies of the parameterizations called before (the surface fluxes) implicitly over the column?"
    couple_surface_fluxes::Bool = true
end

Adapt.@adapt_structure BulkRichardsonDiffusion

# generator function
BulkRichardsonDiffusion(SG::SpectralGrid; kwargs...) = BulkRichardsonDiffusion{SG.NF}(; kwargs...)

# diffusivities are defined on the interface k+½ below layer k, the surface interface (k = nlayers) is zero
variables(::BulkRichardsonDiffusion) = (
    ParameterizationVariable(:boundary_layer_height, Grid2D(), desc = "Boundary layer height above the surface", units = "m"),
    ParameterizationVariable(:vertical_diffusivity_momentum, GridXYZ(), desc = "Vertical diffusivity of momentum at the interface below each layer", units = "m^2/s"),
    ParameterizationVariable(:vertical_diffusivity_heat, GridXYZ(), desc = "Vertical diffusivity of heat and humidity at the interface below each layer", units = "m^2/s"),
    ScratchVariable(:density_over_distance, GridXYZ(), desc = "Air density over the distance between layers at the interface below each layer", units = "kg/m^4", namespace = :grid),
    ScratchVariable(:tridiagonal_upper, GridXYZ(), desc = "Upper diagonal of the forward-eliminated tridiagonal system", units = "1", namespace = :grid),
    ScratchVariable(:tridiagonal_rhs, GridXYZ(), desc = "Right-hand side of the forward-eliminated tridiagonal system", units = "?", namespace = :grid),
)

function initialize!(::BulkRichardsonDiffusion, model::PrimitiveEquation)
    isnothing(model.boundary_layer) && throw(
        ArgumentError(
            "BulkRichardsonDiffusion uses the surface drag and wind of the boundary layer, " *
                "model.boundary_layer must not be nothing."
        )
    )
    return nothing
end

# function barrier
@propagate_inbounds parameterization!(ij, vars, diffusion::BulkRichardsonDiffusion, model) =
    vertical_diffusion!(ij, vars, diffusion, model)

"""$(TYPEDSIGNATURES)
Vertical diffusion tendencies of u, v, temperature (via dry static energy) and humidity at grid point
`ij`, solved implicitly with the diffusivities from `vertical_diffusivity!`."""
@propagate_inbounds function vertical_diffusion!(ij, vars, diffusion::BulkRichardsonDiffusion, model)
    (; diffuse_momentum, diffuse_static_energy, diffuse_humidity, couple_surface_fluxes) = diffusion

    # escape immediately if all diffusions disabled
    any((diffuse_momentum, diffuse_static_energy, diffuse_humidity)) || return nothing

    vertical_diffusivity!(ij, vars, diffusion, model)

    (; time_stepping) = model
    Δt = default_time_step(time_stepping)       # 2Δt for leapfrog
    g = model.planet.gravity
    (; vertical_coordinates) = model.geometry
    pₛ = vars.parameterizations.surface_pressure[ij]
    K_momentum = vars.parameterizations.vertical_diffusivity_momentum
    K_heat = vars.parameterizations.vertical_diffusivity_heat
    (; density_over_distance, tridiagonal_upper, tridiagonal_rhs) = vars.scratch.grid
    solver_arguments = (density_over_distance, tridiagonal_upper, tridiagonal_rhs, pₛ, Δt, g, vertical_coordinates, couple_surface_fluxes)

    if diffuse_momentum
        u = get_prognostic_step(vars.grid.u, time_stepping, diffusion)
        v = get_prognostic_step(vars.grid.v, time_stepping, diffusion)
        u_tend = get_tendency_step(vars.tendencies.grid.u, time_stepping, diffusion)
        v_tend = get_tendency_step(vars.tendencies.grid.v, time_stepping, diffusion)
        implicit_vertical_diffusion!(ij, u_tend, u, nothing, K_momentum, solver_arguments...)
        implicit_vertical_diffusion!(ij, v_tend, v, nothing, K_momentum, solver_arguments...)
    end

    if diffuse_static_energy
        # diffuse the dry static energy s = cₚT + Φ (divided by cₚ) with Φ held fixed, so that the
        # tendency of s/cₚ is the temperature tendency
        T = get_prognostic_step(vars.grid.temperature, time_stepping, diffusion)
        T_tend = get_tendency_step(vars.tendencies.grid.temperature, time_stepping, diffusion)
        geopotential_term = (vars.dynamics.geopotential, inv(model.atmosphere.heat_capacity))
        implicit_vertical_diffusion!(ij, T_tend, T, geopotential_term, K_heat, solver_arguments...)
    end

    if diffuse_humidity && haskey(vars.grid, :humidity)
        q = get_prognostic_step(vars.grid.humidity, time_stepping, diffusion)
        q_tend = get_tendency_step(vars.tendencies.grid.humidity, time_stepping, diffusion)
        implicit_vertical_diffusion!(ij, q_tend, q, nothing, K_heat, solver_arguments...)
    end
    return nothing
end

# humidity of layer k, zero for models without humidity
@inline column_humidity(q, ij, k) = q[ij, k]
@inline column_humidity(::Nothing, ij, k) = false     # strong zero

"""$(TYPEDSIGNATURES)
Vertical diffusivities [m²/s] of momentum and heat at the interfaces k+½ below every layer k of column
`ij`, written into `vars.parameterizations.vertical_diffusivity_momentum` and `_heat`, the boundary layer
height into `vars.parameterizations.boundary_layer_height` and the air density over the distance between
the layers at the interfaces into `vars.scratch.grid.density_over_distance`. The diffusivity is
K = min(max(K_profile, K_local, K_background), K_max) with

- K_profile the boundary-layer profile of Frierson et al. 2006, eqs. 16–20, from the friction velocity
  of the surface drag, a stability correction from the bulk Richardson number between surface and
  lowermost layer, and the boundary layer height h where the bulk Richardson number relative to the
  lowermost layer reaches `critical_Richardson`, interpolated linearly between layers;
- K_local = l² |∂V/∂z| f(Ri) with the gradient Richardson number Ri between adjacent layers, the
  mixing length 1/l = 1/(κz) + 1/λ, long-tail stability functions of Louis et al. 1982 for Ri > 0 and
  f = √(1 - 18Ri) for Ri < 0 (Holtslag and Boville, 1993) with the buoyancy-driven part using a mixing
  length of at least `convective_mixing_length`;
- K_background = `background_diffusivity` ⋅ exp(-z/`background_height`)."""
@propagate_inbounds function vertical_diffusivity!(ij, vars, diffusion::BulkRichardsonDiffusion, model)
    (; time_stepping, atmosphere, geometry) = model
    (; nlayers, vertical_coordinates) = geometry
    g = model.planet.gravity
    cₚ = atmosphere.heat_capacity
    R_dry = atmosphere.R_dry
    κ = diffusion.von_Karman
    Ri_c = diffusion.critical_Richardson
    fb = diffusion.surface_layer_fraction
    (; mixing_length, free_mixing_length, mixing_length_height, convective_mixing_length) = diffusion
    (; background_diffusivity, background_height, max_diffusivity) = diffusion

    u = get_prognostic_step(vars.grid.u, time_stepping, diffusion)
    v = get_prognostic_step(vars.grid.v, time_stepping, diffusion)
    T = get_prognostic_step(vars.grid.temperature, time_stepping, diffusion)
    q = haskey(vars.grid, :humidity) ? get_prognostic_step(vars.grid.humidity, time_stepping, diffusion) : nothing
    Φ = vars.dynamics.geopotential      # includes the orography
    pₛ = vars.parameterizations.surface_pressure[ij]
    Φₛ = g * model.orography.orography[ij]     # surface geopotential

    K_momentum = vars.parameterizations.vertical_diffusivity_momentum
    K_heat = vars.parameterizations.vertical_diffusivity_heat
    (; density_over_distance) = vars.scratch.grid
    NF = eltype(K_momentum)

    # SURFACE LAYER: friction velocity from the drag of the boundary layer, surface Richardson number
    # from the surface temperatures (land and ocean, weighted by the land fraction)
    Φ_N = Φ[ij, nlayers]
    Ri_land, Ri_ocean = bulk_richardson_surface(ij, Φ_N - Φₛ, vars, atmosphere, diffusion, time_stepping)
    land_fraction = model.land_sea_mask.land_fraction[ij]
    Ri_surface = land_fraction * Ri_land + (1 - land_fraction) * Ri_ocean
    drag = vars.parameterizations.boundary_layer_drag_momentum[ij]
    friction_velocity = sqrt(drag) * vars.parameterizations.surface_wind_speed[ij]

    # roughness length for the stability correction, from the surface roughness scheme if present
    z₀ = haskey(vars.parameterizations, :momentum_roughness) ?
        vars.parameterizations.momentum_roughness[ij] : diffusion.roughness_length
    z₀ = ifelse(z₀ > 0, z₀, diffusion.roughness_length)

    # BOUNDARY LAYER HEIGHT h [m] above the surface
    Tᵥ_N = virtual_temperature(T[ij, nlayers], column_humidity(q, ij, nlayers), atmosphere)
    h = boundary_layer_height(ij, u, v, T, q, Φ, Φₛ, Tᵥ_N, g, cₚ, Ri_c, atmosphere)
    vars.parameterizations.boundary_layer_height[ij] = h

    # Stable correction of Frierson 2006, eq. 20, [1 + Ri/Ri_c log(z/z₀)/(1 - Ri/Ri_c)]⁻¹ with the surface
    # Richardson number, here written as (1 - r)/(1 - r + r log(z/z₀)), r = Ri/Ri_c ∈ [0, 1], to avoid a
    # division by zero at Ri = Ri_c. It is 1 for unstable conditions (Ri ≤ 0) and 0 for Ri ≥ Ri_c
    r = clamp(Ri_surface / Ri_c, 0, 1)
    surface_layer_height = fb * h
    stability_surface_layer = stable_correction(r, surface_layer_height, z₀)

    Tᵥ_above = virtual_temperature(T[ij, 1], column_humidity(q, ij, 1), atmosphere)
    for k in 1:(nlayers - 1)
        # interface k+½ between layer k (above) and k+1 (below)
        Tᵥ_below = virtual_temperature(T[ij, k + 1], column_humidity(q, ij, k + 1), atmosphere)
        ΔΦ = Φ[ij, k] - Φ[ij, k + 1]                    # geopotential difference between layers
        Δz = ΔΦ / g                                     # distance between layers [m]
        z = ((Φ[ij, k] + Φ[ij, k + 1]) / 2 - Φₛ) / g    # height of the interface above the surface [m]
        Tᵥ = (Tᵥ_above + Tᵥ_below) / 2                  # virtual temperature at the interface

        # air density at the interface over the distance between layers [kg/m⁴]
        p = pressure_half(k + 1, pₛ, vertical_coordinates)
        density_over_distance[ij, k] = p / (R_dry * Tᵥ * Δz)

        # 1. BOUNDARY LAYER PROFILE, Frierson 2006, eqs. 16-20, Kb(z) = κ u⋆ z in the surface layer
        # z < fb*h, cubic decay to zero at z = h above, zero above the boundary layer
        in_surface_layer = z < surface_layer_height
        stability = ifelse(in_surface_layer, stable_correction(r, z, z₀), stability_surface_layer)
        above_surface_layer = (z > surface_layer_height) & (z < h)
        ξ = ifelse(above_surface_layer, z - surface_layer_height, zero(z)) /
            ifelse(above_surface_layer, (1 - fb) * h, one(z))
        K_profile = κ * friction_velocity * z * stability * (1 - ξ)^2
        K_profile = ifelse(z < h, K_profile, zero(K_profile))

        # 2. LOCAL DIFFUSIVITY, l² S f(Ri) with S = |∂V/∂z|, squared buoyancy frequency
        # N² = g/(cₚTᵥ) ∂sᵥ/∂z from the virtual dry static energy sᵥ = cₚTᵥ + Φ
        ΔV² = (u[ij, k] - u[ij, k + 1])^2 + (v[ij, k] - v[ij, k + 1])^2
        S² = max(ΔV², NF(MINIMUM_WIND_SHEAR²)) / Δz^2
        S = sqrt(S²)
        N² = g * (cₚ * (Tᵥ_above - Tᵥ_below) + ΔΦ) / (cₚ * Tᵥ * Δz)
        Ri = max(N² / S², 0)      # gradient Richardson number, stable part

        # asymptotic mixing length λ, constant below mixing_length_height, then decaying towards
        # free_mixing_length (Holtslag and Boville, 1993, eq. 3.6), mixing length 1/l = 1/(κz) + 1/λ
        λ = free_mixing_length + (mixing_length - free_mixing_length) * exp(min(1 - z / mixing_length_height, 0))
        l = κ * z * λ / (κ * z + λ)

        # unstable: l²S√(1 - 18Ri) = l²√(S² - 18N²) (Holtslag and Boville, 1993, eq. 3.7), split into a shear
        # part l²S and a buoyancy part l_c²(√(S² - 18N²) - S) whose mixing length l_c uses at least the
        # convective_mixing_length, so that superadiabatic layers are mixed within a few time steps also in
        # the free atmosphere where λ is small. Continuous at N² = 0, l_c = l reproduces eq. 3.7
        λ_convective = max(λ, convective_mixing_length)
        l_convective = κ * z * λ_convective / (κ * z + λ_convective)
        unstable = N² < 0
        K_unstable = l^2 * S + l_convective^2 * (sqrt(S² - 18 * min(N², 0)) - S)
        K_local_momentum = ifelse(unstable, K_unstable, l^2 * S * stability_function_momentum(Ri))
        K_local_heat = ifelse(unstable, K_unstable, l^2 * S * stability_function_heat(Ri))

        # 3. BACKGROUND DIFFUSIVITY
        K_background = background_diffusivity * exp(-z / background_height)

        K_momentum[ij, k] = min(max(K_profile, K_local_momentum, K_background), max_diffusivity)
        K_heat[ij, k] = min(max(K_profile, K_local_heat, K_background), max_diffusivity)

        Tᵥ_above = Tᵥ_below
    end

    # no diffusive flux through the surface, which is handled by the surface fluxes
    K_momentum[ij, nlayers] = 0
    K_heat[ij, nlayers] = 0
    density_over_distance[ij, nlayers] = 0
    return nothing
end

# minimum wind speed [m/s] squared for the bulk Richardson number (Holtslag and Boville, 1993)
const MINIMUM_WIND_SPEED² = 1

# minimum wind speed difference between layers [m/s] squared for the shear in the local diffusivity
const MINIMUM_WIND_SHEAR² = 0.01

"""$(TYPEDSIGNATURES)
Stability correction [1 + Ri/Ri_c log(z/z₀)/(1 - Ri/Ri_c)]⁻¹ of Frierson 2006, eq. 20, written with
r = Ri/Ri_c ∈ [0, 1] as (1 - r)/(1 - r + r log(z/z₀)) to be finite at r = 1. z is capped one order of
magnitude above the roughness length `z₀` as in `neutral_transfer_coefficient`."""
@inline stable_correction(r, z, z₀) = (1 - r) / (1 - r + r * log(max(z, 10z₀) / z₀))

"""$(TYPEDSIGNATURES)
Long-tail stability function for momentum of Louis et al. 1982, 1/(1 + 2bRi/√(1 + dRi)) with b = d = 5,
for the gradient Richardson number `Ri` ≥ 0."""
@inline stability_function_momentum(Ri) = 1 / (1 + 10Ri / sqrt(1 + 5Ri))

"""$(TYPEDSIGNATURES)
Long-tail stability function for heat of Louis et al. 1982, 1/(1 + 3bRi√(1 + dRi)) with b = d = 5,
for the gradient Richardson number `Ri` ≥ 0."""
@inline stability_function_heat(Ri) = 1 / (1 + 15Ri * sqrt(1 + 5Ri))

"""$(TYPEDSIGNATURES)
Boundary layer height h [m] above the surface of column `ij`: the height where the bulk Richardson
number relative to the lowermost layer N, Ri(z) = gz (sᵥ(z) - sᵥ_N) / (cₚTᵥ_N |V(z)|²) with the virtual
dry static energy sᵥ = cₚTᵥ + Φ (Frierson et al. 2006), first reaches `Ri_c` from below, interpolated
linearly between the layers that bracket Ri_c with Ri = 0 at the lowermost layer. The height of the
uppermost layer if Ri < Ri_c everywhere. Fixed trip count, no early exit."""
@propagate_inbounds function boundary_layer_height(ij, u, v, T, q, Φ, Φₛ, Tᵥ_N, g, cₚ, Ri_c, atmosphere)
    nlayers = size(T, 2)
    sᵥ_N = cₚ * Tᵥ_N + Φ[ij, nlayers]
    h = (Φ[ij, 1] - Φₛ) / g                 # top layer height if no crossing of Ri_c
    z_below = (Φ[ij, nlayers] - Φₛ) / g     # height of the layer below k, starting with N
    Ri_below = zero(h)                      # Ri relative to the lowermost layer is zero there
    found = false

    for k in (nlayers - 1):-1:1
        z = (Φ[ij, k] - Φₛ) / g
        sᵥ = cₚ * virtual_temperature(T[ij, k], column_humidity(q, ij, k), atmosphere) + Φ[ij, k]
        V² = max(u[ij, k]^2 + v[ij, k]^2, MINIMUM_WIND_SPEED²)
        Ri = g * z * (sᵥ - sᵥ_N) / (cₚ * Tᵥ_N * V²)

        # first crossing of Ri_c from below, Ri_below < Ri_c ≤ Ri, guard division for other layers
        crossing = !found & (Ri >= Ri_c)
        weight = ifelse(crossing, Ri_c - Ri_below, zero(Ri)) / ifelse(crossing, Ri - Ri_below, one(Ri))
        h = ifelse(crossing, z_below + weight * (z - z_below), h)
        found |= crossing

        z_below = z
        Ri_below = Ri
    end
    return h
end

# value of the variable in layer k that is diffused, with or without a geopotential term
@inline diffused_value(variable, ::Nothing, ij, k) = variable[ij, k]
@inline diffused_value(variable, (Φ, scale)::Tuple, ij, k) = variable[ij, k] + scale * Φ[ij, k]

"""$(TYPEDSIGNATURES)
Implicit (backward Euler) vertical diffusion of `variable` (plus `scale * Φ` with `geopotential_term =
(Φ, scale)`, or `nothing`) in column `ij` with diffusivity `K` [m²/s] on the interfaces below every layer.
With the diffusion operator L as a flux divergence on the pressure thickness Δp of the layers,

    (L x)ₖ = g/Δpₖ [aₖ₊½ (xₖ₊₁ - xₖ) - aₖ₋½ (xₖ - xₖ₋₁)],     aₖ₊½ = ρₖ₊½ Kₖ₊½ / Δzₖ₊½,

and no flux through the top and the surface, the tendency y = (x⁺ - x⁻)/Δt is the solution of the
tridiagonal system (1 - ΔtL) y = L x⁻ + S, solved with the Thomas algorithm. With `couple = true` the
source S is the `tendency` accumulated so far (e.g. the surface fluxes into the lowermost layer), which
is then spread over the column within the time step and `tendency` is replaced by y. Otherwise S = 0 and
y is added to `tendency`. Conserves the mass-weighted column integral Σ xₖ Δpₖ exactly."""
@propagate_inbounds function implicit_vertical_diffusion!(
        ij,
        tendency,                   # tendency of variable, accumulate into or replace
        variable,                   # variable to diffuse
        geopotential_term,          # (Φ, scale) to diffuse variable + scale*Φ, or nothing
        K,                          # diffusivity [m²/s] at the interface below every layer
        density_over_distance,      # ρ/Δz [kg/m⁴] at the interface below every layer
        tridiagonal_upper,          # scratch, upper diagonal after forward elimination
        tridiagonal_rhs,            # scratch, right-hand side after forward elimination
        pₛ,                         # surface pressure [Pa]
        Δt,                         # time step [s]
        g,                          # gravity [m/s²]
        vertical_coordinates,
        couple::Bool,               # use tendency as source in the implicit solve?
    )
    nlayers = size(variable, 2)
    NF = eltype(tendency)

    # FORWARD ELIMINATION, from top to surface
    a_above = zero(NF)              # conductance [kg/m²/s] at the interface above layer k, no flux through top
    x_above = diffused_value(variable, geopotential_term, ij, 1)
    upper_above = zero(NF)          # upper diagonal and rhs of the eliminated layer above
    rhs_above = zero(NF)

    for k in 1:nlayers
        x = diffused_value(variable, geopotential_term, ij, k)
        x_below = diffused_value(variable, geopotential_term, ij, min(k + 1, nlayers))
        a_below = ifelse(k < nlayers, density_over_distance[ij, k] * K[ij, k], zero(NF))

        g_Δp = g / pressure_thickness(k, pₛ, vertical_coordinates)
        explicit_tendency = g_Δp * (a_below * (x_below - x) - a_above * (x - x_above))
        source = ifelse(couple, tendency[ij, k], zero(NF))

        # row k of (1 - ΔtL): lower, diagonal and upper entry, diagonally dominant
        α = Δt * g_Δp
        lower = -α * a_above
        upper = -α * a_below
        diagonal = 1 + α * (a_above + a_below)

        pivot = diagonal - lower * upper_above          # ≥ 1
        upper_above = upper / pivot
        rhs_above = (explicit_tendency + source - lower * rhs_above) / pivot
        tridiagonal_upper[ij, k] = upper_above
        tridiagonal_rhs[ij, k] = rhs_above

        a_above = a_below
        x_above = x
    end

    # BACK SUBSTITUTION, from surface to top
    y_below = zero(NF)
    for k in nlayers:-1:1
        y = tridiagonal_rhs[ij, k] - tridiagonal_upper[ij, k] * y_below
        tendency[ij, k] = ifelse(couple, y, tendency[ij, k] + y)
        y_below = y
    end
    return nothing
end
