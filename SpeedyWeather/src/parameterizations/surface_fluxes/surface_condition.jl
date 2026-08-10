@parameterized @kwdef struct ConstantWindSlowdown{NF} <: AbstractParameterization
    "[OPTION] Ratio of near-surface wind to lowest-level wind [1]"
    @param wind_slowdown::NF = 0.95 (bounds = 0 .. 1,)
end

Adapt.@adapt_structure ConstantWindSlowdown
ConstantWindSlowdown(SG::SpectralGrid; kwargs...) = ConstantWindSlowdown{SG.NF}(; kwargs...)
initialize!(::ConstantWindSlowdown, ::PrimitiveEquation) = nothing

@propagate_inbounds function get_slowdown!(ij, vars, ws::ConstantWindSlowdown, model)
    return ws.wind_slowdown
end

"""Calculates the near-surface wind speed from the lowest model
level wind speed using surface roughness. Fields are $(TYPEDFIELDS)"""
@kwdef struct WindSlowdown{NF} <: AbstractParameterization

end

Adapt.@adapt_structure WindSlowdown
WindSlowdown(SG::SpectralGrid; kwargs...) = WindSlowdown{SG.NF}(; kwargs...)
initialize!(::WindSlowdown, ::PrimitiveEquation) = nothing

@propagate_inbounds function get_slowdown!(ij, vars, ws::WindSlowdown, model)
    g = model.planet.gravity
    (; nlayers) = model.geometry

    Φ₁ = vars.dynamics.geopotential[ij, nlayers]
    z_surf = model.orography.orography[ij]

    z₁ = (Φ₁ / g) - z_surf
    z₀ = vars.parameterizations.momentum_roughness[ij]

    slowdown = log(10 / z₀) / log(z₁ / z₀)

    return slowdown
end

export SurfaceCondition

"""Surface condition parameterization that calculates near-surface atmospheric
variables needed for surface flux calculations. Computes surface wind speed
including sub-grid scale gusts, surface air density, and surface air temperature
by extrapolating from the lowest model level to the surface using standard
atmospheric relationships. Fields are $(TYPEDFIELDS)"""
@kwdef struct SurfaceCondition{NF, WS} <: AbstractParameterization
    "[OPTION] Ratio of near-surface wind to lowest-level wind [1]"
    wind_slowdown::WS

    "[OPTION] Wind speed of sub-grid scale gusts [m/s]"
    gust_speed::NF = 1

    "[OPTION] Beljaars (1994) free convection parameter β [1]"
    β::NF = 1.2

    "[OPTION] Boundary layer depth for free convection scaling [m]"
    z_i::NF = 1000.0

    "[OPTION] Typical heat transfer coefficient for scaling [1]"
    C_H::NF = 1.0e-3

    "[OPTION] Minimum background gustiness [m/s]"
    min_gust::NF = 0.5
end

Adapt.@adapt_structure SurfaceCondition

function SurfaceCondition(
        SG::SpectralGrid;
        wind_slowdown = WindSlowdown(SG),
        kwargs...
    )
    NF = SG.NF
    WS = typeof(wind_slowdown)
    return SurfaceCondition{NF, WS}(; wind_slowdown = wind_slowdown, kwargs...)
end

initialize!(::SurfaceCondition, ::PrimitiveEquation) = nothing

# function barrier
@propagate_inbounds function parameterization!(ij, vars, sc::SurfaceCondition, model)
    return surface_condition!(ij, vars, sc, model)
end

# @propagate_inbounds function surface_condition!(ij, vars, surface_condition::SurfaceCondition, model)

#     (; wind_slowdown, gust_speed) = surface_condition
#     (; nlayers) = model.geometry
#     coord = model.geometry.vertical_coordinates
#     (; atmosphere) = model

#     # Fortran SPEEDY documentation eq. 49 but use previous time step for numerical stability
#     u_grid = get_prognostic_step(vars.grid.u, model.time_stepping, surface_condition)
#     v_grid = get_prognostic_step(vars.grid.v, model.time_stepping, surface_condition)

#     wind_slowdown_ratio = get_slowdown!(ij, vars, wind_slowdown, model)

#     uₛ = wind_slowdown_ratio * u_grid[ij, nlayers]
#     vₛ = wind_slowdown_ratio * v_grid[ij, nlayers]

#     # Fortran SPEEDY documentation eq. 50, sqrt(u² + v² + gust_speed²)
#     surface_wind_speed = sqrt(muladd(uₛ, uₛ, muladd(vₛ, vₛ, gust_speed^2)))
#     vars.parameterizations.surface_wind_speed[ij] = surface_wind_speed

#     # Surface air density
#     (; surface_air_density) = vars.parameterizations
#     temperature = get_prognostic_step(vars.grid.temperature, model.time_stepping, surface_condition)
#     pₛ = vars.parameterizations.surface_pressure[ij] # surface pressure [Pa]
#     (; R_dry, κ) = model.atmosphere

#     σ = pressure(nlayers, pₛ, coord) / pₛ           # σ vertical coordinate at lowest model level
#     T = temperature[ij, nlayers]                    # virtual temperature at lowest model level [K]
#     q = haskey(vars.grid, :humidity) ?              # specific humidity at lowest model level [kg/kg]
#         get_prognostic_step(vars.grid.humidity, model.time_stepping, surface_condition)[ij, nlayers] : zero(T)
#     Tᵥ = virtual_temperature(T, q, atmosphere)      # virtual temperature at lowest model level [K]
#     σ⁻ᵏ = σ^(-κ)                                    # precalculate
#     Tᵥ *= σ⁻ᵏ                                       # lower to surface assuming dry adiabatic lapse rate
#     ρ = pₛ / (R_dry * Tᵥ)                           # surface air density [kg/m³] from ideal gas law
#     surface_air_density[ij] = ρ                     # store for surface temp/humidity fluxes

#     # Surface air temperature
#     (; surface_air_temperature) = vars.parameterizations
#     T *= σ⁻ᵏ                                        # lower to surface assuming dry adiabatic lapse rate
#     surface_air_temperature[ij] = T                 # store for surface temp/humidity fluxes
#     return nothing
# end

@propagate_inbounds function surface_condition!(ij, vars, surface_condition::SurfaceCondition, model)

    (; wind_slowdown, β, z_i, C_H, min_gust) = surface_condition
    (; nlayers) = model.geometry
    coord = model.geometry.vertical_coordinates
    (; atmosphere, planet, time_stepping, land_sea_mask) = model

    # 1. Lower air temperature to surface assuming dry adiabatic lapse rate
    temperature = get_prognostic_step(vars.grid.temperature, time_stepping, surface_condition)
    pₛ = vars.parameterizations.surface_pressure[ij]
    (; R_dry, κ) = atmosphere

    σ = pressure(nlayers, pₛ, coord) / pₛ
    T = temperature[ij, nlayers]
    q = haskey(vars.grid, :humidity) ?
        get_prognostic_step(vars.grid.humidity, time_stepping, surface_condition)[ij, nlayers] : zero(T)

    Tᵥ = virtual_temperature(T, q, atmosphere)
    σ⁻ᵏ = σ^(-κ)
    Tᵥ *= σ⁻ᵏ
    T *= σ⁻ᵏ

    vars.parameterizations.surface_air_temperature[ij] = T
    vars.parameterizations.surface_air_density[ij] = pₛ / (R_dry * Tᵥ)

    # 2. Get underlying surface skin temperature
    land_fraction = land_sea_mask.land_fraction[ij]
    SST = get_prognostic_step(vars.prognostic.ocean.sea_surface_temperature, time_stepping, surface_condition)[ij]
    T_land = vars.prognostic.land.soil_temperature[ij, 1]

    # Area-weighted skin temperature
    T_skin = SST * (1 - land_fraction) + T_land * land_fraction

    # 3. Compute skin virtual temperature
    # (Assuming saturated at surface to calculate maximum buoyancy potential)
    q_skin = saturation_humidity(T_skin, pₛ, atmosphere)
    Tᵥ_skin = virtual_temperature(T_skin, q_skin, atmosphere)

    # 4. Explicit Free Convection Gustiness (IFS Beljaars Bulk Approximation)
    ΔTᵥ = max(Tᵥ_skin - Tᵥ, 0) # Only unstable conditions generate convective gusts

    g = planet.gravity
    W_gust = sqrt((β^3) * C_H * z_i * (g / max(Tᵥ, 200)) * ΔTᵥ)
    W_gust = max(W_gust, min_gust) # Maintain a small background turbulence

    # 5. Final Effective Wind Speed Calculation
    u_grid = get_prognostic_step(vars.grid.u, time_stepping, surface_condition)
    v_grid = get_prognostic_step(vars.grid.v, time_stepping, surface_condition)

    wind_slowdown_ratio = get_slowdown!(ij, vars, wind_slowdown, model)

    uₛ = wind_slowdown_ratio * u_grid[ij, nlayers]
    vₛ = wind_slowdown_ratio * v_grid[ij, nlayers]

    # U_eff = sqrt(u² + v² + W_gust²)
    vars.parameterizations.surface_wind_speed[ij] = sqrt(muladd(uₛ, uₛ, muladd(vₛ, vₛ, W_gust^2)))

    return nothing
end
