@parameterized @kwdef struct ConstantWindSlowdown{NF} <: AbstractParameterization
    "[OPTION] Ratio of near-surface wind to lowest-level wind [1]"
    @param wind_slowdown::NF = 0.95 (bounds = 0 .. 1,)
end

Adapt.@adapt_structure ConstantWindSlowdown
ConstantWindSlowdown(SG::SpectralGrid; kwargs...) = ConstantWindSlowdown{SG.NF}(; kwargs...)
initialize!(::ConstantWindSlowdown, ::PrimitiveEquation) = nothing

@propagate_inbounds function calc_wind_slowdown(ij, vars, ws::ConstantWindSlowdown, model)
    return ws.wind_slowdown
end

"""Calculates the near-surface wind speed from the lowest model
level wind speed using surface roughness. Fields are $(TYPEDFIELDS)"""
@kwdef struct LearnedWindSlowdown{NF} <: AbstractParameterization

end

Adapt.@adapt_structure LearnedWindSlowdown
LearnedWindSlowdown(SG::SpectralGrid; kwargs...) = LearnedWindSlowdown{SG.NF}(; kwargs...)
initialize!(::LearnedWindSlowdown, ::PrimitiveEquation) = nothing

@propagate_inbounds function calc_wind_slowdown(ij, vars, ws::LearnedWindSlowdown, model)
    g = model.planet.gravity
    (; nlayers) = model.geometry

    Φ₁ = vars.dynamics.geopotential[ij, nlayers]
    z_surf = model.orography.orography[ij]

    z₁ = (Φ₁ / g) - z_surf
    z₀ = vars.parameterizations.momentum_roughness[ij]

    # TODO: add multiplier based on profile from lowest model level to surface
    slowdown = log(convert(typeof(z₀), 10) / z₀) / log(z₁ / z₀) # * φ

    return slowdown
end

@parameterized @kwdef struct ConstantWindGust{NF} <: AbstractParameterization
    "[OPTION] Wind speed of sub-grid scale gusts [m/s]"
    @param wind_gust::NF = 0.5 (bounds = 0 .. 10,)
end

Adapt.@adapt_structure ConstantWindGust
ConstantWindGust(SG::SpectralGrid; kwargs...) = ConstantWindGust{SG.NF}(; kwargs...)
initialize!(::ConstantWindGust, ::PrimitiveEquation) = nothing

@propagate_inbounds function calc_wind_gust(ij, vars, wg::ConstantWindGust, model, surface_condition, pₛ, Tᵥ)
    return wg.wind_gust
end

@parameterized @kwdef struct BeljaarsWindGust{NF} <: AbstractParameterization
    "[OPTION] Beljaars (1994) free convection parameter β [1]"
    @param β::NF = 1.2

    "[OPTION] Boundary layer depth for free convection scaling [m]"
    @param z_i::NF = 1000.0

    "[OPTION] Typical heat transfer coefficient for scaling [1]"
    @param C_H::NF = 1.0e-3

    "[OPTION] Minimum background gustiness [m/s]"
    @param min_gust::NF = 0.5
end

Adapt.@adapt_structure BeljaarsWindGust
BeljaarsWindGust(SG::SpectralGrid; kwargs...) = BeljaarsWindGust{SG.NF}(; kwargs...)
initialize!(::BeljaarsWindGust, ::PrimitiveEquation) = nothing

@propagate_inbounds function calc_wind_gust(ij, vars, wg::BeljaarsWindGust, model, surface_condition, pₛ, Tᵥ)
    (; β, z_i, C_H, min_gust) = wg
    (; gravity) = model.planet
    (; land_sea_mask, time_stepping, atmosphere) = model

    # Calculate skin temperature #TODO implement actual skin temperature parameterization
    land_fraction = land_sea_mask.land_fraction[ij]
    SST = get_prognostic_step(vars.prognostic.ocean.sea_surface_temperature, time_stepping, surface_condition)[ij]
    T_land = vars.prognostic.land.soil_temperature[ij, 1]
    T_skin = SST * (1 - land_fraction) + T_land * land_fraction

    # Calculate skin virtual temperature
    q_skin = saturation_humidity(T_skin, pₛ, atmosphere)
    Tᵥ_skin = virtual_temperature(T_skin, q_skin, atmosphere)
    # Virtual temperature gradient
    ΔTᵥ = max(Tᵥ_skin - Tᵥ, zero(Tᵥ)) # Only unstable conditions generate convective gusts

    # Wind gust approximation based on Beljaars (1994) free convection scaling
    gust = sqrt((β^3) * C_H * z_i * (gravity / max(Tᵥ, convert(typeof(Tᵥ), 200))) * ΔTᵥ)
    gust = max(gust, min_gust) # Maintain a small background turbulence

    return gust
end

export SurfaceCondition

"""Surface condition parameterization that calculates near-surface atmospheric
variables needed for surface flux calculations. Computes surface wind speed
including sub-grid scale gusts, surface air density, and surface air temperature
by extrapolating from the lowest model level to the surface using standard
atmospheric relationships. Fields are $(TYPEDFIELDS)"""
@kwdef struct SurfaceCondition{NF, WS, WG} <: AbstractParameterization
    "[OPTION] Ratio of near-surface wind to lowest-level wind [1]"
    wind_slowdown::WS

    "[OPTION] Wind speed of sub-grid scale gusts [m/s]"
    wind_gust::WG
end

Adapt.@adapt_structure SurfaceCondition

function SurfaceCondition(
        SG::SpectralGrid;
        wind_slowdown = LearnedWindSlowdown(SG),
        wind_gust = BeljaarsWindGust(SG),
        kwargs...
    )
    NF = SG.NF
    WS = typeof(wind_slowdown)
    WG = typeof(wind_gust)
    return SurfaceCondition{NF, WS, WG}(
        ;
        wind_slowdown = wind_slowdown,
        wind_gust = wind_gust,
        kwargs...
    )
end

initialize!(::SurfaceCondition, ::PrimitiveEquation) = nothing

# function barrier
@propagate_inbounds function parameterization!(ij, vars, sc::SurfaceCondition, model)
    return surface_condition!(ij, vars, sc, model)
end

@propagate_inbounds function surface_condition!(ij, vars, surface_condition::SurfaceCondition, model)

    (; wind_slowdown, wind_gust) = surface_condition
    (; nlayers) = model.geometry
    coord = model.geometry.vertical_coordinates
    (; atmosphere, time_stepping) = model

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

    u_grid = get_prognostic_step(vars.grid.u, time_stepping, surface_condition)
    v_grid = get_prognostic_step(vars.grid.v, time_stepping, surface_condition)

    wind_slowdown_ratio = calc_wind_slowdown(ij, vars, wind_slowdown, model)

    uₛ = wind_slowdown_ratio * u_grid[ij, nlayers]
    vₛ = wind_slowdown_ratio * v_grid[ij, nlayers]

    gust = calc_wind_gust(ij, vars, wind_gust, model, surface_condition, pₛ, Tᵥ)

    vars.parameterizations.surface_wind_speed[ij] = sqrt(muladd(uₛ, uₛ, muladd(vₛ, vₛ, gust^2)))

    return nothing
end
