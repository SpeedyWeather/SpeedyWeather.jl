export ConstantWindSlowdown, LearnedWindSlowdown, ConstantWindGust, BeljaarsWindGust, SurfaceCondition

@parameterized @kwdef struct ConstantWindSlowdown{NF} <: AbstractParameterization
    "[OPTION] Ratio of near-surface wind to lowest-level wind [1]"
    @param wind_slowdown::NF = 0.95 (bounds = 0 .. 1,)
end

Adapt.@adapt_structure ConstantWindSlowdown
ConstantWindSlowdown(SG::SpectralGrid; kwargs...) = ConstantWindSlowdown{SG.NF}(; kwargs...)
initialize!(::ConstantWindSlowdown, ::PrimitiveEquation) = nothing

@propagate_inbounds function calc_wind_slowdown(ij, vars, ws::ConstantWindSlowdown, model, surface_condition, U₁)
    return ws.wind_slowdown
end

"""Calculates the near-surface wind speed from the lowest model
level wind speed using surface roughness. Fields are $(TYPEDFIELDS)"""
@kwdef struct LearnedWindSlowdown{NF} <: AbstractParameterization
    # Ocean slowdown constants
    @param c1::NF = 0.14511
    @param c2::NF = 1.9883
    @param c3::NF = -10.904
    @param c4::NF = 0.34863

    # Land slowdown constants
    @param c5::NF = 4.2599
    @param c6::NF = 1.9678
    @param c7::NF = -11.839
    @param c8::NF = 9.6687
    @param c9::NF = -0.29968
    @param c10::NF = 2.3145
    @param c11::NF = 10.215
end

function Adapt.adapt_structure(to, obj::LearnedWindSlowdown{NF}) where {NF}
    return LearnedWindSlowdown{NF}(;
        (
            f => Adapt.adapt(to, getfield(obj, f))
                for f in fieldnames(typeof(obj))
        )...
    )
end

LearnedWindSlowdown(SG::SpectralGrid; kwargs...) = LearnedWindSlowdown{SG.NF}(; kwargs...)
initialize!(::LearnedWindSlowdown, ::PrimitiveEquation) = nothing

@propagate_inbounds function calc_wind_slowdown(ij, vars, ws::LearnedWindSlowdown, model, surface_condition, U₁)
    (; gravity) = model.planet
    (; nlayers) = model.geometry
    (; land_sea_mask, atmosphere, time_stepping) = model
    cₚ = atmosphere.heat_capacity
    Φ₁ = vars.dynamics.geopotential[ij, nlayers]
    z_surf = model.orography.orography[ij]
    land_fraction = land_sea_mask.land_fraction[ij]

    z₁ = (Φ₁ / gravity) - z_surf
    z₀_raw_land = vars.parameterizations.land.momentum_roughness[ij]
    z₀_raw_ocean = vars.parameterizations.ocean.momentum_roughness[ij]

    # Handle initialisation independently for land and ocean masks
    z₀_land = z₀_raw_land > 0 ? z₀_raw_land : convert(typeof(z₀_raw_land), 0.1)
    z₀_ocean = z₀_raw_ocean > 0 ? z₀_raw_ocean : convert(typeof(z₀_raw_ocean), 0.0001)

    T₁ = get_prognostic_step(vars.grid.temperature, time_stepping, ws)[ij, nlayers]
    q₁ = hasproperty(vars.grid, :humidity) ?
        get_prognostic_step(vars.grid.humidity, time_stepping, ws)[ij, nlayers] : zero(T₁)

    Tᵥ₁ = virtual_temperature(T₁, q₁, atmosphere)

    T_skin_ocean = get_prognostic_step(vars.prognostic.ocean.sea_surface_temperature, time_stepping, ws)[ij]
    T_skin_land = vars.prognostic.land.soil_temperature[ij, 1]

    pₛ = vars.parameterizations.surface_pressure[ij]
    q_skin_land = saturation_humidity(T_skin_land, pₛ, atmosphere)
    q_skin_ocean = saturation_humidity(T_skin_ocean, pₛ, atmosphere)

    Tᵥ_skin_land = virtual_temperature(T_skin_land, q_skin_land, atmosphere)
    Tᵥ_skin_ocean = virtual_temperature(T_skin_ocean, q_skin_ocean, atmosphere)
    Δθᵥ_land = (Tᵥ₁ + (gravity / cₚ) * z₁) - Tᵥ_skin_land
    Δθᵥ_ocean = (Tᵥ₁ + (gravity / cₚ) * z₁) - Tᵥ_skin_ocean

    log_z₀_ocean = log(z₀_ocean)
    log_z₀_land = log(z₀_land)

    T_grad_ocean = vars.parameterizations.surface_air_temperature[ij] - T_skin_ocean

    # Unpack constants for learned slowdown parameterization
    (; c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11) = ws

    signed_sqrt(x) = x >= 0 ? sqrt(x) : -sqrt(-x)
    ocean_slowdown = ((atan((log_z₀_ocean + ((log_z₀_ocean + log_z₀_ocean) + log_z₀_ocean)) + ((z₁ - exp(T_grad_ocean)) * c1)) - c2) * (c3 - log_z₀_ocean)) - signed_sqrt(c4 * (log_z₀_ocean - atan((Δθᵥ_ocean + c3) + (T_grad_ocean + U₁))))
    land_slowdown = c5 + ((c6 * atan(exp((c7 - signed_sqrt(U₁)) - (signed_sqrt(Δθᵥ_land + (c8 + log_z₀_land)) - (Δθᵥ_land * ((Δθᵥ_land + (z₁ + Δθᵥ_land))^c9)))) - (log_z₀_land * c10)))) * (c11 + log_z₀_land)

    slowdown = ocean_slowdown * (1 - land_fraction) + land_slowdown * land_fraction

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
    @param min_gust::NF = 0.0
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

    U₁ = sqrt(muladd(u_grid[ij, nlayers], u_grid[ij, nlayers], v_grid[ij, nlayers]^2))

    wind_slowdown = calc_wind_slowdown(ij, vars, wind_slowdown, model, surface_condition, U₁)

    uₛ = wind_slowdown * u_grid[ij, nlayers]
    vₛ = wind_slowdown * v_grid[ij, nlayers]

    gust = calc_wind_gust(ij, vars, wind_gust, model, surface_condition, pₛ, Tᵥ)

    vars.parameterizations.surface_wind_speed[ij] = sqrt(muladd(uₛ, uₛ, muladd(vₛ, vₛ, gust^2)))

    return nothing
end
