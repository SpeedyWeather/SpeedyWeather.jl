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

"""
Learned near-surface wind slowdown parameterization using a symbolic
regression expression.

The slowdown factor is calculated as

    f = clamp(1 - s, 0, 1)

where

    s = c1*U₁ - c2*U₁/z₁ + c3*Δθᵥ +
        c4*log(z₀) + c5*log(z₁) + c6

Here U₁ is the lowest-level wind speed, z₁ is its height above
the surface, z₀ is the momentum roughness length, and Δθᵥ is
the virtual potential temperature difference between the lowest
model level and the surface.

Fields are $(TYPEDFIELDS)
"""
@parameterized @kwdef struct LearnedWindSlowdown{NF} <: AbstractParameterization
    @param c1::NF = 0.0059787630088485833
    @param c2::NF = 0.23547083241598257
    @param c3::NF = 0.0099666424032772266
    @param c4::NF = 0.03703253707251629
    @param c5::NF = 0.02252190141813104
    @param c6::NF = 0.3397981681230896
end

Adapt.@adapt_structure LearnedWindSlowdown

LearnedWindSlowdown(SG::SpectralGrid; kwargs...) = LearnedWindSlowdown{SG.NF}(; kwargs...)
initialize!(::LearnedWindSlowdown, ::PrimitiveEquation) = nothing

@inline function wind_slowdown(ws::LearnedWindSlowdown, U₁, z₁, log_z₀, Δθᵥ)
    (; c1, c2, c3, c4, c5, c6) = ws

    # Prevent division by zero and invalid logarithms
    z = max(z₁, one(z₁))

    # Learned slowdown expression using fused multiply-add
    s = muladd(c1, U₁, c6)
    s = muladd(-c2, U₁ / z, s)
    s = muladd(c3, Δθᵥ, s)
    s = muladd(c4, log_z₀, s)
    s = muladd(c5, log(z), s)

    return clamp(one(s) - s, zero(s), one(s))
end

@propagate_inbounds function calc_wind_slowdown(
        ij, vars, ws::LearnedWindSlowdown,
        model, surface_condition, U₁
    )
    (; gravity) = model.planet
    (; nlayers) = model.geometry
    (; land_sea_mask, atmosphere, time_stepping) = model

    land_fraction = land_sea_mask.land_fraction[ij]
    cₚ = atmosphere.heat_capacity

    # Reference height
    Φ₁ = vars.dynamics.geopotential[ij, nlayers]
    z_surf = model.orography.orography[ij]
    z₁ = Φ₁ / gravity - z_surf

    # Combined momentum roughness
    z₀_raw = vars.parameterizations.momentum_roughness[ij]

    z₀_default = land_fraction * oftype(z₀_raw, 0.1) + (1 - land_fraction) * oftype(z₀_raw, 0.0001)
    z₀ = z₀_raw > 0 ? z₀_raw : z₀_default

    # Lowest model-level temperature and humidity
    T₁ = get_prognostic_step(
        vars.grid.temperature, time_stepping, ws
    )[ij, nlayers]

    q₁ = haskey(vars.grid, :humidity) ?
        get_prognostic_step(
            vars.grid.humidity, time_stepping, ws
        )[ij, nlayers] : zero(T₁)

    Tᵥ₁ = virtual_temperature(T₁, q₁, atmosphere)

    T_skin = zero(T₁)
    if land_fraction > 0
        land_temp = get_prognostic_step(
            vars.prognostic.land.soil_temperature,
            time_stepper(time_stepping, :land), ws,
        )
        T_skin += land_fraction * land_temp[ij, 1]
    end
    if land_fraction < 1
        ocean_temp = get_prognostic_step(
            vars.prognostic.ocean.sea_surface_temperature,
            time_stepper(time_stepping, :ocean), ws,
        )
        T_skin += (one(land_fraction) - land_fraction) * ocean_temp[ij]
    end

    # Surface virtual temperature
    pₛ = vars.parameterizations.surface_pressure[ij]
    q_skin = saturation_humidity(T_skin, pₛ, atmosphere)
    Tᵥ_skin = virtual_temperature(T_skin, q_skin, atmosphere)

    # Virtual potential temperature difference
    Δθᵥ = Tᵥ₁ + (gravity / cₚ) * z₁ - Tᵥ_skin

    return wind_slowdown(ws, U₁, z₁, log(z₀), Δθᵥ)

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
    @param C_H::NF = 1.0e-3 #TODO: use real heat transfer coefficient calculated in boundary_layer

    "[OPTION] Minimum background gustiness [m/s]"
    @param min_gust::NF = 0.1
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
@parameterized @kwdef struct SurfaceCondition{WS, WG} <: AbstractParameterization
    "[OPTION] Ratio of near-surface wind to lowest-level wind [1]"
    @component wind_slowdown::WS

    "[OPTION] Wind speed of sub-grid scale gusts [m/s]"
    @component wind_gust::WG
end

Adapt.@adapt_structure SurfaceCondition

function SurfaceCondition(
        SG::SpectralGrid;
        wind_slowdown = LearnedWindSlowdown(SG),
        wind_gust = BeljaarsWindGust(SG)
    )
    return SurfaceCondition(wind_slowdown, wind_gust)
end

initialize!(::SurfaceCondition, ::PrimitiveEquation) = nothing

# function barrier
@propagate_inbounds function parameterization!(ij, vars, sc::SurfaceCondition, model)
    return surface_condition!(ij, vars, sc, model)
end

variables(::SurfaceCondition) = (
    ParameterizationVariable(:surface_wind_u, GridXY(), desc = "Mean zonal surface wind", units = "m/s"),
    ParameterizationVariable(:surface_wind_v, GridXY(), desc = "Mean meridional surface wind", units = "m/s"),
)

@propagate_inbounds function surface_condition!(ij, vars, surface_condition::SurfaceCondition, model)

    (; wind_slowdown, wind_gust) = surface_condition
    (; nlayers) = model.geometry
    coord = model.geometry.vertical_coordinates
    (; atmosphere, time_stepping) = model

    # 1. Lower air temperature to surface assuming dry adiabatic lapse rate
    temperature = get_prognostic_step(vars.grid.temperature, time_stepping, surface_condition)
    pₛ = vars.parameterizations.surface_pressure[ij]
    (; R_dry, κ) = atmosphere

    σ = pressure(nlayers, pₛ, coord) / pₛ           # σ vertical coordinate at lowest model level
    T = temperature[ij, nlayers]                    # virtual temperature at lowest model level [K]
    q = haskey(vars.grid, :humidity) ?              # specific humidity at lowest model level [kg/kg]
        get_prognostic_step(vars.grid.humidity, model.time_stepping, surface_condition)[ij, nlayers] : zero(T)

    Tᵥ = virtual_temperature(T, q, atmosphere)      # virtual temperature at lowest model level [K]
    σ⁻ᵏ = σ^(-κ)                                    # precalculate lapse rate factor
    Tᵥ *= σ⁻ᵏ                                       # lower to surface assuming dry adiabatic lapse rate
    T *= σ⁻ᵏ
    vars.parameterizations.surface_air_temperature[ij] = T              # store for flux computation
    vars.parameterizations.surface_air_density[ij] = pₛ / (R_dry * Tᵥ)  # surface air density [kg/m³] from ideal gas law

    # 2. Surface wind
    u_grid = get_prognostic_step(vars.grid.u, time_stepping, surface_condition)
    v_grid = get_prognostic_step(vars.grid.v, time_stepping, surface_condition)

    # lowermost-layer wind speed sqrt(u^2 + v^2) but written with fused muladd
    Uₙ = sqrt(muladd(u_grid[ij, nlayers], u_grid[ij, nlayers], v_grid[ij, nlayers]^2))

    wind_slowdown = calc_wind_slowdown(ij, vars, wind_slowdown, model, surface_condition, Uₙ)

    uₛ = wind_slowdown * u_grid[ij, nlayers]        # surface wind components
    vₛ = wind_slowdown * v_grid[ij, nlayers]
    vars.parameterizations.surface_wind_u[ij] = uₛ  # store for fluxes
    vars.parameterizations.surface_wind_v[ij] = vₛ

    # gusts are added for surface wind speed but not on individual components
    gust = calc_wind_gust(ij, vars, wind_gust, model, surface_condition, pₛ, Tᵥ)

    vars.parameterizations.surface_wind_speed[ij] = sqrt(muladd(uₛ, uₛ, muladd(vₛ, vₛ, gust^2)))

    return nothing
end
