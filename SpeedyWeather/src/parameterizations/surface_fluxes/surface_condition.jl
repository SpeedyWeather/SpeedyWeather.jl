abstract type AbstractSurfaceCondition <: AbstractParameterization end

export SurfaceCondition
export OceanNeutralWindSpeed

@kwdef struct OceanNeutralWindSpeed{NF} <: AbstractParameterization end

Adapt.@adapt_structure OceanNeutralWindSpeed

OceanNeutralWindSpeed(SG::SpectralGrid; kwargs...) = OceanNeutralWindSpeed{SG.NF}(; kwargs...)

function variables(::OceanNeutralWindSpeed)
    return (
        ParameterizationVariable(:neutral_wind_speed, Grid2D(), desc = "Neutral surface wind speed", units = "m/s"),
    )
end

initialize!(::OceanNeutralWindSpeed, ::PrimitiveEquation) = nothing

@propagate_inbounds function parameterization!(ij, vars, scheme::OceanNeutralWindSpeed{NF}, model) where {NF}
    (; land_fraction) = model.land_sea_mask
    (land_fraction[ij] < 1) || return nothing
    return neutral_wind_speed(ij, vars, scheme, model)
end

"""Surface condition parameterization that calculates near-surface atmospheric
variables needed for surface flux calculations. Computes surface wind speed
including sub-grid scale gusts, surface air density, and surface air temperature
by extrapolating from the lowest model level to the surface using standard
atmospheric relationships. Fields are $(TYPEDFIELDS)"""
@kwdef struct SurfaceCondition{NF} <: AbstractSurfaceCondition
    "[OPTION] Ratio of near-surface wind to lowest-level wind [1]"
    wind_slowdown::NF = 0.95

    "[OPTION] Wind speed of sub-grid scale gusts [m/s]"
    gust_speed::NF = 5

    "[OPTION] Calculate the neutral wind speed [m/s]"
    neutral_wind::OceanNeutralWindSpeed{NF}
end

Adapt.@adapt_structure SurfaceCondition

function SurfaceCondition(
        SG::SpectralGrid;
        neutral_wind = OceanNeutralWindSpeed(SG),
        kwargs...
    )
    return SurfaceCondition{SG.NF}(; neutral_wind, kwargs...)
end

function variables(SC::AbstractSurfaceCondition)
    return (
        ParameterizationVariable(:surface_wind_speed, Grid2D(), desc = "Surface wind speed", units = "m/s"),
        ParameterizationVariable(:surface_air_density, Grid2D(), desc = "Surface air density", units = "kg/m³"),
        ParameterizationVariable(:surface_air_temperature, Grid2D(), desc = "Surface air temperature", units = "K"),
    )
end

initialize!(::SurfaceCondition, ::PrimitiveEquation) = nothing

# function barrier
@propagate_inbounds function parameterization!(ij, vars, sc::SurfaceCondition, model)
    return surface_condition!(ij, vars, sc, model)
end

@propagate_inbounds function surface_condition!(ij, vars, surface_condition::SurfaceCondition, model)

    (; wind_slowdown, gust_speed) = surface_condition
    (; nlayers) = model.geometry
    coord = model.geometry.vertical_coordinates
    (; atmosphere) = model

    # Fortran SPEEDY documentation eq. 49 but use previous time step for numerical stability
    u_grid = get_prognostic_step(vars.grid.u, model.time_stepping, surface_condition)
    v_grid = get_prognostic_step(vars.grid.v, model.time_stepping, surface_condition)
    uₛ = wind_slowdown * u_grid[ij, nlayers]
    vₛ = wind_slowdown * v_grid[ij, nlayers]

    # Fortran SPEEDY documentation eq. 50, sqrt(u² + v² + gust_speed²)
    surface_wind_speed = sqrt(muladd(uₛ, uₛ, muladd(vₛ, vₛ, gust_speed^2)))
    vars.parameterizations.surface_wind_speed[ij] = surface_wind_speed

    # Surface air density
    (; surface_air_density) = vars.parameterizations
    temperature = get_prognostic_step(vars.grid.temperature, model.time_stepping, surface_condition)
    pₛ = vars.parameterizations.surface_pressure[ij] # surface pressure [Pa]
    (; R_dry, κ) = model.atmosphere

    σ = pressure(nlayers, pₛ, coord) / pₛ           # σ vertical coordinate at lowest model level
    T = temperature[ij, nlayers]                    # virtual temperature at lowest model level [K]
    q = haskey(vars.grid, :humidity) ?              # specific humidity at lowest model level [kg/kg]
        get_prognostic_step(vars.grid.humidity, model.time_stepping, surface_condition)[ij, nlayers] : zero(T)
    Tᵥ = virtual_temperature(T, q, atmosphere)      # virtual temperature at lowest model level [K]
    σ⁻ᵏ = σ^(-κ)                                    # precalculate
    Tᵥ *= σ⁻ᵏ                                       # lower to surface assuming dry adiabatic lapse rate
    ρ = pₛ / (R_dry * Tᵥ)                           # surface air density [kg/m³] from ideal gas law
    surface_air_density[ij] = ρ                     # store for surface temp/humidity fluxes

    # Surface air temperature
    (; surface_air_temperature) = vars.parameterizations
    T *= σ⁻ᵏ                                        # lower to surface assuming dry adiabatic lapse rate
    surface_air_temperature[ij] = T                 # store for surface temp/humidity fluxes

    return nothing
end


"""Ocean-based neutral wind speed calculation from actual wind speed, 
derived from ERA5 data via symbolic regression."""
@propagate_inbounds function neutral_wind_speed(ij, vars, scheme::OceanNeutralWindSpeed{NF}, model) where {NF}
    (; surface_wind_speed) = vars.parameterizations
    (; surface_air_temperature) = vars.parameterizations

    c1 = NF(-0.039317116)
    c2 = NF(-2.9858496)
    c3 = NF(2.0046231e-10)
    c4 = NF(1.0768474)
    c5 = NF(0.20268184)
    c6 = NF(1.2684147)
    c7 = NF(-0.94933933)
    c8 = NF(0.041551278)
    c9 = NF(5.8649142)

    sst = vars.prognostic.ocean.sea_surface_temperature[ij]
    t_diff = surface_air_temperature[ij] - sst # TODO: replace SST with ocean skin temperature
    ws_safe = max(surface_wind_speed[ij], NF(1.0e-6))
    log_arg = max(c1 * ws_safe * t_diff + exp(t_diff), NF(1.0e-8))

    numerator = 2 * t_diff + c8 * exp(t_diff) - c3 * (c4^surface_air_temperature[ij])
    denominator = t_diff * (log(log_arg) + c2) + c5 * (c6^ws_safe) + c9 * (ws_safe^c7) + ws_safe

    vars.parameterizations.neutral_wind_speed[ij] = max(surface_wind_speed[ij] - (numerator / denominator), 0)

    return nothing
end
