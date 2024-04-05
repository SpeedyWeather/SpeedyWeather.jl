abstract type AbstractSurfaceThermodynamics <: AbstractParameterization end
abstract type AbstractSurfaceWind <: AbstractParameterization end
abstract type AbstractSurfaceHeatFlux <: AbstractParameterization end
abstract type AbstractSurfaceEvaporation <: AbstractParameterization end

# defines the order in which they are called und unpacks to dispatch
function surface_fluxes!(column::ColumnVariables, model::PrimitiveEquation)

    # get temperature, humidity and density at surface
    surface_thermodynamics!(column, model.surface_thermodynamics, model)

    # also calculates surface wind speed necessary for other fluxes too
    surface_wind_stress!(column, model.surface_wind, model)

    # now call other heat (wet and dry) and humidity fluxes (PrimitiveWet only)
    surface_heat_flux!(column, model.surface_heat_flux, model)
    model isa PrimitiveWet && surface_evaporation!(column, model.surface_evaporation, model)
end

## SURFACE THERMODYNAMICS
export SurfaceThermodynamicsConstant
struct SurfaceThermodynamicsConstant <: AbstractSurfaceThermodynamics end
SurfaceThermodynamicsConstant(SG::SpectralGrid) = SurfaceThermodynamicsConstant()
initialize!(::SurfaceThermodynamicsConstant,::PrimitiveEquation) = nothing

function surface_thermodynamics!(   column::ColumnVariables,
                                    ::SurfaceThermodynamicsConstant,
                                    model::PrimitiveWet)

    # surface value is same as lowest model level
    column.surface_temp = column.temp[end]      # todo use constant POTENTIAL temperature
    column.surface_humid = column.humid[end]    # humidity at surface is the same as 

    # surface air density via virtual temperature
    (; R_dry) = model.atmosphere
    Tᵥ = column.temp_virt[column.nlev]
    column.surface_air_density = column.pres[end]/(R_dry*Tᵥ)
end

function surface_thermodynamics!(   column::ColumnVariables,
                                    ::SurfaceThermodynamicsConstant,
                                    model::PrimitiveDry)
    (; R_dry) = model.atmosphere
    # surface value is same as lowest model level
    column.surface_temp = column.temp[end]   # todo use constant POTENTIAL temperature
    column.surface_air_density = column.pres[end]/(R_dry*column.surface_temp)
end

## WIND STRESS
export NoSurfaceWind
struct NoSurfaceWind <: AbstractSurfaceWind end
NoSurfaceWind(::SpectralGrid) = NoSurfaceWind()
initialize!(::NoSurfaceWind, ::PrimitiveEquation) = nothing
surface_wind_stress!(::ColumnVariables, ::NoSurfaceWind, ::PrimitiveEquation) = nothing

export SurfaceWind
Base.@kwdef struct SurfaceWind{NF<:AbstractFloat} <: AbstractSurfaceWind
    "Ratio of near-surface wind to lowest-level wind [1]"
    f_wind::NF = 0.95

    "Wind speed of sub-grid scale gusts [m/s]"
    V_gust::NF = 5

    "Use (possibly) flow-dependent column.boundary_layer_drag coefficient"
    use_boundary_layer_drag::Bool = true

    "Otherwise, drag coefficient over land (orography = 0) [1]"
    drag_land::NF = 2.4e-3
    
    "Otherwise, Drag coefficient over sea [1]"
    drag_sea::NF = 1.8e-3

    "Flux limiter to cap the max of surface momentum fluxes [kg/m/s²]"
    max_flux::NF = 0.1
end

SurfaceWind(SG::SpectralGrid; kwargs...) = SurfaceWind{SG.NF}(; kwargs...)
initialize!(::SurfaceWind, ::PrimitiveEquation) = nothing

function surface_wind_stress!(  column::ColumnVariables,
                                surface_wind::SurfaceWind,
                                model::PrimitiveEquation)

    (; land_fraction) = column
    (; f_wind, V_gust, drag_land, drag_sea) = surface_wind
    (; max_flux) = surface_wind

    # SPEEDY documentation eq. 49
    column.surface_u = f_wind*column.u[end] 
    column.surface_v = f_wind*column.v[end]
    (; surface_u, surface_v) = column

    # SPEEDY documentation eq. 50
    column.surface_wind_speed = sqrt(surface_u^2 + surface_v^2 + V_gust^2)

    # drag coefficient either from SurfaceEvaporation or from a central drag coefficient
    drag_sea, drag_land = surface_wind.use_boundary_layer_drag ?
                                (column.boundary_layer_drag, column.boundary_layer_drag) : 
                                (drag_sea, drag_land)
    
    # surface wind stress: quadratic drag, fractional land-sea mask
    ρ = column.surface_air_density
    V₀ = column.surface_wind_speed
    drag = land_fraction*drag_land + (1-land_fraction)*drag_sea

    # add flux limiter to avoid heavy drag in (initial) shock
    u_flux = ρ*drag*V₀*surface_u
    v_flux = ρ*drag*V₀*surface_v
    column.flux_u_upward[end] -= clamp(u_flux, -max_flux, max_flux)
    column.flux_v_upward[end] -= clamp(v_flux, -max_flux, max_flux)

    # SPEEDY documentation eq. 52, 53, accumulate fluxes with +=
    # column.flux_u_upward[end] -= ρ*drag*V₀*surface_u
    # column.flux_v_upward[end] -= ρ*drag*V₀*surface_v
    
    return nothing
end

## SENSIBLE HEAT FLUX
export NoSurfaceHeatFlux
struct NoSurfaceHeatFlux <: AbstractSurfaceHeatFlux end
NoSurfaceHeatFlux(::SpectralGrid) = NoSurfaceHeatFlux()
initialize!(::NoSurfaceHeatFlux, ::PrimitiveEquation) = nothing
surface_heat_flux!(::ColumnVariables, ::NoSurfaceHeatFlux, ::PrimitiveEquation) = nothing

export SurfaceHeatFlux
Base.@kwdef struct SurfaceHeatFlux{NF<:AbstractFloat} <: AbstractSurfaceHeatFlux
    
    "Use (possibly) flow-dependent column.boundary_layer_drag coefficient"
    use_boundary_layer_drag::Bool = true
    
    "Otherwise, use the following drag coefficient for heat fluxes over land"
    heat_exchange_land::NF = 1.2e-3    # for neutral stability

    "Otherwise, use the following drag coefficient for heat fluxes over sea"
    heat_exchange_sea::NF = 0.9e-3

    "Flux limiter for surface heat fluxes [W/m²]"
    max_flux::NF = 100
end

SurfaceHeatFlux(SG::SpectralGrid; kwargs...) = SurfaceHeatFlux{SG.NF}(; kwargs...)
initialize!(::SurfaceHeatFlux, ::PrimitiveEquation) = nothing

function surface_heat_flux!(   
    column::ColumnVariables,
    heat_flux::SurfaceHeatFlux,
    model::PrimitiveEquation,
)   
    cₚ = model.atmosphere.heat_capacity
    (; heat_exchange_land, heat_exchange_sea, max_flux) = heat_flux

    ρ = column.surface_air_density
    V₀ = column.surface_wind_speed
    T_skin_sea = column.skin_temperature_sea
    T_skin_land = column.skin_temperature_land
    T = column.surface_temp
    land_fraction = column.land_fraction

    # drag coefficient
    drag_sea, drag_land = heat_flux.use_boundary_layer_drag ?
                        (column.boundary_layer_drag, column.boundary_layer_drag) : 
                        (heat_exchange_sea, heat_exchange_land)

    # SPEEDY documentation Eq. 54 and 56
    flux_land = ρ*drag_land*V₀*cₚ*(T_skin_land - T)
    flux_sea  = ρ*drag_sea*V₀*cₚ*(T_skin_sea  - T)

    # flux limiter
    flux_land = clamp(flux_land, -max_flux, max_flux)
    flux_sea = clamp(flux_sea, -max_flux, max_flux)

    # mix fluxes for fractional land-sea mask
    land_available = isfinite(T_skin_land)
    sea_available = isfinite(T_skin_sea)

    if land_available && sea_available
        column.flux_temp_upward[end] += land_fraction*flux_land + (1-land_fraction)*flux_sea

    # but in case only land or sea are available use those ones only
    elseif land_available
        column.flux_temp_upward[end] += flux_land

    elseif sea_available
        column.flux_temp_upward[end] += flux_sea

    # or no flux in case none is defined (shouldn't happen with default surface data)
    # else   # doesn't have to be executed because fluxes are accumulated
    #   column.flux_temp_upward[end] += 0
    end
    return nothing
end

## SURFACE EVAPORATION
export NoSurfaceEvaporation
struct NoSurfaceEvaporation <: AbstractSurfaceHeatFlux end
NoSurfaceEvaporation(::SpectralGrid) = NoSurfaceEvaporation()
initialize!(::NoSurfaceEvaporation, ::PrimitiveEquation) = nothing
surface_evaporation!(::ColumnVariables, ::NoSurfaceEvaporation, ::PrimitiveEquation) = nothing

export SurfaceEvaporation

"""
Surface evaporation following a bulk formula with wind from model.surface_wind 
$(TYPEDFIELDS)"""
Base.@kwdef struct SurfaceEvaporation{NF<:AbstractFloat} <: AbstractSurfaceEvaporation
    
    "Use column.boundary_layer_drag coefficient"
    use_boundary_layer_drag::Bool = true

    "Otherwise, use the following drag coefficient for evaporation over land"
    moisture_exchange_land::NF = 1.2e-3

    "Or drag coefficient for evaporation over sea"
    moisture_exchange_sea::NF = 0.9e-3
end

SurfaceEvaporation(SG::SpectralGrid; kwargs...) = SurfaceEvaporation{SG.NF}(; kwargs...)
initialize!(::SurfaceEvaporation, ::PrimitiveWet) = nothing

function surface_evaporation!(  column::ColumnVariables,
                                evaporation::SurfaceEvaporation,
                                model::PrimitiveWet)

    (; skin_temperature_sea, skin_temperature_land, pres) = column
    (; moisture_exchange_land, moisture_exchange_sea) = evaporation
    α = column.soil_moisture_availability

    # SATURATION HUMIDITY OVER LAND AND OCEAN
    surface_pressure = pres[end]
    sat_humid_land = saturation_humidity(skin_temperature_land, surface_pressure, model.clausius_clapeyron)
    sat_humid_sea = saturation_humidity(skin_temperature_sea, surface_pressure, model.clausius_clapeyron)

    ρ = column.surface_air_density
    V₀ = column.surface_wind_speed
    land_fraction = column.land_fraction
    (; surface_humid) = column

    # drag coefficient either from SurfaceEvaporation or from a central drag coefficient
    drag_sea, drag_land = evaporation.use_boundary_layer_drag ?
                        (column.boundary_layer_drag, column.boundary_layer_drag) : 
                        (moisture_exchange_sea, moisture_exchange_land)

    flux_sea = ρ*drag_sea*V₀*max(sat_humid_sea  - surface_humid, 0)
    flux_land = ρ*drag_land*V₀*α*max(sat_humid_land  - surface_humid, 0)

    # mix fluxes for fractional land-sea mask
    land_available = isfinite(skin_temperature_land) && isfinite(α)
    sea_available = isfinite(skin_temperature_sea)

    if land_available && sea_available
        column.flux_humid_upward[end] += land_fraction*flux_land + (1-land_fraction)*flux_sea

    # but in case only land or sea are available use those ones only
    elseif land_available
        column.flux_humid_upward[end] += flux_land

    elseif sea_available
        column.flux_humid_upward[end] += flux_sea

    # or no flux in case none is defined (shouldn't happen with default surface data)
    # else   # doesn't have to be executed because fluxes are accumulated
    #   column.flux_temp_upward[end] += 0
    end
end