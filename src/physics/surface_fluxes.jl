abstract type AbstractSurfaceWind <: AbstractParameterization end
abstract type AbstractSurfaceThermodynamics <: AbstractParameterization end
abstract type AbstractSurfaceHeat <: AbstractParameterization end
abstract type AbstractEvaporation <: AbstractParameterization end

function surface_fluxes!(column::ColumnVariables, model::PrimitiveEquation)

    # get temperature, humidity and density at surface
    surface_thermodynamics!(column, model.surface_thermodynamics, model.atmosphere, model)

    # also calculates surface wind speed necessary for other fluxes too
    surface_wind_stress!(column, model.surface_wind)

    # now call other heat and humidity fluxes
    sensible_heat_flux!(column, model.surface_heat_flux, model.atmosphere)
    evaporation!(column, model)
end

export SurfaceThermodynamicsConstant
struct SurfaceThermodynamicsConstant{NF<:AbstractFloat} <: AbstractSurfaceThermodynamics end
SurfaceThermodynamicsConstant(SG::SpectralGrid) = SurfaceThermodynamicsConstant{SG.NF}()

function surface_thermodynamics!(   column::ColumnVariables,
                                    ::SurfaceThermodynamicsConstant,
                                    atmosphere::AbstractAtmosphere,
                                    model::PrimitiveWet)

    # surface value is same as lowest model level
    column.surface_temp = column.temp[end]      # todo use constant POTENTIAL temperature
    column.surface_humid = column.humid[end]    # humidity at surface is the same as 

    # surface air density via virtual temperature
    (; R_dry) = atmosphere
    Tᵥ = column.temp_virt[column.nlev]
    column.surface_air_density = column.pres[end]/(R_dry*Tᵥ)
end

function surface_thermodynamics!(   column::ColumnVariables,
                                    ::SurfaceThermodynamicsConstant,
                                    atmosphere::AbstractAtmosphere,
                                    model::PrimitiveDry)
    (; R_dry) = atmosphere
    # surface value is same as lowest model level
    column.surface_temp = column.temp[end]   # todo use constant POTENTIAL temperature
    column.surface_air_density = column.pres[end]/(R_dry*column.surface_temp)
end

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

function surface_wind_stress!(  column::ColumnVariables{NF},
                                surface_wind::SurfaceWind) where NF

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

export SurfaceSensibleHeat
Base.@kwdef struct SurfaceSensibleHeat{NF<:AbstractFloat} <: AbstractSurfaceHeat
    
    "Use (possibly) flow-dependent column.boundary_layer_drag coefficient"
    use_boundary_layer_drag::Bool = true
    
    "Otherwise, use the following drag coefficient for heat fluxes over land"
    heat_exchange_land::NF = 1.2e-3    # for neutral stability

    "Otherwise, use the following drag coefficient for heat fluxes over sea"
    heat_exchange_sea::NF = 0.9e-3

    "Flux limiter for surface heat fluxes [W/m²]"
    max_flux::NF = 100
end

SurfaceSensibleHeat(SG::SpectralGrid; kwargs...) = SurfaceSensibleHeat{SG.NF}(; kwargs...)

function sensible_heat_flux!(   
    column::ColumnVariables,
    heat_flux::SurfaceSensibleHeat,
    atmosphere::AbstractAtmosphere
)   
    cₚ = atmosphere.heat_capacity
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

export SurfaceEvaporation

"""
Surface evaporation following a bulk formula with wind from model.surface_wind 
$(TYPEDFIELDS)"""
Base.@kwdef struct SurfaceEvaporation{NF<:AbstractFloat} <: AbstractEvaporation
    
    "Use column.boundary_layer_drag coefficient"
    use_boundary_layer_drag::Bool = true

    "Otherwise, use the following drag coefficient for evaporation over land"
    moisture_exchange_land::NF = 1.2e-3

    "Or drag coefficient for evaporation over sea"
    moisture_exchange_sea::NF = 0.9e-3
end

SurfaceEvaporation(SG::SpectralGrid; kwargs...) = SurfaceEvaporation{SG.NF}(; kwargs...)

# don't do anything for dry core
function evaporation!(  column::ColumnVariables,
                        model::PrimitiveDry)
    return nothing
end

# function barrier
function evaporation!(  column::ColumnVariables,
                        model::PrimitiveWet)
    evaporation!(column, model.evaporation, model.clausius_clapeyron)
end

function evaporation!(  column::ColumnVariables{NF},
                        evaporation::SurfaceEvaporation,
                        clausius_clapeyron::AbstractClausiusClapeyron) where NF

    (; skin_temperature_sea, skin_temperature_land, pres) = column
    (; moisture_exchange_land, moisture_exchange_sea) = evaporation
    α = column.soil_moisture_availability

    # SATURATION HUMIDITY OVER LAND AND OCEAN
    surface_pressure = pres[end]
    sat_humid_land = saturation_humidity(skin_temperature_land, surface_pressure, clausius_clapeyron)
    sat_humid_sea = saturation_humidity(skin_temperature_sea, surface_pressure, clausius_clapeyron)

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