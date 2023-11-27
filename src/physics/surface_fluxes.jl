function surface_fluxes!(column::ColumnVariables,model::PrimitiveEquation)

    # get temperature, humidity and density at surface
    surface_thermodynamics!(column,model.surface_thermodynamics,model.constants,model)

    # also calculates surface wind speed necessary for other fluxes too
    surface_wind_stress!(column,model.surface_wind)

    # now call other heat and humidity fluxes
    sensible_heat_flux!(column,model.surface_heat_flux,model.constants)
    evaporation!(column,model)
end

struct SurfaceThermodynamicsConstant{NF<:AbstractFloat} <: AbstractSurfaceThermodynamics{NF} end
struct SurfaceThermodynamicsExtrapolate{NF<:AbstractFloat} <: AbstractSurfaceThermodynamics{NF} end

SurfaceThermodynamicsConstant(SG::SpectralGrid) = SurfaceThermodynamicsConstant{SG.NF}()
SurfaceThermodynamicsExtrapolate(SG::SpectralGrid) = SurfaceThermodynamicsExtrapolate{SG.NF}()

function surface_thermodynamics!(   column::ColumnVariables,
                                    ::SurfaceThermodynamicsConstant,
                                    C::DynamicsConstants,
                                    model::PrimitiveWet)

    # surface value is same as lowest model level
    column.temp[end] = column.temp[end-1]       # todo use constant POTENTIAL temperature
    column.humid[end] = column.humid[end-1]     # humidity at surface is the same as 

    # surface air density via virtual temperature
    (;R_dry,μ_virt_temp) = C
    T = column.temp[end]        # surface air temperature
    q = column.humid[end]       # surface humidity
    Tᵥ = T*(1 + μ_virt_temp*q)  # virtual temperature
    column.surface_air_density = column.pres[end]/R_dry/Tᵥ
end

function surface_thermodynamics!(   column::ColumnVariables,
                                    ::SurfaceThermodynamicsConstant,
                                    C::DynamicsConstants,
                                    model::PrimitiveDry)

    # surface value is same as lowest model level
    column.temp[end] = column.temp[end-1]       # todo use constant POTENTIAL temperature
    column.surface_air_density = column.pres[end]/C.R_dry/column.temp[end]
end

Base.@kwdef struct SurfaceWind{NF<:AbstractFloat} <: AbstractSurfaceWind{NF}
    "Ratio of near-surface wind to lowest-level wind [1]"
    f_wind::NF = 0.95

    "Wind speed of sub-grid scale gusts [m/s]"
    V_gust::NF = 5

    # TODO make this orography dependent
    "Drag coefficient over land (orography = 0) [1]"
    drag_land::NF = 2.4e-3

    "Drag coefficient over sea [1]"
    drag_sea::NF = 1.8e-3

    "Flux limiter [kg/m/s²]"
    max_flux::NF = 1.5
end

SurfaceWind(SG::SpectralGrid;kwargs...) = SurfaceWind{SG.NF}(;kwargs...)

function surface_wind_stress!(  column::ColumnVariables{NF},
                                surface_wind::SurfaceWind) where NF

    (;u,v,land_fraction) = column
    (;f_wind, V_gust, drag_land, drag_sea, max_flux) = surface_wind     

    # SPEEDY documentation eq. 49
    u[end] = f_wind*u[end-1] 
    v[end] = f_wind*v[end-1]

    # SPEEDY documentation eq. 50
    column.surface_wind_speed = sqrt(u[end]^2 + v[end]^2 + V_gust^2)

    # surface wind stress: quadratic drag, fractional land-sea mask
    ρ = column.surface_air_density
    V₀ = column.surface_wind_speed
    drag = land_fraction*drag_land + (1-land_fraction)*drag_sea

    # SPEEDY documentation eq. 52, 53, accumulate fluxes with +=
    # add flux limiter to avoid heavy drag in initial shock
    # u_flux = sign(u[end])*min(abs(ρdragV₀*u[end]),max_flux)
    # v_flux = sign(v[end])*min(abs(ρdragV₀*v[end]),max_flux)
    # column.flux_u_upward[end] -= u_flux
    # column.flux_v_upward[end] -= v_flux

    column.flux_u_upward[end] -= ρ*drag*V₀*u[end]
    column.flux_v_upward[end] -= ρ*drag*V₀*v[end]
    
    return nothing
end

Base.@kwdef struct SurfaceSensibleHeat{NF<:AbstractFloat} <: AbstractSurfaceHeat{NF}
    heat_exchange_land::Float64 = 1.2e-3    # for neutral stability
    heat_exchange_sea::Float64 = 0.9e-3
end

SurfaceSensibleHeat(SG::SpectralGrid;kwargs...) = SurfaceSensibleHeat{SG.NF}(;kwargs...)

function sensible_heat_flux!(   column::ColumnVariables{NF},
                                heat_flux::SurfaceSensibleHeat,
                                C::DynamicsConstants) where NF
    (;cₚ) = C
    heat_exchange_land = convert(NF,heat_flux.heat_exchange_land)
    heat_exchange_sea  = convert(NF,heat_flux.heat_exchange_sea)

    ρ = column.surface_air_density
    V₀ = column.surface_wind_speed
    T_skin_sea = column.skin_temperature_sea
    T_skin_land = column.skin_temperature_land
    T = column.temp[end]
    land_fraction = column.land_fraction

    # SPEEDY documentation Eq. 54 and 56
    flux_land = ρ*heat_exchange_land*V₀*cₚ*(T_skin_land - T)
    flux_sea  = ρ*heat_exchange_sea*V₀*cₚ*(T_skin_sea  - T)

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

Base.@kwdef struct SurfaceEvaporation{NF<:AbstractFloat} <: AbstractEvaporation{NF}
    moisture_exchange_land::Float64 = 1.2e-3    # for neutral stability
    moisture_exchange_sea::Float64 = 0.9e-3
end

SurfaceEvaporation(SG::SpectralGrid;kwargs...) = SurfaceEvaporation{SG.NF}(;kwargs...)

# don't do anything for dry core
function evaporation!(  column::ColumnVariables,
                        model::PrimitiveDry)
    return nothing
end

# function barrier
function evaporation!(  column::ColumnVariables,
                        model::PrimitiveWet)
    evaporation!(column,model.evaporation,model.thermodynamics)
end

function evaporation!(  column::ColumnVariables{NF},
                        evaporation::SurfaceEvaporation,
                        thermodynamics::Thermodynamics) where NF

    (;skin_temperature_sea, skin_temperature_land, pres, humid) = column
    moisture_exchange_land = convert(NF,evaporation.moisture_exchange_land)
    moisture_exchange_sea  = convert(NF,evaporation.moisture_exchange_sea)
    α = column.soil_moisture_availability

    # SATURATION HUMIDITY OVER LAND AND OCEAN
    surface_pressure = pres[end]
    sat_humid_land = saturation_humidity(skin_temperature_land,surface_pressure,thermodynamics)
    sat_humid_sea = saturation_humidity(skin_temperature_sea,surface_pressure,thermodynamics)

    ρ = column.surface_air_density
    V₀ = column.surface_wind_speed
    land_fraction = column.land_fraction
    flux_sea = ρ*moisture_exchange_sea*V₀*max(sat_humid_sea  - humid[end],0)
    flux_land = ρ*moisture_exchange_land*V₀*α*max(sat_humid_land  - humid[end],0)

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