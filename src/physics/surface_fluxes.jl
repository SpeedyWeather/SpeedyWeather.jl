function surface_fluxes!(column::ColumnVariables,model::PrimitiveEquation)

    # get temperature, humidity and density at surface
    surface_thermodynamics!(column,model.surface_thermodynamics,model.constants,model)

    # also calculates surface wind speed necessary for other fluxes too
    surface_wind_stress!(column,model.surface_wind)

    # now call other heat and humidity fluxes
    # soil_moisture!(column,model.vegetation)
    sensible_heat_flux!(column,model.surface_heat_flux,model.constants)
    # evaporation!(column,model.surface_evaporation,model)
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
    f_wind::Float64 = 0.95

    "Wind speed of sub-grid scale gusts [m/s]"
    V_gust::Float64 = 5

    # TODO make this orography dependent
    "Drag coefficient over land (orography = 0) [1]"
    drag_land::Float64 = 2.4e-3

    "Drag coefficient over sea [1]"
    drag_sea::Float64 = 1.8e-3
end

SurfaceWind(SG::SpectralGrid;kwargs...) = SurfaceWind{SG.NF}(;kwargs...)

function surface_wind_stress!(  column::ColumnVariables{NF},
                                surface_wind::SurfaceWind) where NF

    (;u,v,land_fraction) = column
    f_wind = convert(NF,surface_wind.f_wind)
    V_gust = convert(NF,surface_wind.V_gust)            
    drag_land = convert(NF,surface_wind.drag_land)            
    drag_sea = convert(NF,surface_wind.drag_sea)            

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
    column.flux_u_upward[end] += -ρ*drag*V₀*u[end]
    column.flux_v_upward[end] += -ρ*drag*V₀*v[end]

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
    # flux_land = ρ*heat_exchange_land*V₀*cₚ*(T_skin_land - T)
    flux_land = 0
    flux_sea  = ρ*heat_exchange_sea *V₀*cₚ*(T_skin_sea  - T)

    # mix fluxes for fractional land-sea mask
    land_available = ~isnan(T_skin_land)
    sea_available = ~isnan(T_skin_sea)

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

    nlat_half::Int

    moisture_exchange_land::Float64 = 1.2e-3    # for neutral stability
    moisture_exchange_sea::Float64 = 0.9e-3

    "Depth of top soil layer [m]"
    D_top::Float64 = 0.07

    "Depth of root layer [m]"
    D_root::Float64 = 0.21

    "Soil wetness at field capacity [volume fraction]"
    W_cap::Float64 = 0.3

    "Soil wetness at wilting point [volume fraction]"
    W_wilt::Float64 = 0.17

    # FILE OPTIONS
    "path to the folder containing the soil moisture and vegetation file, pkg path default"
    path::String = "SpeedyWeather.jl/input_data"

    "filename of that netcdf file"
    file::String = "soil_moisture.nc"

    # "Grid the soil and vegetation file comes on"
    # file_Grid::Type{<:AbstractGrid} = FullGaussianGrid

    # # arrays to be initialised
    # "Soil moisture availability index"
    # soil_moisture_availability::Grid = zeros(Grid{NF},nlat_half)
end

function initialize!(evaporation::SurfaceEvaporation)

    # LOAD NETCDF FILE
    if evaporation.path == "SpeedyWeather.jl/input_data"
        path = joinpath(@__DIR__,"../../input_data",evaporation.file)
    else
        path = joinpath(evaporation.path,evaporation.file)
    end
    ncfile = NCDataset(path)

    # high resolution land-sea mask
    lsm_highres = evaporation.file_Grid(ncfile["lsm"][:,:])

    # average onto grid cells of the model
    RingGrids.grid_cell_average!(land_sea_mask.land_sea_mask,lsm_highres)
    end


function SurfaceEvaporation(SG::SpectralGrid;kwargs...)
    return SurfaceEvaporation(nlat_half=SG.nlat_half;kwargs...)
end

# don't do anything for dry core
function evaporation!(  column::ColumnVariables,
                        moisture_flux::AbstractEvaporation,
                        model::PrimitiveDry)
    return nothing
end

# function barrier
function evaporation!(  column::ColumnVariables{NF},
                        moisture_flux::SurfaceEvaporation,
                        model::PrimitiveWet) where NF
    evaporation!(column,moisture_flux,model.soil_moisture,model.thermodynamics)
end

function evaporation!(  column::ColumnVariables{NF},
                        evaporation::SurfaceEvaporation,
                        thermodynamics::Thermodynamics) where NF

    (;skin_temperature_sea, skin_temperature_land, pres, humid) = column
    (;e₀, T₀, C₁, C₂, T₁, T₂) = thermodynamics.magnus_coefs
    (;mol_ratio) = thermodynamics      # = mol_mass_vapour/mol_mass_dry_air = 0.622
    α = soil_moisture.availability[column.ij]

    # LAND SKIN TEMPERATURE
    # change coefficients for water (temp > T₀) or ice (else)
    C, T = skin_temperature_land > T₀ ? (C₁, T₁) : (C₂, T₂)
    sat_vap_pres = e₀ * exp(C * (skin_temperature_land - T₀) / (skin_temperature_land - T))
    sat_humid_land = mol_ratio*sat_vap_pres / (pres[end] - (1-mol_ratio)*sat_vap_pres)

    # SEA SKIN TEMPERATURE
    C, T = skin_temperature_sea > T₀ ? (C₁, T₁) : (C₂, T₂)
    sat_vap_pres = e₀ * exp(C * (skin_temperature_sea - T₀) / (skin_temperature_sea - T))
    sat_humid_sea = mol_ratio*sat_vap_pres / (pres[end] - (1-mol_ratio)*sat_vap_pres)

end