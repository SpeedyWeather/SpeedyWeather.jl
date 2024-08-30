## LAND TEMPERATURE
abstract type AbstractLand{NF, Grid} end

function Base.show(io::IO, L::AbstractLand)
    println(io, "$(typeof(L)) <: AbstractLand")
    keys = propertynames(L)
    print_fields(io, L, keys)
end

export SeasonalLandTemperature
Base.@kwdef struct SeasonalLandTemperature{NF, Grid<:AbstractGrid{NF}} <: AbstractLand{NF, Grid}

    "number of latitudes on one hemisphere, Equator included"
    nlat_half::Int

    # OPTIONS
    "Time step used to update land surface temperatures"
    Δt::Dates.Day = Dates.Day(3)

    "path to the folder containing the land temperature file, pkg path default"
    path::String = "SpeedyWeather.jl/input_data"

    "filename of land surface temperatures"
    file::String = "land_surface_temperature.nc"

    "variable name in netcdf file"
    varname::String = "lst"

    "Grid the land surface temperature file comes on"
    file_Grid::Type{<:AbstractGrid} = FullGaussianGrid

    "The missing value in the data respresenting ocean"
    missing_value::NF = NF(NaN)

    # to be filled from file
    "Monthly land surface temperatures [K], interpolated onto Grid"
    monthly_temperature::Vector{Grid} = [zeros(Grid, nlat_half) for _ in 1:12]
end

# generator function
function SeasonalLandTemperature(SG::SpectralGrid; kwargs...)
    (; NF, Grid, nlat_half) = SG
    return SeasonalLandTemperature{NF, Grid{NF}}(; nlat_half, kwargs...)
end

function initialize!(land::SeasonalLandTemperature, model::PrimitiveEquation)
    load_monthly_climatology!(land.monthly_temperature, land)
end

function initialize!(   land::PrognosticVariablesLand,
                        time::DateTime,
                        model::PrimitiveEquation)
    land_timestep!(land, time, model.land, initialize=true)
    model isa PrimitiveWet && soil_timestep!(land, time, model.soil)
end

# function barrier
function land_timestep!(land::PrognosticVariablesLand,
                        time::DateTime,
                        model::PrimitiveEquation;
                        initialize::Bool = false)
    # the time step is dictated by the land "model" 
    executed = land_timestep!(land, time, model.land; initialize)
    executed && model isa PrimitiveWet && soil_timestep!(land, time, model.soil)
end

function land_timestep!(land::PrognosticVariablesLand{NF},
                        time::DateTime,
                        land_model::SeasonalLandTemperature;
                        initialize::Bool = false) where NF

    # escape immediately if Δt of land model hasn't passed yet
    # unless the land hasn't been initialized yet
    initialize || (time - land.time) < land_model.Δt && return false    # = executed

    # otherwise update land prognostic variables:
    land.time = time
    this_month = Dates.month(time)
    next_month = (this_month % 12) + 1      # mod for dec 12 -> jan 1

    # linear interpolation weight between the two months
    # TODO check whether this shifts the climatology by 1/2 a month
    weight = convert(NF, Dates.days(time-Dates.firstdayofmonth(time))/Dates.daysinmonth(time))

    (; monthly_temperature) = land_model
    @. land.land_surface_temperature = (1-weight) * monthly_temperature[this_month] +
                                          weight  * monthly_temperature[next_month]

    return true # = executed
end

## SOIL MOISTURE
abstract type AbstractSoil{NF, Grid} end

function Base.show(io::IO, S::AbstractSoil)
    println(io, "$(typeof(S)) <: AbstractSoil")
    keys = propertynames(S)
    print_fields(io, S, keys)
end

export SeasonalSoilMoisture
Base.@kwdef struct SeasonalSoilMoisture{NF, Grid<:AbstractGrid{NF}} <: AbstractSoil{NF, Grid}

    "number of latitudes on one hemisphere, Equator included"
    nlat_half::Int

    # OPTIONS
    "Depth of top soil layer [m]"
    D_top::Float64 = 0.07

    "Depth of root layer [m]"
    D_root::Float64 = 0.21

    "Soil wetness at field capacity [volume fraction]"
    W_cap::Float64 = 0.3

    "Soil wetness at wilting point [volume fraction]"
    W_wilt::Float64 = 0.17

    # READ CLIMATOLOGY FROM FILE
    "path to the folder containing the soil moisture file, pkg path default"
    path::String = "SpeedyWeather.jl/input_data"

    "filename of soil moisture"
    file::String = "soil_moisture.nc"

    "variable name in netcdf file"
    varname_layer1::String = "swl1"
    varname_layer2::String = "swl2"

    "Grid the soil moisture file comes on"
    file_Grid::Type{<:AbstractGrid} = FullGaussianGrid

    "The missing value in the data respresenting ocean"
    missing_value::NF = NF(NaN)

    # to be filled from file
    "Monthly soil moisture volume fraction [1], top layer, interpolated onto Grid"
    monthly_soil_moisture_layer1::Vector{Grid} = [zeros(Grid, nlat_half) for _ in 1:12]

    "Monthly soil moisture volume fraction [1], 2nd layer, interpolated onto Grid"
    monthly_soil_moisture_layer2::Vector{Grid} = [zeros(Grid, nlat_half) for _ in 1:12]
end

# generator function
function SeasonalSoilMoisture(SG::SpectralGrid; kwargs...)
    (; NF, Grid, nlat_half) = SG
    return SeasonalSoilMoisture{NF, Grid{NF}}(; nlat_half, kwargs...)
end

function initialize!(soil::SeasonalSoilMoisture, model::PrimitiveWet)
    load_monthly_climatology!(soil.monthly_soil_moisture_layer1, soil, varname=soil.varname_layer1)
    load_monthly_climatology!(soil.monthly_soil_moisture_layer2, soil, varname=soil.varname_layer2)
end

function soil_timestep!(land::PrognosticVariablesLand{NF},
                        time::DateTime,
                        soil_model::SeasonalSoilMoisture) where NF

    this_month = Dates.month(time)
    next_month = (this_month % 12) + 1      # mod for dec 12 -> jan 1

    # linear interpolation weight between the two months
    # TODO check whether this shifts the climatology by 1/2 a month
    weight = convert(NF, Dates.days(time-Dates.firstdayofmonth(time))/Dates.daysinmonth(time))

    (; monthly_soil_moisture_layer1, monthly_soil_moisture_layer2) = soil_model
    @. land.soil_moisture_layer1 = (1-weight) * monthly_soil_moisture_layer1[this_month] +
                                        weight  * monthly_soil_moisture_layer1[next_month]
    @. land.soil_moisture_layer2 = (1-weight) * monthly_soil_moisture_layer2[this_month] +
                                        weight  * monthly_soil_moisture_layer2[next_month]

    return nothing
end

## SOIL MOISTURE
abstract type AbstractVegetation{NF, Grid} <: AbstractParameterization end

export VegetationClimatology
Base.@kwdef struct VegetationClimatology{NF, Grid<:AbstractGrid{NF}} <: AbstractVegetation{NF, Grid}

    "number of latitudes on one hemisphere, Equator included"
    nlat_half::Int

    # OPTIONS
    "Combine high and low vegetation factor, a in high + a*low [1]"
    low_veg_factor::Float64 = 0.8

    "path to the folder containing the soil moisture file, pkg path default"
    path::String = "SpeedyWeather.jl/input_data"

    "filename of soil moisture"
    file::String = "vegetation.nc"

    "variable name in netcdf file"
    varname_vegh::String = "vegh"
    varname_vegl::String = "vegl"

    "Grid the soil moisture file comes on"
    file_Grid::Type{<:AbstractGrid} = FullGaussianGrid

    "The missing value in the data respresenting ocean"
    missing_value::NF = NF(NaN)

    # to be filled from file
    "High vegetation cover [1], interpolated onto Grid"
    high_cover::Grid = zeros(Grid, nlat_half)

    "Low vegetation cover [1], interpolated onto Grid"
    low_cover::Grid = zeros(Grid, nlat_half)
end

# generator function
function VegetationClimatology(SG::SpectralGrid; kwargs...)
    (; NF, Grid, nlat_half) = SG
    return VegetationClimatology{NF, Grid{NF}}(; nlat_half, kwargs...)
end

function initialize!(vegetation::VegetationClimatology, model::PrimitiveEquation)

    # LOAD NETCDF FILE
    if vegetation.path == "SpeedyWeather.jl/input_data"
        path = joinpath(@__DIR__, "../../input_data", vegetation.file)
    else
        path = joinpath(vegetation.path, vegetation.file)
    end
    ncfile = NCDataset(path)

    # high and low vegetation cover
    vegh = vegetation.file_Grid(ncfile[vegetation.varname_vegh].var[:, :], input_as=Matrix)
    vegl = vegetation.file_Grid(ncfile[vegetation.varname_vegl].var[:, :], input_as=Matrix)

    # interpolate onto grid
    high_vegetation_cover = vegetation.high_cover
    low_vegetation_cover = vegetation.low_cover
    interpolator = RingGrids.interpolator(Float32, high_vegetation_cover, vegh)
    interpolate!(high_vegetation_cover, vegh, interpolator)
    interpolate!(low_vegetation_cover, vegl, interpolator)
end

# function barrier
function soil_moisture_availability!(
    diagn::SurfaceVariables,
    land::PrognosticVariablesLand,
    model::PrimitiveWet)
    soil_moisture_availability!(diagn, land, model.soil, model.vegetation)
end

function soil_moisture_availability!(
    ::SurfaceVariables,
    ::PrognosticVariablesLand,
    ::PrimitiveDry)
    return nothing
end

function soil_moisture_availability!(
    diagn::SurfaceVariables{NF},
    land::PrognosticVariablesLand,
    soil::AbstractSoil,
    vegetation::AbstractVegetation) where NF

    (; soil_moisture_availability) = diagn
    (; soil_moisture_layer1, soil_moisture_layer2) = land
    (; high_cover, low_cover) = vegetation
    
    D_top = convert(NF, soil.D_top)
    D_root = convert(NF, soil.D_root)
    W_cap = convert(NF, soil.W_cap)
    W_wilt = convert(NF, soil.W_wilt)
    low_veg_factor = convert(NF, vegetation.low_veg_factor)

    # precalculate
    r = 1/(D_top*W_cap + D_root*(W_cap - W_wilt))

    @inbounds for ij in eachgridpoint(soil_moisture_availability,
                            soil_moisture_layer1,
                            soil_moisture_layer2,
                            high_cover, low_cover)
        
        # Fortran SPEEDY source/land_model.f90 line 111 origin unclear
        veg = max(0, high_cover[ij] + low_veg_factor*low_cover[ij])

        # Fortran SPEEDY documentation eq. 51
        soil_moisture_availability[ij] = r*(D_top*soil_moisture_layer1[ij] +
            veg*D_root*max(soil_moisture_layer2[ij]- W_wilt, 0))
    end
end
