abstract type AbstractOcean end

function Base.show(io::IO, O::AbstractOcean)
    println(io, "$(typeof(O)) <: AbstractOcean")
    keys = propertynames(O)
    print_fields(io, O, keys)
end

# function barrier for all oceans
function initialize!(   ocean::PrognosticVariablesOcean,
                        time::DateTime,
                        model::PrimitiveEquation)
    initialize!(ocean, time, model.ocean, model)
end

# function barrier for all oceans
function ocean_timestep!(   ocean::PrognosticVariablesOcean,
                            time::DateTime,
                            model::PrimitiveEquation)
    ocean_timestep!(ocean, time, model.ocean)
end


## SEASONAL OCEAN CLIMATOLOGY

export SeasonalOceanClimatology
Base.@kwdef struct SeasonalOceanClimatology{NF, Grid<:AbstractGrid{NF}} <: AbstractOcean

    "number of latitudes on one hemisphere, Equator included"
    nlat_half::Int

    "[OPTION] Time step used to update sea surface temperatures"
    Δt::Dates.Day = Dates.Day(3)

    "[OPTION] Path to the folder containing the sea surface temperatures, pkg path default"
    path::String = "SpeedyWeather.jl/input_data"

    "[OPTION] Filename of sea surface temperatures"
    file::String = "sea_surface_temperature.nc"

    "[OPTION] Variable name in netcdf file"
    varname::String = "sst"

    "[OPTION] Grid the sea surface temperature file comes on"
    file_Grid::Type{<:AbstractGrid} = FullGaussianGrid

    "[OPTION] The missing value in the data respresenting land"
    missing_value::NF = NF(NaN)

    # to be filled from file
    "Monthly sea surface temperatures [K], interpolated onto Grid"
    monthly_temperature::Vector{Grid} = [zeros(Grid, nlat_half) for _ in 1:12]
end

# generator function
function SeasonalOceanClimatology(SG::SpectralGrid; kwargs...)
    (; NF, Grid, nlat_half) = SG
    return SeasonalOceanClimatology{NF, Grid{NF}}(; nlat_half, kwargs...)
end

function initialize!(ocean::SeasonalOceanClimatology, model::PrimitiveEquation)
    load_monthly_climatology!(ocean.monthly_temperature, ocean)
end

function load_monthly_climatology!( 
    monthly::Vector{Grid},
    scheme;
    varname::String = scheme.varname
) where {Grid<:AbstractGrid}

    # LOAD NETCDF FILE
    if scheme.path == "SpeedyWeather.jl/input_data"
        path = joinpath(@__DIR__, "../../input_data", scheme.file)
    else
        path = joinpath(scheme.path, scheme.file)
    end
    ncfile = NCDataset(path)

    # create interpolator from grid in file to grid used in model
    nx, ny = ncfile.dim["lon"], ncfile.dim["lat"]
    npoints = nx*ny
    NF_file = typeof(ncfile[varname].attrib["_FillValue"])
    grid = scheme.file_Grid(zeros(NF_file, npoints))
    interp = RingGrids.interpolator(Float32, monthly[1], grid)

    # interpolate and store in monthly
    for month in 1:12
        this_month = ncfile[varname][:, :, month]
        ij = 0
        for j in 1:ny
            for i in 1:nx
                ij += 1
                x = this_month[i, j]
                grid[ij] = ismissing(x) ? scheme.missing_value : x
            end
        end
        interpolate!(monthly[month], grid, interp)
    end
    return nothing
end

function initialize!(   
    ocean::PrognosticVariablesOcean,
    time::DateTime,
    ocean_model::SeasonalOceanClimatology,
    model::PrimitiveEquation,
)
    ocean.time = time   # set initial time
    interpolate_monthly!(   ocean.sea_surface_temperature,
                            ocean_model.monthly_temperature,
                            time)
end
    
function ocean_timestep!(   ocean::PrognosticVariablesOcean,
                            time::DateTime,
                            ocean_model::SeasonalOceanClimatology)

    # escape immediately if Δt of ocean model hasn't passed yet
    (time - ocean.time) < ocean_model.Δt && return nothing

    # otherwise update ocean prognostic variables:
    ocean.time = time
    interpolate_monthly!(   ocean.sea_surface_temperature,
                            oean_model.monthly_temperature,
                            time)
    return nothing
end

function interpolate_monthly!(
    grid::Grid,
    monthly::Vector{Grid},
    time::DateTime,
) where Grid

    this_month = Dates.month(time)
    next_month = (this_month % 12) + 1      # mod for dec 12 -> jan 1

    # linear interpolation weight between the two months
    # TODO check whether this shifts the climatology by 1/2 a month
    weight = convert(eltype(grid), Dates.days(time-Dates.firstdayofmonth(time))/Dates.daysinmonth(time))

    @. grid = (1-weight) * monthly[this_month] +
                 weight  * monthly[next_month]

    return nothing
end

## CONSTANT OCEAN CLIMATOLOGY

export ConstantOceanClimatology
Base.@kwdef struct ConstantOceanClimatology <: AbstractOcean
    "[OPTION] path to the folder containing the land-sea mask file, pkg path default"
    path::String = "SpeedyWeather.jl/input_data"

    "[OPTION] filename of sea surface temperatures"
    file::String = "sea_surface_temperature.nc"

    "[OPTION] Variable name in netcdf file"
    varname::String = "sst"

    "[OPTION] Grid the sea surface temperature file comes on"
    file_Grid::Type{<:AbstractGrid} = FullGaussianGrid

    "[OPTION] The missing value in the data respresenting land"
    missing_value::Float64 = NaN
end

# generator function, just pass on kwargs
function ConstantOceanClimatology(SG::SpectralGrid; kwargs...)
    ConstantOceanClimatology(; kwargs...)
end

# nothing to initialize for model.ocean
initialize!(::ConstantOceanClimatology, ::PrimitiveEquation) = nothing

# initialize 
function initialize!(   
    ocean::PrognosticVariablesOcean,
    time::DateTime,
    ocean_model::ConstantOceanClimatology,
    model::PrimitiveEquation,
)
    Grid = typeof(ocean.sea_surface_temperature)
    (; nlat_half) = ocean.sea_surface_temperature
    monthly = [zeros(Grid, nlat_half) for _ in 1:12]
    load_monthly_climatology!(monthly, ocean_model)
    interpolate_monthly!(ocean.sea_surface_temperature, monthly, time)
end

function ocean_timestep!(
    ocean::PrognosticVariablesOcean,
    time::DateTime,
    ocean_model::ConstantOceanClimatology
)
    return nothing
end

## CONSTANT OCEAN CLIMATOLOGY

export AquaPlanet
Base.@kwdef struct AquaPlanet{NF} <: AbstractOcean
    "Number of latitude rings"
    nlat::Int

    "[OPTION] Temperature on the Equator [K]"
    temp_equator::NF = 302

    "[OPTION] Temperature at the poles [K]"
    temp_poles::NF = 273

    "Latitudinal temperature profile [K]"
    temp_lat::Vector{NF} = zeros(NF, nlat)
end

# generator function
function AquaPlanet(SG::SpectralGrid; kwargs...)
    (; NF, nlat) = SG
    AquaPlanet{NF}(; nlat, kwargs...)
end

# nothing to initialize for model.ocean
function initialize!(ocean_model::AquaPlanet, model::PrimitiveEquation)
    (; coslat²) = model.geometry
    (; temp_lat, temp_equator, temp_poles) = ocean_model
    @. temp_lat = (temp_equator-temp_poles)*coslat² + temp_poles
end

# initialize 
function initialize!(   
    ocean::PrognosticVariablesOcean,
    time::DateTime,
    ocean_model::AquaPlanet,
    model::PrimitiveEquation,
)
    (; sea_surface_temperature) = ocean
    for (j, ring) in enumerate(eachring(sea_surface_temperature))
        for ij in ring
            sea_surface_temperature[ij] = ocean_model.temp_lat[j]
        end
    end
end

function ocean_timestep!(
    ocean::PrognosticVariablesOcean,
    time::DateTime,
    ocean_model::AquaPlanet,
)
    return nothing
end