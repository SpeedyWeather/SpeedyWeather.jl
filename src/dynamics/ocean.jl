abstract type AbstractOcean{NF,Grid} end

function Base.show(io::IO,O::AbstractOcean)
    println(io,"$(typeof(O)) <: AbstractOcean")
    keys = propertynames(O)
    print_fields(io,O,keys)
end

Base.@kwdef struct SeasonalOceanClimatology{NF,Grid<:AbstractGrid{NF}} <: AbstractOcean{NF,Grid}

    "number of latitudes on one hemisphere, Equator included"
    nlat_half::Int

    # OPTIONS
    "Time step used to update sea surface temperatures"
    Δt::Dates.Day = Dates.Day(3)

    "path to the folder containing the land-sea mask file, pkg path default"
    path::String = "SpeedyWeather.jl/input_data"

    "filename of sea surface temperatures"
    file::String = "sea_surface_temperature.nc"

    "Grid the sea surface temperature file comes on"
    file_Grid::Type{<:AbstractGrid} = FullGaussianGrid

    "The missing value in the data respresenting land"
    missing_value::NF = NF(NaN)

    # to be filled from file
    "Monthly sea surface temperatures [K], interpolated onto Grid"
    monthly_temperature::Vector{Grid} = [zeros(Grid,nlat_half) for _ in 1:12]
end

# generator function
function SeasonalOceanClimatology(SG::SpectralGrid;kwargs...)
    (;NF,Grid,nlat_half) = SG
    return SeasonalOceanClimatology{NF,Grid{NF}}(;nlat_half,kwargs...)
end

function initialize!(ocean::SeasonalOceanClimatology{NF,Grid}) where {NF,Grid}
    # LOAD NETCDF FILE
    if ocean.path == "SpeedyWeather.jl/input_data"
        path = joinpath(@__DIR__,"../../input_data",ocean.file)
    else
        path = joinpath(ocean.path,ocean.file)
    end
    ncfile = NCDataset(path)

    # create interpolator from grid in file to grid used in model
    nx, ny = ncfile.dim["lon"], ncfile.dim["lat"]
    npoints = nx*ny
    NF_file = typeof(ncfile["sst"].attrib["_FillValue"])
    sst = ocean.file_Grid(zeros(NF_file,npoints))
    interp = RingGrids.interpolator(NF,ocean.monthly_temperature[1],sst)

    # interpolate and store in ocean
    for month in 1:12
        sst_this_month = ncfile["sst"][:,:,month]
        ij = 0
        for j in 1:ny
            for i in 1:nx
                ij += 1
                x = sst_this_month[i,j]
                sst[ij] = ismissing(x) ? ocean.missing_value : x
            end
        end
        interpolate!(ocean.monthly_temperature[month],sst,interp)
    end
    return nothing
end

function initialize!(   ocean::PrognosticVariablesOcean,
                        time::DateTime,
                        model::PrimitiveEquation)
    ocean_timestep!(ocean,time,model,initialize=true)
end

# function barrier
function ocean_timestep!(   ocean::PrognosticVariablesOcean,
                            time::DateTime,
                            model::PrimitiveEquation;
                            initialize::Bool = false)
    ocean_timestep!(ocean,time,model.ocean;initialize)
end

function ocean_timestep!(   ocean::PrognosticVariablesOcean{NF},
                            time::DateTime,
                            ocean_model::SeasonalOceanClimatology;
                            initialize::Bool = false) where NF

    # escape immediately if Δt of ocean model hasn't passed yet
    # unless the ocean hasn't been initialized yet
    initialize || (time - ocean.time) < ocean_model.Δt && return nothing

    # otherwise update ocean prognostic variables:
    ocean.time = time
    this_month = Dates.month(time)
    next_month = (this_month % 12) + 1      # mod for dec 12 -> jan 1

    # linear interpolation weight between the two months
    # TODO check whether this shifts the climatology by 1/2 a month
    weight = convert(NF,Dates.days(time-Dates.firstdayofmonth(time))/Dates.daysinmonth(time))

    (;monthly_temperature) = ocean_model
    @. ocean.sea_surface_temperature = (1-weight) * monthly_temperature[this_month] +
                                          weight  * monthly_temperature[next_month]

    return nothing
end
