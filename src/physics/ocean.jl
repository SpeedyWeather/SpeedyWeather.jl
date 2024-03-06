abstract type AbstractOcean{NF, Grid} end

function Base.show(io::IO, O::AbstractOcean)
    println(io, "$(typeof(O)) <: AbstractOcean")
    keys = propertynames(O)
    print_fields(io, O, keys)
end

export SeasonalOceanClimatology
Base.@kwdef struct SeasonalOceanClimatology{NF, Grid<:AbstractGrid{NF}} <: AbstractOcean{NF, Grid}

    "number of latitudes on one hemisphere, Equator included"
    nlat_half::Int

    # OPTIONS
    "Time step used to update sea surface temperatures"
    Δt::Dates.Day = Dates.Day(3)

    "path to the folder containing the land-sea mask file, pkg path default"
    path::String = "SpeedyWeather.jl/input_data"

    "filename of sea surface temperatures"
    file::String = "sea_surface_temperature.nc"

    "Variable name in netcdf file"
    varname::String = "sst"

    "Grid the sea surface temperature file comes on"
    file_Grid::Type{<:AbstractGrid} = FullGaussianGrid

    "The missing value in the data respresenting land"
    missing_value::NF = NF(NaN)

    # to be filled from file
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

function initialize!(   ocean::PrognosticVariablesOcean,
                        time::DateTime,
                        model::PrimitiveEquation)
    ocean_timestep!(ocean, time, model, initialize=true)
end

# function barrier
function ocean_timestep!(   ocean::PrognosticVariablesOcean,
                            time::DateTime,
                            model::PrimitiveEquation;
                            initialize::Bool = false)
    ocean_timestep!(ocean, time, model.ocean; initialize)
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
    weight = convert(NF, Dates.days(time-Dates.firstdayofmonth(time))/Dates.daysinmonth(time))

    (; monthly_temperature) = ocean_model
    @. ocean.sea_surface_temperature = (1-weight) * monthly_temperature[this_month] +
                                          weight  * monthly_temperature[next_month]

    return nothing
end
