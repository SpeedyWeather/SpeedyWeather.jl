abstract type AbstractLand{NF,Grid} end

function Base.show(io::IO,L::AbstractLand)
    println(io,"$(typeof(L)) <: AbstractLand")
    keys = propertynames(L)
    print_fields(io,L,keys)
end

Base.@kwdef struct SeasonalLandClimatology{NF,Grid<:AbstractGrid{NF}} <: AbstractLand{NF,Grid}

    "number of latitudes on one hemisphere, Equator included"
    nlat_half::Int

    # OPTIONS
    "Time step used to update land surface temperatures"
    Δt::Dates.Day = Dates.Day(3)

    "path to the folder containing the land-sea mask file, pkg path default"
    path::String = "SpeedyWeather.jl/input_data"

    "filename of sea surface temperatures"
    file::String = "land_surface_temperature.nc"

    "Grid the sea surface temperature file comes on"
    file_Grid::Type{<:AbstractGrid} = FullGaussianGrid

    "The missing value in the data respresenting land"
    missing_value::NF = NF(NaN)

    # to be filled from file
    "Monthly land surface temperatures [K], interpolated onto Grid"
    monthly_temperature::Vector{Grid} = [zeros(Grid,nlat_half) for _ in 1:12]
end

# generator function
function SeasonalLandClimatology(SG::SpectralGrid;kwargs...)
    (;NF,Grid,nlat_half) = SG
    return SeasonalLandClimatology{NF,Grid{NF}}(;nlat_half,kwargs...)
end

function initialize!(land::SeasonalLandClimatology{NF,Grid}) where {NF,Grid}
    # LOAD NETCDF FILE
    if land.path == "SpeedyWeather.jl/input_data"
        path = joinpath(@__DIR__,"../../input_data",land.file)
    else
        path = joinpath(land.path,land.file)
    end
    ncfile = NCDataset(path)

    # create interpolator from grid in file to grid used in model
    nx, ny = ncfile.dim["lon"], ncfile.dim["lat"]
    npoints = nx*ny
    NF_file = typeof(ncfile["lst"].attrib["_FillValue"])
    lst = land.file_Grid(zeros(NF_file,npoints))
    interp = RingGrids.interpolator(NF,land.monthly_temperature[1],lst)

    # interpolate and store in land
    for month in 1:12
        lst_this_month = ncfile["lst"][:,:,month]
        ij = 0
        for j in 1:ny
            for i in 1:nx
                ij += 1
                x = lst_this_month[i,j]
                lst[ij] = ismissing(x) ? land.missing_value : x
            end
        end
        interpolate!(land.monthly_temperature[month],lst,interp)
    end
    return nothing
end

function initialize!(   land::PrognosticVariablesLand,
                        time::DateTime,
                        model::PrimitiveEquation)
    land_timestep!(land,time,model,initialize=true)
end

# function barrier
function land_timestep!(land::PrognosticVariablesLand,
                        time::DateTime,
                        model::PrimitiveEquation;
                        initialize::Bool = false)
    land_timestep!(land,time,model.land;initialize)
end

function land_timestep!(land::PrognosticVariablesLand{NF},
                        time::DateTime,
                        land_model::SeasonalLandClimatology;
                        initialize::Bool = false) where NF

    # escape immediately if Δt of land model hasn't passed yet
    # unless the land hasn't been initialized yet
    initialize || (time - land.time) < land_model.Δt && return nothing

    # otherwise update land prognostic variables:
    land.time = time
    this_month = Dates.month(time)
    next_month = (this_month % 12) + 1      # mod for dec 12 -> jan 1

    # linear interpolation weight between the two months
    # TODO check whether this shifts the climatology by 1/2 a month
    weight = convert(NF,Dates.days(time-Dates.firstdayofmonth(time))/Dates.daysinmonth(time))

    (;monthly_temperature) = land_model
    @. land.land_surface_temperature = (1-weight) * monthly_temperature[this_month] +
                                          weight  * monthly_temperature[next_month]

    return nothing
end
