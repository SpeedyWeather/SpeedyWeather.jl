abstract type AbstractLandTemperature <: AbstractParameterization end

export SeasonalLandTemperature
@kwdef struct SeasonalLandTemperature{NF, Grid} <: AbstractLand

    "number of latitudes on one hemisphere, Equator included"
    nlat_half::Int

    # OPTIONS
    "[OPTION] path to the folder containing the land temperature file, pkg path default"
    path::String = "SpeedyWeather.jl/input_data"

    "[OPTION] filename of land surface temperatures"
    file::String = "land_surface_temperature.nc"

    "[OPTION] variable name in netcdf file"
    varname::String = "lst"

    "[OPTION] Grid the land surface temperature file comes on"
    file_Grid::Type{<:AbstractGrid} = FullGaussianGrid

    "[OPTION] The missing value in the data respresenting ocean"
    missing_value::NF = NF(NaN)

    # to be filled from file
    "Monthly land surface temperatures [K], interpolated onto Grid"
    monthly_temperature::Grid = zeros(Grid, nlat_half, 12)
end

# generator function
function SeasonalLandTemperature(SG::SpectralGrid; kwargs...)
    (; NF, GridVariable3D, nlat_half) = SG
    return SeasonalLandTemperature{NF, GridVariable3D}(; nlat_half, kwargs...)
end

function initialize!(land::SeasonalLandTemperature, model::PrimitiveEquation)
    (; monthly_temperature) = land

    # LOAD NETCDF FILE
    if land.path == "SpeedyWeather.jl/input_data"
        path = joinpath(@__DIR__, "../../../input_data", land.file)
    else
        path = joinpath(land.path, land.file)
    end
    ncfile = NCDataset(path)

    # read out netCDF data
    nx, ny, nt = ncfile.dim["lon"], ncfile.dim["lat"], ncfile.dim["time"]
    npoints = nx*ny
    fill_value = ncfile[land.varname].attrib["_FillValue"]
    lst = land.file_Grid(ncfile[land.varname].var[:, :, :], input_as=Matrix)
    lst[lst .=== fill_value] .= land.missing_value      # === to include NaN
    
    @boundscheck grids_match(monthly_temperature, lst, vertical_only=true) || throw(DimensionMismatch(monthly_temperature, lst))

    # create interpolator from grid in file to grid used in model
    interp = RingGrids.interpolator(Float32, monthly_temperature, lst)
    interpolate!(monthly_temperature, lst, interp)
    return nothing
end

function initialize!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    land::SeasonalLandTemperature,
    model::PrimitiveEquation,
)
    # initialize land temperature by "running" the step at the current time
    timestep!(progn, diagn, land, model)
end

function timestep!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    land::SeasonalLandTemperature,
    model::PrimitiveEquation,
)
    (; time) = progn.clock

    this_month = Dates.month(time)
    next_month = (this_month % 12) + 1      # mod for dec 12 -> jan 1

    # linear interpolation weight between the two months
    # TODO check whether this shifts the climatology by 1/2 a month
    (; monthly_temperature) = land
    (; soil_temperature) = progn.land
    NF = eltype(soil_temperature)
    weight = convert(NF, Dates.days(time-Dates.firstdayofmonth(time))/Dates.daysinmonth(time))

    k0 = 1   # first layer only
    for ij in eachgridpoint(soil_temperature)
        soil_temperature[ij, k0] = (1-weight) * monthly_temperature[ij, this_month] +
                                                weight  * monthly_temperature[ij, next_month]
    end

    # set other layers to the same temperature
    for k in eachgrid(soil_temperature)
        if k != k0
            soil_temperature[:, k] .= soil_temperature[:, k0]
        end
    end

    return nothing
end

## CONSTANT LAND CLIMATOLOGY
export ConstantLandTemperature
@kwdef struct ConstantLandTemperature{NF} <: AbstractLand
    "[OPTION] Globally constant temperature"
    temperature::NF = 285
end

# generator function
ConstantLandTemperature(SG::SpectralGrid; kwargs...) = ConstantLandTemperature{SG.NF}(; kwargs...)

initialize!(land::ConstantLandTemperature, model::PrimitiveEquation) = nothing
function initialize!(   
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    land::ConstantLandTemperature,
    model::PrimitiveEquation,
)
    set!(progn.land.surface_temperature, land.temperature)
end

# temperature is constant so do nothing during land timestep
timestep!(progn::PrognosticVariables, diagn::DiagnosticVariables, land::ConstantLandTemperature, args...) = nothing