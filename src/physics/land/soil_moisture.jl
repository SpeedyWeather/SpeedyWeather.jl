abstract type AbstractSoil <: AbstractParameterization end

export SeasonalSoilMoisture
@kwdef struct SeasonalSoilMoisture{NF, Grid} <: AbstractSoil

    "number of latitudes on one hemisphere, Equator included"
    nlat_half::Int

    "number of soil layers"
    nlayers::Int

    # OPTIONS
    "[OPTION] Depth of top soil layer [m]"
    D_top::NF = 0.07

    "[OPTION] Depth of root layer [m]"
    D_root::NF = 0.21

    "[OPTION] Soil wetness at field capacity [volume fraction]"
    W_cap::NF = 0.3

    "[OPTION] Soil wetness at wilting point [volume fraction]"
    W_wilt::NF = 0.17

    # READ CLIMATOLOGY FROM FILE
    "[OPTION] path to the folder containing the soil moisture file, pkg path default"
    path::String = "SpeedyWeather.jl/input_data"

    "[OPTION] filename of soil moisture"
    file::String = "soil_moisture.nc"

    "[OPTION] variable name in netcdf file"
    varname_layer1::String = "swl1"
    varname_layer2::String = "swl2"

    "[OPTION] Grid the soil moisture file comes on"
    file_Grid::Type{<:AbstractGrid} = FullGaussianGrid

    "[OPTION] The missing value in the data respresenting ocean"
    missing_value::NF = NF(NaN)

    # to be filled from file
    "Monthly soil moisture volume fraction [1], interpolated onto Grid"
    monthly_soil_moisture::Grid = zeros(Grid, nlat_half, nlayers, 12)
end

# generator function
function SeasonalSoilMoisture(SG::SpectralGrid; kwargs...)
    (; NF, GridVariable4D, nlayers_soil, nlat_half) = SG
    return SeasonalSoilMoisture{NF, GridVariable4D}(; nlat_half, nlayers=nlayers_soil, kwargs...)
end

function initialize!(soil::SeasonalSoilMoisture, model::PrimitiveWet)
    (; monthly_soil_moisture) = soil

    # LOAD NETCDF FILE
    if soil.path == "SpeedyWeather.jl/input_data"
        path = joinpath(@__DIR__, "../../../input_data", soil.file)
    else
        path = joinpath(soil.path, soil.file)
    end
    ncfile = NCDataset(path)

    # read out netCDF data
    nx, ny, nt = ncfile.dim["lon"], ncfile.dim["lat"], ncfile.dim["time"]
    nlat_half = RingGrids.get_nlat_half(soil.file_Grid, nx*ny)

    # the soil moisture from file but wrapped into a grid
    NF = eltype(monthly_soil_moisture)
    soil_moisture_file = zeros(soil.file_Grid{NF}, nlat_half, soil.nlayers, nt)

    for l in 1:nt   # read out monthly to swap dimensions to horizontal - vertical - time
        soil_moisture_file[:, 1, l] .= vec(ncfile[soil.varname_layer1].var[:, :, l])
        soil_moisture_file[:, 2, l] .= vec(ncfile[soil.varname_layer2].var[:, :, l])
    end

    fill_value1 = NF(ncfile[soil.varname_layer1].attrib["_FillValue"])
    fill_value2 = NF(ncfile[soil.varname_layer2].attrib["_FillValue"])
    fill_value1 === fill_value2 || @warn "Fill values are different for the two soil layers, use only from layer 1"
    soil_moisture_file[soil_moisture_file .=== fill_value1] .= soil.missing_value      # === to include NaN
    
    @boundscheck grids_match(monthly_soil_moisture, soil_moisture_file, vertical_only=true) ||
        throw(DimensionMismatch(monthly_soil_moisture, soil_moisture_file))

    # create interpolator from grid in file to grid used in model
    interp = RingGrids.interpolator(Float32, monthly_soil_moisture, soil_moisture_file)
    interpolate!(monthly_soil_moisture, soil_moisture_file, interp)
    return nothing
end

function initialize!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    soil::SeasonalSoilMoisture,
    model::PrimitiveEquation,
)
    # initialize land temperature by "running" the step at the current time
    timestep!(progn, diagn, soil, model)
end

function timestep!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    soil::SeasonalSoilMoisture,
    model::PrimitiveEquation,
)
    (; time) = progn.clock

    this_month = Dates.month(time)
    next_month = (this_month % 12) + 1      # mod for dec 12 -> jan 1

    # linear interpolation weight between the two months
    # TODO check whether this shifts the climatology by 1/2 a month
    NF = eltype(progn.land.soil_moisture)
    weight = convert(NF, Dates.days(time-Dates.firstdayofmonth(time))/Dates.daysinmonth(time))

    (; monthly_soil_moisture) = soil
    (; soil_moisture) = progn.land

    for k in eachgrid(soil_moisture)
        for ij in eachgridpoint(soil_moisture)
            soil_moisture[ij, k] = (1-weight) * monthly_soil_moisture[ij, k, this_month] +
                                    weight  * monthly_soil_moisture[ij, k, next_month]
        end
    end

    return nothing
end