abstract type AbstractVegetation <: AbstractParameterization end

export NoVegetation
struct NoVegetation <: AbstractVegetation end
NoVegetation(SG::SpectralGrid) = NoVegetation()
initialize!(vegetation::NoVegetation, model::PrimitiveEquation) = nothing

function initialize!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    vegetation::NoVegetation,
    model::PrimitiveWet)
    # initialize by running a "timestep"
    timestep!(diagn, progn, vegetation, model.land.soil_moisture, model)
end

function timestep!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    vegetation::NoVegetation,
    model::PrimitiveWet)
    # a "timestep" of no vegetation is just to calculate the soil moisture availability
    soil_moisture_availability!(diagn, progn, vegetation, model.land.temperature, model)
end

function soil_moisture_availability!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    vegetation::NoVegetation,
    land::AbstractLandTemperature,
    model::PrimitiveWet,
)

    (; soil_moisture) = progn.land
    (; W_cap, W_wilt) = vegetation
    D_top = land.z₁
    D_root = land.z₂

    r = 1/(D_top*W_cap + D_root*(W_cap - W_wilt))

    for ij in eachgridpoint(soil_moisture_availability, high_cover, low_cover)
        soil_moisture_availability[ij] = r*D_top*soil_moisture[ij, 1]
    end

    return nothing
end

export VegetationClimatology
@kwdef struct VegetationClimatology{NF, Grid} <: AbstractVegetation

    "number of latitudes on one hemisphere, Equator included"
    nlat_half::Int

    # OPTIONS
    "[OPTION] Combine high and low vegetation factor, a in high + a*low [1]"
    low_veg_factor::NF = 0.8

    "[OPTION] Soil wetness at field capacity [volume fraction]"
    W_cap::NF = 0.3

    "[OPTION] Soil wetness at wilting point [volume fraction]"
    W_wilt::NF = 0.17

    "[OPTION] path to the folder containing the soil moisture file, pkg path default"
    path::String = "SpeedyWeather.jl/input_data"

    "[OPTION] filename of soil moisture"
    file::String = "vegetation.nc"

    "[OPTION] variable name in netcdf file"
    varname_vegh::String = "vegh"
    varname_vegl::String = "vegl"

    "[OPTION] Grid the soil moisture file comes on"
    file_Grid::Type{<:AbstractGrid} = FullGaussianGrid

    "[OPTION] The missing value in the data respresenting ocean"
    missing_value::NF = NF(NaN)

    # to be filled from file
    "High vegetation cover [1], interpolated onto Grid"
    high_cover::Grid = zeros(Grid, nlat_half)

    "Low vegetation cover [1], interpolated onto Grid"
    low_cover::Grid = zeros(Grid, nlat_half)
end

# generator function
function VegetationClimatology(SG::SpectralGrid; kwargs...)
    (; NF, GridVariable2D, nlat_half) = SG
    return VegetationClimatology{NF, GridVariable2D}(; nlat_half, kwargs...)
end

function initialize!(vegetation::VegetationClimatology, model::PrimitiveEquation)

    # LOAD NETCDF FILE
    if vegetation.path == "SpeedyWeather.jl/input_data"
        path = joinpath(@__DIR__, "../../../input_data", vegetation.file)
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

function initialize!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    veg::VegetationClimatology,
    model::PrimitiveEquation,
)
    # initialize land temperature by "running" the step at the current time
    timestep!(progn, diagn, veg, model)
end

# function barrier
function timestep!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    vegetation::VegetationClimatology,
    model::PrimitiveWet)

    # a "timestep" of vegetation climatology is just to calculate the soil moisture availability
    soil_moisture_availability!(diagn, progn, vegetation, model.land.temperature, model)
end

function soil_moisture_availability!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    vegetation::AbstractVegetation,
    land::AbstractLandTemperature,
    model::PrimitiveWet,
)
    (; soil_moisture_availability) = diagn.physics
    (; soil_moisture) = progn.land
    (; high_cover, low_cover, low_veg_factor) = vegetation
    (; W_cap, W_wilt) = vegetation
    D_top = land.z₁
    D_root = land.z₂

    @boundscheck grids_match(high_cover, low_cover, soil_moisture_availability) || throws(BoundsError)
    @boundscheck grids_match(soil_moisture, soil_moisture_availability, horizontal_only=true) || throws(BoundsError)
    @boundscheck size(soil_moisture, 2) >= 2    # defined for two layers

    # precalculate
    r = 1/(D_top*W_cap + D_root*(W_cap - W_wilt))

    for ij in eachgridpoint(soil_moisture_availability, high_cover, low_cover)
        
        # Fortran SPEEDY source/land_model.f90 line 111 origin unclear
        veg = max(0, high_cover[ij] + low_veg_factor*low_cover[ij])

        # Fortran SPEEDY documentation eq. 51
        soil_moisture_availability[ij] = r*(D_top*soil_moisture[ij, 1] +
            veg*D_root*max(soil_moisture[ij, 2]- W_wilt, 0))
    end
end