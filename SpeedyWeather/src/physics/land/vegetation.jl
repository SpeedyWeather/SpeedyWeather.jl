abstract type AbstractVegetation <: AbstractLandComponent end

export NoVegetation
struct NoVegetation <: AbstractVegetation end
NoVegetation(SG::SpectralGrid, geometry::LandGeometryOrNothing = nothing) = NoVegetation()
initialize!(vegetation::NoVegetation, model::PrimitiveEquation) = nothing

function initialize!(
        progn::PrognosticVariables,
        diagn::DiagnosticVariables,
        vegetation::NoVegetation,
        model::PrimitiveEquation
    )
    # initialize by running a "timestep"
    return timestep!(progn, diagn, vegetation, model)
end

function timestep!(
        progn::PrognosticVariables,
        diagn::DiagnosticVariables,
        vegetation::NoVegetation,
        model::PrimitiveEquation
    )
    # a "timestep" of no vegetation is just to calculate the soil moisture availability
    return soil_moisture_availability!(diagn, progn, vegetation, model)
end

function soil_moisture_availability!(
        diagn::DiagnosticVariables,
        progn::PrognosticVariables,
        vegetation::AbstractVegetation,
        model::PrimitiveDry,
    )
    return nothing
end

function soil_moisture_availability!(
        diagn::DiagnosticVariables,
        progn::PrognosticVariables,
        vegetation::NoVegetation,
        model::PrimitiveWet,
    )
    # view on the top layer of soil moisture
    soil_moisture_top = field_view(progn.land.soil_moisture, :, 1)
    (; soil_moisture_availability) = diagn.physics.land

    # Fortran SPEEDY documentation eq. 51 with vegetation = 0
    W_cap = model.land.thermodynamics.field_capacity
    W_wilt = model.land.thermodynamics.wilting_point
    D_top = model.land.geometry.layer_thickness[1]
    D_root = model.land.geometry.layer_thickness[2]

    soil_moisture_availability .= D_top * soil_moisture_top * W_cap /
        (D_top * W_cap + D_root * (W_cap - W_wilt))
    return nothing
end

function variables(::NoVegetation)
    return (
        PrognosticVariable(name = :soil_moisture, dims = Grid3D(), desc = "Soil moisture content (fraction of capacity)", units = "1", namespace = :land),
        DiagnosticVariable(name = :soil_moisture_availability, dims = Grid2D(), desc = "Soil moisture availability for evaporation", units = "1", namespace = :land),
    )
end

export VegetationClimatology
@kwdef struct VegetationClimatology{NF, GridVariable2D} <: AbstractVegetation
    "[OPTION] Combine high and low vegetation factor, a in high + a*low [1]"
    low_veg_factor::NF = 0.8

    "[OPTION] filename of soil moisture"
    file::String = "vegetation.nc"

    "[OPTION] variable name in netcdf file for high vegetation"
    varname_vegh::String = "vegh"

    "[OPTION] variable name in netcdf file for low vegetation"
    varname_vegl::String = "vegl"

    "[OPTION] Grid the soil moisture file comes on"
    file_Grid::Type{<:AbstractGrid} = FullGaussianGrid

    "[OPTION] The missing value in the data respresenting ocean"
    missing_value::NF = NaN
    # to be filled from file
    "High vegetation cover [1], interpolated onto Grid"
    high_cover::GridVariable2D

    "Low vegetation cover [1], interpolated onto Grid"
    low_cover::GridVariable2D
end

# TODO to adapt create a ManualVegetationClimatology component like AlbedoClimatology is adapted to ManualAlbedo
# do all vegetations need a low_veg_factor?
# Adapt.adapt_structure(to, veg::VegetationClimatology) = adapt(to, ManualVegetationClimatology(veg.high_cover, veg.low_cover))

# generator function
function VegetationClimatology(SG::SpectralGrid, geometry::LandGeometryOrNothing = nothing; kwargs...)
    (; NF, GridVariable2D, grid) = SG
    high_cover = zeros(GridVariable2D, grid)
    low_cover = zeros(GridVariable2D, grid)
    return VegetationClimatology{NF, GridVariable2D}(; high_cover, low_cover, kwargs...)
end

function initialize!(vegetation::VegetationClimatology, model::PrimitiveEquation)

    # LOAD NETCDF FILE
    path = get_asset("data", vegetation.file)
    ncfile = NCDataset(path)

    # high and low vegetation cover
    vegh = vegetation.file_Grid(ncfile[vegetation.varname_vegh].var[:, :], input_as = Matrix)
    vegl = vegetation.file_Grid(ncfile[vegetation.varname_vegl].var[:, :], input_as = Matrix)
    vegh = on_architecture(model.architecture, vegh)
    vegl = on_architecture(model.architecture, vegl)

    # interpolate onto grid
    high_vegetation_cover = vegetation.high_cover
    low_vegetation_cover = vegetation.low_cover
    interpolator = RingGrids.interpolator(high_vegetation_cover, vegh, NF = Float32)
    interpolate!(high_vegetation_cover, vegh, interpolator)
    return interpolate!(low_vegetation_cover, vegl, interpolator)
end

function initialize!(
        progn::PrognosticVariables,
        diagn::DiagnosticVariables,
        veg::VegetationClimatology,
        model::PrimitiveEquation,
    )
    # initialize land temperature by "running" the step at the current time
    return timestep!(progn, diagn, veg, model)
end

# function barrier
function timestep!(
        progn::PrognosticVariables,
        diagn::DiagnosticVariables,
        vegetation::VegetationClimatology,
        model::PrimitiveEquation
    )

    # a "timestep" of vegetation climatology is just to calculate the soil moisture availability
    return soil_moisture_availability!(diagn, progn, vegetation, model)
end

function soil_moisture_availability!(
        diagn::DiagnosticVariables,
        progn::PrognosticVariables,
        vegetation::VegetationClimatology,
        model::PrimitiveDry,
    )
    return nothing
end

function soil_moisture_availability!(
        diagn::DiagnosticVariables,
        progn::PrognosticVariables,
        vegetation::VegetationClimatology,
        model::PrimitiveWet,
    )
    (; vegetation_high, vegetation_low, soil_moisture_availability) = diagn.physics.land
    (; soil_moisture) = progn.land
    (; low_veg_factor) = vegetation

    # copy over vegetation fields into diagnostic variables
    vegetation_high .= vegetation.high_cover
    vegetation_low .= vegetation.low_cover

    W_cap = model.land.thermodynamics.field_capacity
    W_wilt = model.land.thermodynamics.wilting_point
    D_top = model.land.geometry.layer_thickness[1]
    D_root = model.land.geometry.layer_thickness[2]

    @boundscheck fields_match(vegetation_high, vegetation_low, soil_moisture_availability) || throws(BoundsError)
    @boundscheck fields_match(soil_moisture, soil_moisture_availability, horizontal_only = true) || throws(BoundsError)
    @boundscheck size(soil_moisture, 2) >= 2                # defined for two layers
    @boundscheck size(soil_moisture_availability, 2) == 1   # 2D only

    # precalculate denominator
    r = 1 / (D_top * W_cap + D_root * (W_cap - W_wilt))

    # pack parameters into named tuple for kernel
    params = (; low_veg_factor, r, W_cap, W_wilt, D_top, D_root)

    launch!(
        architecture(soil_moisture_availability), LinearWorkOrder,
        size(soil_moisture_availability), soil_moisture_availability_kernel!,
        soil_moisture_availability, soil_moisture, vegetation_high, vegetation_low, params
    )

    return nothing
end

@kernel inbounds = true function soil_moisture_availability_kernel!(
        soil_moisture_availability, soil_moisture, vegetation_high, vegetation_low, params
    )

    ij = @index(Global, Linear)    # every grid point ij

    (; low_veg_factor, r, W_cap, W_wilt, D_top, D_root) = params

    # Fortran SPEEDY source/land_model.f90 line 111 origin unclear
    veg = max(0, vegetation_high[ij] + low_veg_factor * vegetation_low[ij])

    # Fortran SPEEDY documentation eq. 51, original formulation
    # soil_moisture_availability[ij] = r*(D_top*soil_moisture[ij, 1] +
    #     veg*D_root*max(soil_moisture[ij, 2] - W_wilt, 0))
    # Soil moisture is defined as volume fraction wrt to field capacity
    # so multiply by W_cap here (not done in Fortran SPEEDY)
    soil_moisture_availability[ij] = r * (
        D_top * soil_moisture[ij, 1] * W_cap +
            veg * D_root * max(soil_moisture[ij, 2] * W_cap - W_wilt, 0)
    )
end

function variables(::VegetationClimatology)
    return (
        DiagnosticVariable(name = :vegetation_high, dims = Grid2D(), desc = "Vegetation high cover", units = "1", namespace = :land),
        DiagnosticVariable(name = :vegetation_low, dims = Grid2D(), desc = "Vegetation low cover", units = "1", namespace = :land),
        DiagnosticVariable(name = :soil_moisture_availability, dims = Grid2D(), desc = "Soil moisture availability for evaporation", units = "1", namespace = :land),
        PrognosticVariable(name = :soil_moisture, dims = Grid3D(), desc = "Soil moisture content (fraction of capacity)", units = "1", namespace = :land),
    )
end
