# define all variables for output
include("dynamics.jl")        # collected as DynamicsOutput()
include("precipitation.jl")   # PrecipitationOutput()
include("boundaries.jl")      # BoundaryOutput()
include("radiation.jl")       # RadiationOutput()
include("stochastic.jl")      # RandomPatternOutput()
include("surface_fluxes.jl")  # SurfaceFluxesOutput()
include("land.jl")            # LandOutput()
include("ocean.jl")           # SeaSurfaceTemperatureOutput()
include("tracers.jl")         # TracerOutput()
include("boundary_layer.jl")  # BoundaryLayerOutput()

# collect all together for conveneince
AllOutputVariables() = (
    DynamicsOutput()...,
    PrecipitationOutput()...,
    BoundaryOutput()...,
    RadiationOutput()...,
    RandomPatternOutput(),
    SurfaceFluxesOutput()...,
    LandOutput()...,
    OceanOutput()...,
    BoundaryLayerOutput()...,
)

get_indices(i, variable::AbstractOutputVariable) = get_indices(i, Val.(variable.dims_xyzt)...)
get_indices(i, x::Val{true}, y::Val{true}, z::Val{true}, t::Val{true}) = (:, :, :, i)   # 3D + time
get_indices(i, x::Val{true}, y::Val{true}, z::Val{true}, t::Val{false}) = (:, :, :)     # 3D
get_indices(i, x::Val{true}, y::Val{true}, z::Val{false}, t::Val{true}) = (:, :, i)     # 2D + time
get_indices(i, x::Val{true}, y::Val{true}, z::Val{false}, t::Val{false}) = (:, :)       # 2D

is3D(variable::AbstractOutputVariable) = variable.dims_xyzt[3]
is_land(variable::AbstractOutputVariable) = hasproperty(variable, :is_land) ? variable.is_land : false
hastime(variable::AbstractOutputVariable) = variable.dims_xyzt[4]
