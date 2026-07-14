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

"""$(TYPEDSIGNATURES)
Ensemble-aware variant: appends the `ensemble_index` as the outermost (last, Julia
column-major) index when ensemble output is on (`ensemble_index > 0`), otherwise falls
back to the regular indices. Works for both time-varying and static variables (static
variables ignore `i`). Used by the ensemble-capable `ZarrOutput` backend."""
get_indices(i, variable::AbstractOutputVariable, ensemble_index::Integer) =
    ensemble_index > 0 ? (get_indices(i, variable)..., ensemble_index) : get_indices(i, variable)

get_indices(i, x::Val{true}, y::Val{true}, z::Val{true}, t::Val{true}) = (:, :, :, i)   # 3D + time
get_indices(i, x::Val{true}, y::Val{true}, z::Val{true}, t::Val{false}) = (:, :, :)     # 3D
get_indices(i, x::Val{true}, y::Val{true}, z::Val{false}, t::Val{true}) = (:, :, i)     # 2D + time
get_indices(i, x::Val{true}, y::Val{true}, z::Val{false}, t::Val{false}) = (:, :)       # 2D

is2D(variable::AbstractOutputVariable) = ~variable.dims_xyzt[3]
is3D(variable::AbstractOutputVariable) = variable.dims_xyzt[3]
is_land(variable::AbstractOutputVariable) = hasproperty(variable, :is_land) ? variable.is_land : false
hastime(variable::AbstractOutputVariable) = variable.dims_xyzt[4]

"""$(TYPEDSIGNATURES)
Name of the vertical dimension of `variable` in the output file or store:
the atmospheric "layer" dimension, or "soil_layer" for land variables.
Extend for custom output variables written on their own vertical dimension,
together with [`get_nlayers`](@ref) and [`define_dimension!`](@ref)."""
vertical_dimension(variable::AbstractOutputVariable) = is_land(variable) ? "soil_layer" : "layer"

"""$(TYPEDSIGNATURES)
Number of vertical layers `variable` is written on, as allocated in the scratch
fields of `output`. Extend for custom output variables written on their own
vertical dimension, see [`vertical_dimension`](@ref)."""
get_nlayers(output::AbstractOutput, variable::AbstractOutputVariable) =
    is_land(variable) ? size(output.field3Dland, 2) : size(output.field3D, 2)

"""$(TYPEDSIGNATURES)
Hook called by `define_variable!` before `variable` is defined in the output
file or store `dest`; no-op by default as the default dimensions (see
[`vertical_dimension`](@ref)) are defined upfront by `initialize!`. Custom
output variables written on their own vertical dimension extend this to lazily
define that dimension, typically via [`get_dimension`](@ref) and
[`define_coordinate!`](@ref) so that one method covers all output backends."""
define_dimension!(dest, variable::AbstractOutputVariable) = nothing

"""$(TYPEDSIGNATURES) Like `path` but returns `nothing` instead of throwing an error if the variable is not defined in the simulation."""
path_or_nothing(variable::AbstractOutputVariable, simulation) = try path(variable, simulation) catch FieldError; nothing end
