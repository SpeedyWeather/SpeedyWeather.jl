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

"""$(TYPEDSIGNATURES)
Indices to write the output slice of `variable` at output step `i` into a store that keeps
the horizontal dimension *unravelled* into a single flat dimension, as a `Field`'s `data`
already is. The `x` and `y` flags of `variable.dims_xyzt` collapse into that one flat
dimension, so a 3D+time variable is written as `(:, :, i)` where the equivalent
[`get_indices`](@ref) for a rectangular `lon`/`lat` store gives `(:, :, :, i)`.
Used by the flat [`HEALPixOutput`](@ref) backend."""
get_flat_indices(i, variable::AbstractOutputVariable) =
    get_flat_indices(i, Val(OutputVariables.is3D(variable)), Val(OutputVariables.hastime(variable)))

get_flat_indices(i, z::Val{true}, t::Val{true}) = (:, :, i)     # 3D + time
get_flat_indices(i, z::Val{true}, t::Val{false}) = (:, :)       # 3D
get_flat_indices(i, z::Val{false}, t::Val{true}) = (:, i)       # 2D + time
get_flat_indices(i, z::Val{false}, t::Val{false}) = (:,)        # 2D

"""$(TYPEDSIGNATURES)
Name of the vertical dimension of `variable` in the output file or store:
the atmospheric "layer" dimension, or "soil_layer" for land variables.
Extend for custom output variables written on their own vertical dimension,
together with [`get_nlayers`](@ref) and [`define_dimension!`](@ref)."""
vertical_dimension_name(variable::AbstractOutputVariable) = OutputVariables.is_land(variable) ? "soil_layer" : "layer"

"""$(TYPEDSIGNATURES)
Number of vertical layers `variable` is written on, as allocated in the scratch
fields of `output`. Extend for custom output variables written on their own
vertical dimension, see [`vertical_dimension_name`](@ref)."""
get_nlayers(output::AbstractOutput, variable::AbstractOutputVariable) =
    OutputVariables.is_land(variable) ? size(output.field3Dland, 2) : size(output.field3D, 2)

"""$(TYPEDSIGNATURES)
Function called by `define_variable!` before `variable` is defined in the output
file or store `dest`; no-op by default as the default dimensions (see
[`vertical_dimension_name`](@ref)) are defined upfront by `initialize!`. Custom
output variables written on their own vertical dimension extend this to
define that dimension, typically via [`get_dimension_length`](@ref) and
[`define_coordinate!`](@ref) so that one method covers all output backends."""
define_dimension!(dest, variable::AbstractOutputVariable) = nothing

"""$(TYPEDSIGNATURES)
Warn once about the output variables of `output` that don't exist in `simulation` and will
therefore be skipped, both in the file/store schema and on every output step."""
function warn_nonexisting_variables(output::AbstractOutput, simulation)
    skipped = [key for (key, var) in output.variables if !OutputVariables.exists_in_simulation(var, simulation)]
    isempty(skipped) || @warn "Some output.variables do not exist in simulation, skipping: " *
        "$(join(skipped, ", ")). This usually means the output writer was constructed with a " *
        "different Model type than the model it is attached to."
    return nothing
end


# the output variables themselves live in the `OutputVariables` submodule, which
# imports the helpers above, so it is included after them
include("OutputVariables.jl")
include("deprecated.jl")
