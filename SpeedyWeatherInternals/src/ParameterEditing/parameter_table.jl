"""
    ParameterTable{NT<:NamedTuple} <: ModelParameters.AbstractModel

Lightweight wrapper around a `NamedTuple` of parameters following the ModelParameters `AbstractModel` interface.
This provides a table-like interface for interacting with model parameters.
"""
struct ParameterTable{NT <: NamedTuple} <: ModelParameters.AbstractModel
    parent::NT
    ParameterTable(parent::NamedTuple) = new{typeof(parent)}(parent)
    ParameterTable(; params...) = ParameterTable((; params...))
end

"""
    $SIGNATURES

Simple recursive dispatch algorithm for unpacking nested named tuples of `ParameterTable`
"""
unpack_params(x) = x
unpack_params(nt::NamedTuple) = map(unpack_params, nt)
unpack_params(params::ParameterTable) = unpack_params(parent(params))

"""
    $SIGNATURES

Unpacks and strips all `AbstractParam` types from the given `params` object, replacing them with their numeric `val`.
This is a type-stable alternative to the (currently unstable) `stripparams` in `ModelParameters`.
"""
stripparams(param::AbstractParam) = value(param)
stripparams(params::NamedTuple) = map(stripparams, params)
stripparams(params::ParameterTable) = stripparams(unpack_params(params))

# Base overrides for ParameterTable

# parameter subsets
Base.getindex(ps::ParameterTable, param_label::String) = getindex(ps, [param_label])
@inline function Base.getindex(ps::ParameterTable, param_labels::Vector{String})
    # extract labels from ComponentVector
    ls = labels(vec(ps))
    idx = reduce(vcat, map(query -> findall(key -> startswith(key, query), ls), param_labels))
    # check if any names were not found
    length(idx) > 0 || @warn "no parameters matched to labels: $(param_labels)"
    # get flattened parameter tuple
    params = ModelParameters.params(ps)
    # get tuple of selected parameters
    selected_params = params[idx]
    param_subset = _selectrecursive(parent(ps)) do query
        # select only the parameters identical to (===) those in selected_params
        any(map(key -> key === query, selected_params))
    end
    # extract parameters from nested named tuple and reconstruct ParameterTable
    return ParameterTable(param_subset)
end

# selecting columns
Base.getindex(ps::ParameterTable, nm::Symbol) = getindex(ps, :, nm)
@inline function Base.getindex(ps::ParameterTable, ::Colon, nm::Symbol)
    # include fieldname only (no component)
    return if nm == :idx
        1:length(ps)
    elseif nm == :fieldname
        ModelParameters.paramfieldnames(ps)
    elseif any(map(p -> hasproperty(p, nm), ModelParameters.params(ps)))
        map(ModelParameters.params(ps)) do p
            hasproperty(p, nm) ? getproperty(p, nm) : nothing
        end
    else
        error("$nm is not a property of any parameter")
    end
end

# non-mutating setindex!
# consider using a method instead of setindex! to avoid confusion?
Base.setindex!(ps::ParameterTable, x, nm::Symbol) = setindex!(ps, x, :, nm)
@inline function Base.setindex!(ps::ParameterTable, x, ::Colon, nm::Symbol)
    # include fieldname only (no component)
    if nm == :idx
        error("cannot set :idx")
    elseif nm == :fieldname
        error("cannot set :fieldname")
    else
        newparent = if nm in keys(ps)
            ModelParameters._setindex(ModelParameters.parent(ps), Tuple(x), nm)
        else
            ModelParameters._addindex(ModelParameters.parent(ps), Tuple(x), nm)
        end
        return ParameterTable(newparent)
    end
end

Base.keys(params::ParameterTable) = (:idx, :fieldname, union(map(keys, params)...)...)

function Base.vec(params::ParameterTable)
    # recursively unpack params and build component vector
    return ComponentArray(unpack_params(params))
end

# Override internal ModelParameters method _columntypes to condense type names in table schema (just to look nicer)
# TODO: propose this change upstream in ModelParameters.jl
ModelParameters._columntypes(ps::ParameterTable) = map(k -> promote_type(map(typeof, getindex(ps, k))...), keys(ps))
