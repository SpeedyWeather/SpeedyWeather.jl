"""
    NumberParam{NF <: Number} <: AbstractParam

Specialized implementation of `AbstractParam` for Speedy models that defines a fixed
set of fields with concrete types.
"""
@kwdef struct NumberParam{NF <: Number} <: AbstractParam{NF}
    "numeric value of the parameter"
    value::NF

    "numerical domain on which the parameter is defined"
    bounds::Domain = Unbounded

    "human-readable description of the parameter"
    desc::String = ""

    "additional attributes with unspecified types"
    attrs::NamedTuple = (;)

    # default constructor
    NumberParam(value::NF, bounds::Domain, desc::String, attrs::NamedTuple) where {NF <: AbstractFloat} = new{NF}(value, bounds, desc, attrs)
    # convenience constructor
    NumberParam(value::NF; bounds = Unbounded, desc = "", attrs...) where {NF <: AbstractFloat} = new{NF}(value, bounds, desc, NamedTuple(attrs))
    # mandatory constructor from ModelParameters that allows for automated reconstruction
    NumberParam(nt::NamedTuple) = new{typeof(nt.val)}(nt.val, nt.bounds, nt.desc, _attrs(nt))
end

# Mandatory ModelParameters `AbstractParam` interface methods
ModelParameters.parent(param::NumberParam) = (val = value(param), attributes(param)..., bounds = bounds(param), desc = description(param))
ModelParameters.rebuild(param::NumberParam, newvalues) = NumberParam(newvalues)

"""
Returns the value of the given `param`.
"""
value(param::NumberParam) = getfield(param, :value)

"""
Returns the numerical bounds of the parameter `param`.
"""
bounds(param::NumberParam) = getfield(param, :bounds)

"""
Returns a human-readable description of the parameter `param`.
"""
description(param::NumberParam) = getfield(param, :desc)

"""
Returns a `NamedTuple` of any additional attributes defined for `param`.
"""
attributes(param::NumberParam) = getfield(param, :attrs)

# parameters method interface
"""
    $SIGNATURES

Extract parameters from the given `obj` as (possibly nested) named-tuple of `NumberParam`s or some other
`AbstractParam` type. If `obj`
"""
parameters(obj; kwargs...) = (;)
parameters(param::PT; kwargs...) where {PT <: AbstractParam} = parameters(PT, param; kwargs...)
parameters(param::Union{Number, AbstractArray}; kwargs...) = parameters(NumberParam, param; kwargs...)
parameters(::Type{PT}, obj; kwargs...) where {PT <: AbstractParam} = parameters(obj; kwargs...)
parameters(::Type{PT}, param::AbstractParam; kwargs...) where {PT <: AbstractParam} = PT(merge(parent(param), kwargs))
parameters(::Type{PT}, x::Union{Number, AbstractArray}; kwargs...) where {PT <: AbstractParam} = PT(x; kwargs...)

"""
    $SIGNATURES

Convenience method that creates a model parameter from its property with the given `name` and optional extra attributes in `kwargs`.
A parameter attribute `copmonenttype` is automatically added with value `T`.
"""
parameterof(obj::T, ::Val{propname}; kwargs...) where {T, propname} = parameterof(NumberParam, obj, Val{propname}(); kwargs...)
parameterof(::Type{PT}, obj::T, ::Val{propname}; kwargs...) where {PT <: AbstractParam, T, propname} = parameters(PT, getproperty(obj, propname); merge((; kwargs...), (component_type = T,))...)

# reconstruct

"""
    $SIGNATURES

Reconstructs the given (possibly nested) data structure with the given `values`. If `values` is a `NamedTuple`,
the nested structure must match that of `obj`. This function is used to reconstruct model types from `ParameterTable`.
"""
@inline reconstruct(obj::T, value::T) where {T} = value
@inline reconstruct(obj::AbstractParam, value::T) where {T} = ModelParameters.update(obj, Tuple(value))
@inline reconstruct(obj::NamedTuple{keys, V}, values::NamedTuple{keys, V}) where {keys, V <: Tuple} = values
@inline reconstruct(obj, values::ParameterTable) = reconstruct(obj, stripparams(values))
@inline reconstruct(obj, values::AbstractArray) = isempty(values) ? obj : error("Cannot reconstruct $(typeof(obj)) from non-empty array of type $(typeof(values))")
@generated function reconstruct(obj, values::Union{NamedTuple, ComponentArray})
    keysof(::Type{<:NamedTuple{keys}}) where {keys} = keys
    keysof(::Type{<:ComponentArray{T, N, A, Tuple{Axis{coords}}}}) where {T, N, A, coords} = keys(coords)
    recursive_calls = map(k -> :(reconstruct(obj.$k, values.$k)), keysof(values))
    return quote
        # recursively call reconstruct for all keys specified in values
        patchvals = tuple($(recursive_calls...))
        # construct a named tuple with the reconstructed values and apply with setproperties
        patch = NamedTuple{$(keysof(values))}(patchvals)
        return setproperties(obj, patch)
    end
end

# Internal helper functions

"""
Internal helper function for filtering out built-in parameter fields.
"""
_attrs(nt::NamedTuple) = nt[filter(k -> k ∉ (:val, :bounds, :desc), keys(nt))]

"""
    _selectrecursive(selector, nt::NamedTuple)

Recursively filters out all values from a (possibly nested) named tuple `nt` for which `selector`
returns true.
"""
_selectrecursive(selector, ps::ParameterTable) = _selectrecursive(selector, parent(ps))
_selectrecursive(selector, x) = selector(x) ? x : nothing
function _selectrecursive(selector, nt::NamedTuple)
    # recursively apply selector
    new_nt = map(x -> _selectrecursive(selector, x), nt)
    # filter out all nothing values
    selected_keys = filter(k -> !isnothing(new_nt[k]), keys(new_nt))
    # construct named tuple from filtered keys or return nothing if no keys were selected
    return length(selected_keys) > 0 ? NamedTuple{selected_keys}(map(k -> new_nt[k], selected_keys)) : nothing
end
