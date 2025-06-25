using DomainSets: Domain, RealLine

"""
    SpeedyParam{NF<:AbstractFloat} <: ModelParameters.AbstractParam

Specialized implementation of `AbstractParam` for Speedy models that defines a fixed
set of fields with concrete types.
"""
Base.@kwdef struct SpeedyParam{NF<:AbstractFloat} <: ModelParameters.AbstractParam{NF}
    "numeric value of the parameter"
    value::NF = NaN
    
    "numerical domain on which the parameter is defined"
    bounds::Domain = RealLine()

    "human-readable description of the parameter"
    desc::String = ""

    "additional attributes with unspecified types"
    attrs::NamedTuple

    # default constructor
    SpeedyParam(value::NF, bounds::Domain, desc::String, attrs::NamedTuple) where {NF<:AbstractFloat} = new{NF}(value, bounds, desc, attrs)
    # convenience constructor
    SpeedyParam(value::NF; bounds=RealLine(), desc="", attrs...) where {NF<:AbstractFloat} = new{NF}(value, bounds, desc, NamedTuple(attrs))
    # mandatory constructor from ModelParameters that allows for automated reconstruction
    SpeedyParam(nt::NamedTuple) = new{typeof(nt.val)}(nt.val, nt.bounds, nt.desc, _attrs(nt))
end

# Mandatory ModelParameters `AbstractParam` interface methods
ModelParameters.parent(param::SpeedyParam) = (val=value(param), bounds=bounds(param), desc=description(param), attributes(param)...)
ModelParameters.rebuild(param::SpeedyParam, newvalues) = SpeedyParam(newvalues)

"""
Returns the value of the given `param`.
"""
value(param::SpeedyParam) = getfield(param, :value)

"""
Returns the numerical bounds of the parameter `param`.
"""
bounds(param::SpeedyParam) = getfield(param, :bounds)

"""
Returns a human-readable description of the parameter `param`.
"""
description(param::SpeedyParam) = getfield(param, :desc)

"""
Returns a `NamedTuple` of any additional attributes defined for `param`.
"""
attributes(param::SpeedyParam) = getfield(param, :attrs)

# parameterize

"""
    $SIGNATURES

Returns a "parameterized" copy of the given data structure. For numeric types, this simply wraps
the given value in a `Param`-like type. Dispatches can be added for compound data types that
reconstruct the data type with `Param`-like wrappers on all relevant parameters.
"""
parameters(obj; kwargs...) = parameters(SpeedyParam, obj; kwargs...)
parameters(::Type{PT}, obj; kwargs...) where {PT<:ModelParameters.AbstractParam} = (;)
parameters(::Type{PT}, param::ModelParameters.AbstractParam; kwargs...) where {PT<:ModelParameters.AbstractParam} = PT(merge(parent(param), kwargs))
parameters(::Type{PT}, x::Union{Number,AbstractArray}; kwargs...) where {PT<:ModelParameters.AbstractParam} = PT(x; kwargs...)

"""
    $SIGNATURES

Reconstructs the given data structure with the given `values`. If `values` is a `NamedTuple`, the structure
must match that of `obj`. This function is used to reconstruct model types from `SpeedyParams`.
"""
reconstruct(obj::T, value::T) where {T} = value
reconstruct(obj::ModelParameters.AbstractParam, value::T) where {T} = ModelParameters.update(obj, Tuple(value))
reconstruct(obj::NamedTuple{keys,V}, values::NamedTuple{keys,V}) where {keys,V<:Tuple} = values
reconstruct(obj, values::ComponentArray) = reconstruct(obj, NamedTuple(values))
function reconstruct(obj, values::NamedTuple{keys}) where {keys}
    # extract properties from obj as named tuple
    nt = getproperties(obj)
    # recursively call reconstruct for all keys specified in values
    patchvals = map(keys) do key
        reconstruct(nt[key], values[key])
    end
    # construct a named tuple with the reconstructed values and apply with setproperties
    patch = NamedTuple{keys}(patchvals)
    return setproperties(obj, patch)
end

# internal helper function for filtering out built-in parameter fields
_attrs(nt::NamedTuple) = nt[filter(k -> k ∉ (:val,:bounds,:desc), keys(nt))]

