"""
    SpeedyParam{NF<:AbstractFloat} <: AbstractParam

Specialized implementation of `AbstractParam` for Speedy models that defines a fixed
set of fields with concrete types.
"""
Base.@kwdef struct SpeedyParam{NF<:AbstractFloat} <: AbstractParam{NF}
    "numeric value of the parameter"
    value::NF = NaN
    
    "numerical domain on which the parameter is defined"
    bounds::Domain = ℝ

    "human-readable description of the parameter"
    desc::String = ""

    "additional attributes with unspecified types"
    attrs::NamedTuple

    # default constructor
    SpeedyParam(value::NF, bounds::Domain, desc::String, attrs::NamedTuple) where {NF<:AbstractFloat} = new{NF}(value, bounds, desc, attrs)
    # convenience constructor
    SpeedyParam(value::NF; bounds=ℝ, desc="", attrs...) where {NF<:AbstractFloat} = new{NF}(value, bounds, desc, NamedTuple(attrs))
    # mandatory constructor from ModelParameters that allows for automated reconstruction
    SpeedyParam(nt::NamedTuple) = new{typeof(nt.val)}(nt.val, nt.bounds, nt.desc, _attrs(nt))
end

# Mandatory ModelParameters `AbstractParam` interface methods
ModelParameters.parent(param::SpeedyParam) = (val=value(param), attributes(param)..., bounds=bounds(param), desc=description(param))
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

"""
    $SIGNATURES

Returns a "parameterized" copy of the given data structure. For numeric types, this simply wraps
the given value in a `Param`-like type. Dispatches can be added for compound data types that
reconstruct the data type with `Param`-like wrappers on all relevant parameters.
"""
parameters(obj; kwargs...) = parameters(SpeedyParam, obj; kwargs...)
parameters(::Type{PT}, obj; kwargs...) where {PT<:AbstractParam} = (;)
parameters(::Type{PT}, param::AbstractParam; kwargs...) where {PT<:AbstractParam} = PT(merge(parent(param), kwargs))
parameters(::Type{PT}, x::Union{Number,AbstractArray}; kwargs...) where {PT<:AbstractParam} = PT(x; kwargs...)

"""
    $SIGNATURES

Convenience method that creates a model parameter from its property with the given `name` and optional extra attributes in `kwargs`.
A parameter attribute `copmonenttype` is automatically added with value `T`.
"""
parameterof(obj::T, name::Symbol; kwargs...) where {T} = parameterof(SpeedyParam, obj, name; kwargs...)
parameterof(::Type{PT}, obj::T, name::Symbol; kwargs...) where {PT,T} = name => parameters(PT, getproperty(obj, name); componenttype=T, kwargs...)

"""
    SpeedyParams{NT<:NamedTuple} <: ModelParameters.AbstractModel

Lightweight wrapper around a `NamedTuple` of parameters following the ModelParameters `AbstractModel` interface.
This provides a table-like interface for interacting with model parameters.
"""
struct SpeedyParams{NT<:NamedTuple} <: ModelParameters.AbstractModel
    parent::NT
    # This constructor is adapted from ModelParameters.jl
    function SpeedyParams(parent::NamedTuple)
        if ModelParameters.hasparam(parent)
            # add component type and filename attributes to parameters
            expandedpars = ModelParameters._expandkeys(ModelParameters.params(parent))
            parent = ModelParameters.Flatten.reconstruct(parent, expandedpars, ModelParameters.SELECT, ModelParameters.IGNORE)
        else
            @warn "Model has no parameter fields"
        end
        return new{typeof(parent)}(parent)
    end
    SpeedyParams(; params...) = SpeedyParams((; params...))
end

# Base overridese for SpeedyParams

## parameter subsets
Base.getindex(ps::SpeedyParams, param_label::String) = getindex(ps, [param_label])
@inline function Base.getindex(ps::SpeedyParams, param_labels::Vector{String})
    # extract labels from ComponentVector and 
    ls = labels(vec(ps))
    idx = map(i -> findfirst(==(i), ls), param_labels)
    # check if any names were not found
    unmatched_idx = findfirst(isnothing, idx)
    @assert isnothing(unmatched_idx) "$(param_labels[unmatched_idx]) is not a valid parameter index"
    # get flattened parameter tuple
    params = ModelParameters.params(ps)
    # get tuple of selected parameters
    selected_params = params[idx]
    # extract parameters from nested named tuple and reconstruct SpeedyParams
    return SpeedyParams(_selectrecursive(∈(selected_params), parent(ps)))
end

## selecting columns
Base.getindex(ps::SpeedyParams, nm::Symbol) = getindex(ps, :, nm)
@inline function Base.getindex(ps::SpeedyParams, ::Colon, nm::Symbol)
    # include fieldname only (no component)
    if nm == :idx
        1:length(ps)
    elseif nm == :fieldname
        ModelParameters.paramfieldnames(ps)
    else
        map(p -> getindex(p, nm), ModelParameters.params(ps))
    end
end

## non-mutating setindex!
## consider using a method instead of setindex! to avoid confusion?
Base.setindex!(ps::SpeedyParams, x, nm::Symbol) = setindex!(ps, x, :, nm)
@inline function Base.setindex!(ps::SpeedyParams, x, ::Colon, nm::Symbol)
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
        return SpeedyParams(newparent)
    end
end

Base.keys(params::SpeedyParams) = (:idx, :fieldname, keys(first(params))...)
Base.vec(params::SpeedyParams) = ComponentArray(parent(params))

# Override internal ModelParameters method _columntypes to condense type names in table schema (just to look nicer)
# TODO: propose this change upstream in ModelParameters.jl
ModelParameters._columntypes(ps::SpeedyParams) = map(k -> promote_type(map(typeof, getindex(ps, k))...), keys(ps)) 

# reconstruct

"""
    $SIGNATURES

Reconstructs the given data structure with the given `values`. If `values` is a `NamedTuple`, the structure
must match that of `obj`. This function is used to reconstruct model types from `SpeedyParams`.
"""
reconstruct(obj::T, value::T) where {T} = value
reconstruct(obj::AbstractParam, value::T) where {T} = ModelParameters.update(obj, Tuple(value))
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

"""
    $SIGNATURES

Convenience macro for creating "parameterized" `struct`s. Fields can be marked as parameters with `@param`, e.g:

```julia
@parameterized Base.@kwdef struct Foo{T}
    @param x::T = 1.0
    y::T = NaN
end
```
where here `x` will be treated as a parameter field and `y` will not. Parameters will be then automatically
generated in a corresponding `parameters(::Foo)` method definition. Additional parameter attributes can be
supplied as keywords after the parameter, e.g. `@param p::T = 0.5 bounds=UnitInterval()` or
`@param p::T = 0.5 (bounds=UnitInterval(), constant=false)`.
"""
macro parameterized(expr)
    # process struct definition
    ## first declare local variables for capturing metadata
    paraminfo = []
    lastdoc = nothing
    typename = nothing
    ## prewalk (top-down) over expression tree; note that there is some danger of infinite
    ## recursion with prewalk since it will descend into the returned expression. However,
    ## this is not a problem here since we only deconstruct/reduce the @param expressions.
    new_expr = MacroTools.prewalk(expr) do ex
        # case 1: struct definition (top level)
        if isa(ex, Expr) && ex.head == :struct
            # type signature is in second argument;
            # use namify to extract type name (discard type parameters)
            typename = namify(ex.args[2])
            ex
        # case 2: parameterized field definition
        elseif @capture(ex, @param fieldname_::FT_ = defval_; attrs__) || @capture(ex, @param fieldname_::FT_; attrs__)
            # use last seen docstring as description, if present
            desc = isnothing(lastdoc) ? "" : lastdoc
            attrs = if isempty(attrs)
                (;)
            elseif length(attrs) == 1
                # handle both singleton attr=value syntax as well as named tuple syntax (attr1=value1, attr2=value2, ...)
                attrs[1].head == :tuple ? eval(attrs[1]) : eval(:(($(attrs...),)))
            else
                error("invalid syntax for parameter attributes: $(QuoteNode(attrs))")
            end
            # then reset to nothing to prevent duplication
            lastdoc = nothing
            push!(paraminfo, (name=fieldname, type=FT, desc=desc, attrs...))
            :($fieldname::$FT = $defval)
        # case 3: non-parameterized field definition
        elseif @capture(ex, fieldname_::FT_ = defval_) || @capture(ex, fieldname_::FT_)
            # reset lastdoc variable to prevent confusing docstrings between parameteter and non-parameter fields
            lastdoc = nothing
            ex
        # case 4: field documentation
        elseif @capture(ex, docs_String)
            # store in lastdoc field
            lastdoc = docs
            ex
        # catch-all case for all other nodes in the epxression tree
        else
            ex
        end
    end
    # emit parametersof calls for each parsed parameter
    param_constructors = map(paraminfo) do info
        :($(QuoteNode(info.name)) => parameterof(obj, $(QuoteNode(info.name)); desc=$(info.desc), info.attrs..., kwargs...))
    end
    # construct final expression block
    block = Expr(:block)
    ## 1. struct definition
    push!(block.args, esc(new_expr))
    ## 2. parameters method dispatch
    push!(block.args, esc(:(parameters(obj::$(typename); kwargs...) = (; $(param_constructors...)))))
    return block
end

# Internal helper functions

"""
Internal helper function for filtering out built-in parameter fields.
"""
_attrs(nt::NamedTuple) = nt[filter(k -> k ∉ (:val,:bounds,:desc), keys(nt))]

"""
    _selectrecursive(selector, nt::NamedTuple)

Recursively filters out all values from a (possibly nested) named tuple `nt` for which `selector`
returns true.
"""
_selectrecursive(selector, x) = selector(x) ? x : nothing
function _selectrecursive(selector, nt::NamedTuple)
    # recursively apply selector
    new_nt = map(x -> _selectrecursive(selector, x), nt)
    # filter out all nothing values
    selected_keys = filter(k -> !isnothing(new_nt[k]), keys(new_nt))
    # construct filtered named tuple
    return NamedTuple{selected_keys}(map(k -> new_nt[k], selected_keys))
end
