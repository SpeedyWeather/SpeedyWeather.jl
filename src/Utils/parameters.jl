"""
    SpeedyParam{NF<:AbstractFloat} <: AbstractParam

Specialized implementation of `AbstractParam` for Speedy models that defines a fixed
set of fields with concrete types.
"""
Base.@kwdef struct SpeedyParam{NF<:AbstractFloat} <: AbstractParam{NF}
    "numeric value of the parameter"
    value::NF = NaN
    
    "numerical domain on which the parameter is defined"
    bounds::Domain = Unbounded

    "human-readable description of the parameter"
    desc::String = ""

    "additional attributes with unspecified types"
    attrs::NamedTuple

    # default constructor
    SpeedyParam(value::NF, bounds::Domain, desc::String, attrs::NamedTuple) where {NF<:AbstractFloat} = new{NF}(value, bounds, desc, attrs)
    # convenience constructor
    SpeedyParam(value::NF; bounds=Unbounded, desc="", attrs...) where {NF<:AbstractFloat} = new{NF}(value, bounds, desc, NamedTuple(attrs))
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
    SpeedyParams{NT<:NamedTuple} <: ModelParameters.AbstractModel

Lightweight wrapper around a `NamedTuple` of parameters following the ModelParameters `AbstractModel` interface.
This provides a table-like interface for interacting with model parameters.
"""
struct SpeedyParams{NT<:NamedTuple} <: ModelParameters.AbstractModel
    parent::NT
    SpeedyParams(parent::NamedTuple) = new{typeof(parent)}(parent)
    SpeedyParams(; params...) = SpeedyParams((; params...))
end

"""
    $SIGNATURES

Simple recursive dispatch algorithm for unpacking nested named tuples of `SpeedyParams`
"""
unpack_params(x) = x
unpack_params(nt::NamedTuple) = map(unpack_params, nt)
unpack_params(params::SpeedyParams) = unpack_params(parent(params))

"""
    $SIGNATURES

Unpacks and strips all `AbstractParam` types from the given `params` object, replacing them with their numeric `val`.
This is a type-stable alternative to the (currently unstable) `stripparams` in `ModelParameters`.
"""
stripparams(param::AbstractParam) = value(param)
stripparams(params::NamedTuple) = map(stripparams, params)
stripparams(params::SpeedyParams) = stripparams(unpack_params(params))

# Base overrides for SpeedyParams

## parameter subsets
Base.getindex(ps::SpeedyParams, param_label::String) = getindex(ps, [param_label])
@inline function Base.getindex(ps::SpeedyParams, param_labels::Vector{String})
    # extract labels from ComponentVector 
    ls = labels(vec(ps))
    idx = reduce(vcat, map(query -> findall(key -> startswith(key, query), ls), param_labels))
    # check if any names were not found
    length(idx) > 0 || @warn "no parameters matched to labels: $(param_labels)"
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

function Base.vec(params::SpeedyParams)
    # recursively unpack params and build component vector
    return ComponentArray(unpack_params(params))
end

# Override internal ModelParameters method _columntypes to condense type names in table schema (just to look nicer)
# TODO: propose this change upstream in ModelParameters.jl
ModelParameters._columntypes(ps::SpeedyParams) = map(k -> promote_type(map(typeof, getindex(ps, k))...), keys(ps)) 

# parameters method interface
"""
    $SIGNATURES

Extract parameters from the given `obj` as (possibly nested) named-tuple of `SpeedyParam`s or some other
`AbstractParam` type. If `obj`
"""
parameters(obj; kwargs...) = (;)
parameters(param::PT; kwargs...) where {PT<:AbstractParam} = parameters(PT, param; kwargs...)
parameters(param::Union{Number,AbstractArray}; kwargs...) = parameters(SpeedyParam, param; kwargs...)
parameters(::Type{PT}, obj; kwargs...) where {PT<:AbstractParam} = parameters(obj; kwargs...)
parameters(::Type{PT}, param::AbstractParam; kwargs...) where {PT<:AbstractParam} = PT(merge(parent(param), kwargs))
parameters(::Type{PT}, x::Union{Number,AbstractArray}; kwargs...) where {PT<:AbstractParam} = PT(x; kwargs...)

"""
    $SIGNATURES

Convenience method that creates a model parameter from its property with the given `name` and optional extra attributes in `kwargs`.
A parameter attribute `copmonenttype` is automatically added with value `T`.
"""
parameterof(obj::T, ::Val{propname}; kwargs...) where {T,propname} = parameterof(SpeedyParam, obj, Val{propname}(); kwargs...)
parameterof(::Type{PT}, obj::T, ::Val{propname}; kwargs...) where {PT<:AbstractParam,T,propname} = parameters(PT, getproperty(obj, propname); merge((; kwargs...), (componentttype=T,),)...)

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
reconstruct(obj, values::SpeedyParams) = reconstruct(obj, stripparams(values))
@generated function reconstruct(obj, values::NamedTuple{keys}) where {keys}
    recursive_calls = map(k -> :(reconstruct(obj.$k, values.$k)), keys)
    quote
        # recursively call reconstruct for all keys specified in values
        patchvals = tuple($(recursive_calls...))
        # construct a named tuple with the reconstructed values and apply with setproperties
        patch = NamedTuple{keys}(patchvals)
        return setproperties(obj, patch)
    end
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
generated in a corresponding `parameters(::Foo)` method definition with descriptions set to the corresponding
struct field docstring, if present. Additional parameter attributes can be supplied as keywords after the parameter,
e.g. `@param p::T = 0.5 bounds=UnitInterval()` or `@param p::T = 0.5 (bounds=UnitInterval(), constant=false)`.

## Known limitations

This macro will fail to behave correctly in the following known corner cases:

    1. **Untyped struct fields without default values**. The macro will not work if the struct has untyped
    fields without default values specified using the `@kwdef` syntax. This is because untyped fields are
    not distinguishable in the expression tree from symbols in other expressions (such as typed field definitions).
    This problem can be avoided by always specifying types on struct fields (which should always be done anyway!).

    2. **`@kwdef` struct with a single parameter field and no default value**. `@parameterized` has a buggy
    interaction with `@kwdef` in this corner case:
    ```julia
    @parameterized @kwdef struct Foo{T}
        @param x::T
    end
    ```
    where the auto-generated constructors erroneously duplicate `T` as a function and type parameter. This is, however,
    a fairly meaningless corner-case since there is no purpose in marking this struct as `@kwdef` if no default assignments
    are made!

If you encounter any other weird interactions or behaviors, please raise an issue.
"""
macro parameterized(expr)
    function typedef2sig(typedef)
        if MacroTools.@capture(typedef, T_ <: super__)
            T
        else
            typedef
        end
    end
    # process struct definition
    ## first declare local variables for capturing metadata
    params = []
    lastdoc = nothing
    typedef = nothing
    has_kwdef = false
    # defval = nothing
    ## prewalk (top-down) over expression tree; note that there is some danger of infinite
    ## recursion with prewalk since it will descend into the returned expression. However,
    ## this is not a problem here since we only deconstruct/reduce the @param expressions.
    new_expr = MacroTools.prewalk(expr) do ex
        if MacroTools.@capture(ex, @kwdef structdef__) || MacroTools.@capture(ex, Base.@kwdef structdef__)
            has_kwdef = true
            ex
        # struct definition (top level)
        elseif isa(ex, Expr) && ex.head == :struct
            # extract type definition from second argument
            typedef = ex.args[2]
            ex
        # case 1: parameterized field definition
        elseif MacroTools.@capture(ex, @param fieldname_::FT_ = defval_ attrs__) || MacroTools.@capture(ex, @param fieldname_::FT_ attrs__)
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
            paraminfo = (name=fieldname, type=FT, desc=desc, attrs=attrs)
            push!(params, paraminfo)
            # reset lastdoc to nothing to prevent duplication
            lastdoc = nothing
            isnothing(defval) ? :($fieldname::$FT;) : :($fieldname::$FT = $defval;)
        # case 1a: untyped parameter
        elseif MacroTools.@capture(ex, @param fieldname_ = defval_ attrs__) || MacroTools.@capture(ex, @param fieldname_ attrs__)
            lastdoc = nothing
            @warn "ignoring untyped parameter $fieldname"
            isnothing(defval) ? :($fieldname;) : :($fieldname = $defval;)
        # case 2: non-parameter field
        elseif MacroTools.@capture(ex, fieldname_::FT_ = defval_) ||
            MacroTools.@capture(ex, fieldname_::FT_) ||
            MacroTools.@capture(ex, fieldname_ = defval_)
            # reset lastdoc variable to prevent confusing docstrings between parameteter and non-parameter fields
            lastdoc = nothing
            ex
        # case 3: field documentation
        elseif MacroTools.@capture(ex, docs_String)
            # store in lastdoc field
            lastdoc = docs
            ex
        # catch-all case for all other nodes in the epxression tree
        else
            ex
        end
    end
    # extract type info
    typename = MacroTools.namify(typedef)
    typesig = typedef2sig(typedef)
    # emit parameterof calls for each parsed parameter
    param_constructors = map(params) do info
        :($(QuoteNode(info.name)) => SpeedyWeather.parameterof(obj, Val{$(QuoteNode(info.name))}(); desc=$(info.desc), $(info.attrs)..., kwargs...))
    end
    # construct final expression block
    block = Expr(:block)
    ## 1. struct definition
    push!(block.args, esc(new_expr))
    ## 2. parameters method dispatch
    push!(block.args, esc(:(SpeedyWeather.parameters(obj::$(typename); kwargs...) = SpeedyWeather.SpeedyParams((; $(param_constructors...))))))
    ## 3. override ConstructionBase.setproperties if @kwdef is used; we do this so that internal logic in the keyword defaults gets preserved
    if has_kwdef && typesig.head == :curly
        # handle type arguments; first extract argument names (discarding upper type bounds)
        typeargs = map(typedef2sig, typesig.args[2:end])
        # then build method signature
        push!(block.args, esc(:(setproperties(obj::$(typename){$(typeargs...)}, patch::NamedTuple) where {$(typesig.args[2:end]...)} = $(typename){$(typeargs...)}(; patch...))))
    elseif has_kwdef
        # otherwise we can just use the typename
        push!(block.args, esc(:(setproperties(obj::$(typename), patch::NamedTuple) = $(typename)(; patch...))))
    end
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
_selectrecursive(selector, ps::SpeedyParams) = _selectrecursive(selector, parent(ps))
_selectrecursive(selector, x) = selector(x) ? x : nothing
function _selectrecursive(selector, nt::NamedTuple)
    # recursively apply selector
    new_nt = map(x -> _selectrecursive(selector, x), nt)
    # filter out all nothing values
    selected_keys = filter(k -> !isnothing(new_nt[k]), keys(new_nt))
    # construct named tuple from filtered keys or return nothing if no keys were selected
    return length(selected_keys) > 0 ? NamedTuple{selected_keys}(map(k -> new_nt[k], selected_keys)) : nothing
end
