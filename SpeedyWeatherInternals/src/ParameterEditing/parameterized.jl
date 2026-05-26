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
e.g. `@param p::T = 0.5 bounds=UnitInterval` or `@param p::T = 0.5 (bounds=UnitInterval, constant=false)`.

!!! warning "Known limitations"
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

    If you encounter any other weird interactions or behaviors, please raise [an issue](https://github.com/SpeedyWeather/SpeedyWeather.jl/issues/new).
"""
macro parameterized(expr)
    function typedef2sig(typedef)
        return if MacroTools.@capture(typedef, T_ <: super__)
            T
        else
            typedef
        end
    end
    function parse_attributes(attrs)
        return if isempty(attrs)
            (;)
        elseif length(attrs) == 1
            # handle both singleton attr=value syntax as well as named tuple syntax (attr1=value1, attr2=value2, ...)
            attrs[1].head == :tuple ? :($(attrs[1])) : :(($(attrs...),))
        else
            error("invalid syntax for parameter attributes: $(QuoteNode(attrs))")
        end
    end
    # process struct definition
    ## first declare local variables for capturing metadata
    params = []
    lastdoc = nothing
    typedef = nothing
    has_kwdef = false
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
            attrs = parse_attributes(attrs)
            paraminfo = (name = fieldname, type = FT, desc = desc, attrs = attrs)
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
            # case 3: subcomponent
        elseif MacroTools.@capture(ex, @component fieldname_::FT_ = defval_ attrs__) || MacroTools.@capture(ex, @component fieldname_::FT_ attrs__)
            attrs = parse_attributes(attrs)
            # add component name to attributes
            paraminfo = (name = fieldname, type = FT, desc = "", attrs = (component = fieldname, attrs...))
            push!(params, paraminfo)
            # reset lastdoc
            lastdoc = nothing
            isnothing(defval) ? :($fieldname::$FT;) : :($fieldname::$FT = $defval;)
            # case 4: field documentation
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
        :($(QuoteNode(info.name)) => ParameterEditing.parameterof(obj, Val{$(QuoteNode(info.name))}(); desc = $(info.desc), $(info.attrs)..., kwargs...))
    end
    # construct final expression block
    block = Expr(:block)
    struct_block = Expr(:block)
    ## 1. struct definition
    # emit the interpolation placeholder $(Expr(:meta, :doc)) so Julia's doc system
    # will attach any source docstring to the generated struct (same as Base.@kwdef);
    # for some reason, we need to use a dedicated block for the struct def otherwise it won't work
    if !has_kwdef
        push!(struct_block.args, Expr(:meta, :doc))
    end
    push!(struct_block.args, esc(new_expr))
    push!(block.args, struct_block)
    ## 2. parameters method dispatch
    parameters_func = quote
        function ParameterEditing.parameters(obj::$(typename); kwargs...)
            param_nt = (; $(param_constructors...))
            nonempty_keys = filter(name -> !isempty(param_nt[name]), keys(param_nt))
            return ParameterEditing.ParameterTable(NamedTuple{nonempty_keys}(map(name -> param_nt[name], nonempty_keys)))
        end
    end
    push!(block.args, esc(parameters_func))
    ## 3. override ConstructionBase.setproperties if @kwdef is used; we do this so that internal logic in the keyword defaults gets preserved
    if has_kwdef && typesig.head == :curly
        # handle type arguments; first extract argument names (discarding upper type bounds)
        typeargs = map(typedef2sig, typesig.args[2:end])
        # then build method signature - merge original properties with patch to preserve non-parameter fields
        push!(block.args, esc(:(ParameterEditing.setproperties(obj::$(typename){$(typeargs...)}, patch::NamedTuple) where {$(typesig.args[2:end]...)} = $(typename){$(typeargs...)}(; merge(ParameterEditing.getproperties(obj), patch)...))))
    elseif has_kwdef
        # otherwise we can just use the typename
        push!(block.args, esc(:(ParameterEditing.setproperties(obj::$(typename), patch::NamedTuple) = $(typename)(; merge(ParameterEditing.getproperties(obj), patch)...))))
    end
    return block
end
