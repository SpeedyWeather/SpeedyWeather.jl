export PrognosticVariable,
    GridVariable,
    TendencyVariable,
    DynamicsVariable,
    ParameterizationVariable,
    ParticleVariable,
    ScratchVariable

for Var in (
        :PrognosticVariable,            # variable groups to define
        :GridVariable,
        :TendencyVariable,
        :DynamicsVariable,
        :ParameterizationVariable,
        :ParticleVariable,
        :ScratchVariable,
    )

    @eval begin
        """$(TYPEDSIGNATURES) A variable defined through its `name`, dimensions `dims` and optionally
        its `units`, description `desc` and the `namespace` it's sorted under within a
        variable group. This definition is only used to return an object to define that a given
        model component wants to define a variable -- allocation and being sorted into the
        `Variables` tree happens elsewhere."""
        @kwdef struct $Var{D, U} <: AbstractVariable{D}
            name::Symbol
            dims::D
            units::U = ""
            desc::String = ""
            namespace::Symbol = Symbol()    # empty symbol is used for atmosphere
        end

        nonparametric_type(v::$Var) = nonparametric_type(typeof(v))
        nonparametric_type(::Type{<:$Var}) = $Var
        $Var(name, dims; kwargs...) = $Var(; name, dims, kwargs...)
        $Var(name, dims, units; kwargs...) = $Var(; name, dims, units, kwargs...)
        $Var(name, dims, units, desc; kwargs...) = $Var(; name, dims, units, desc, kwargs...)
        $Var(name, dims, units, desc, namespace) = $Var(; name, dims, units, desc, namespace)
    end
end

function Base.show(io::IO, var::AbstractVariable)
    print(
        io, nonparametric_type(var),
        "(:", var.name, ", ", var.dims, ", units = ", var.units,
        ", desc = ", var.desc, ", namespace = ", var.namespace, ")"
    )
    return nothing
end

export Variables
"""$(TYPEDSIGNATURES)
Symbol representing a unique identifier for a variable, combining its type, namespace and name and namespace. 
Used to remove duplicates when extracting variables from the model. E.g. `:Variable_ocean_sea_surface_temperature`."""
identifier(v::AbstractVariable) = Symbol(v.namespace, :_, v.name)

export Variables

"""Struct that holds all variables of a simulation. Variables are split into 
groups corresponding to the fields here $(TYPEDFIELDS) Each group can have
their own namespaces to distinguish between e.g. ocean, land, or tracer variables.
All non-prognostic groups are considered to be diagnostic with no memory between time steps."""
@kwdef struct Variables{Po, G, T, D, Pm, Pt, S} <: AbstractVariables
    "Prognostic variables subject to time stepping."
    prognostic::Po = NamedTuple()

    "Variables defined on the grid, mostly copies of spectral prognostic variables."
    grid::G = NamedTuple()

    "Tendencies of the prognostic variables"
    tendencies::T = NamedTuple()

    "Variables used in the dynamical core"
    dynamics::D = NamedTuple()

    "Variables used in the parameterizations"
    parameterizations::Pm = NamedTuple()

    "Variables used for particle advection"
    particles::Pt = NamedTuple()

    "Scratch variables for temporary storage during calculations with undetermined state (write before read)."
    scratch::S = NamedTuple()
end

Adapt.@adapt_structure Variables

"""$(TYPEDSIGNATURES)
Copy all entries from `src` to `dest` by recursing over the variable groups
and namespaces. Uses `copyto!` for arrays, `copy!` for `Clock`,
and direct assignment for `Ref` values."""
function Base.copy!(dest::Variables, src::Variables)
    for group in propertynames(dest)
        copy_variables!(getfield(dest, group), getfield(src, group))
    end
    return dest
end

"""$(TYPEDSIGNATURES)
Copy variables in `NamedTuples` from `src` into `dest`, only copying keys present in both.
Recurses into nested NamedTuples (namespaces like ocean, land, tracers).
Uses `@generated` to unroll the key iteration at compile time so that Enzyme
can differentiate through this function without runtime reflection more easily."""
@generated function copy_variables!(dest::NamedTuple{K1}, src::NamedTuple{K2}) where {K1, K2}
    exprs = Expr[]
    for key in K1
        key in K2 || continue
        push!(exprs, :(_copy_entry!(getfield(dest, $(QuoteNode(key))), getfield(src, $(QuoteNode(key))))))
    end
    return quote
        $(exprs...)
        return nothing
    end
end

_copy_entry!(dest::AbstractArray, src::AbstractArray) = copyto!(dest, src)
_copy_entry!(dest::NamedTuple, src::NamedTuple) = copy_variables!(dest, src)
_copy_entry!(dest::Base.RefValue, src::Base.RefValue) = (dest[] = src[])

# fallback for mutable structs: copy each field via setfield!
# and immutable structs with array fields (e.g. ScratchMemory): copyto! on arrays.
# uses @generated to avoid runtime reflection (ismutable, fieldnames, isa checks)
# which Enzyme cannot differentiate through.
@generated function _copy_entry!(dest::T, src::T) where {T}
    exprs = Expr[]
    for (i, fname) in enumerate(fieldnames(T))
        ft = fieldtype(T, i)
        if ft <: AbstractArray
            push!(exprs, :(copyto!(getfield(dest, $(QuoteNode(fname))), getfield(src, $(QuoteNode(fname))))))
        elseif ft <: NamedTuple
            push!(exprs, :(copy_variables!(getfield(dest, $(QuoteNode(fname))), getfield(src, $(QuoteNode(fname))))))
        elseif ismutabletype(T)
            push!(exprs, :(setfield!(dest, $(QuoteNode(fname)), getfield(src, $(QuoteNode(fname))))))
        end
        # for immutable structs: skip non-array, non-NamedTuple fields (not copyable in place)
    end
    return quote
        $(exprs...)
        return nothing
    end
end

# pretty printing
function Base.show(io::IO, V::Variables)
    Vsize = prettymemory(Base.summarysize(V))
    print(io, styled"{warning:Variables}", "{@NamedTuple{...}, ...} ", styled"{note:($Vsize)}")
    for (i, p) in enumerate(propertynames(V))
        lasti = i == length(propertynames(V))           # check if last property to choose ending └
        s = lasti ? "└" : "├"                           # choose ending
        psize = prettymemory(Base.summarysize(getfield(V, p)))
        print(io, "\n$s", styled"{info: $p }", styled"{note:($psize)}")
        for (j, k) in enumerate(keys(getfield(V, p)))
            lastj = j == length(keys(getfield(V, p)))   # check if last variable in namespace to choose ending └
            s2 = lastj ? "└" : "├"                      # choose ending └ for last variable
            maybe_bar1 = lasti ? " " : "│"              # if last variable in namespace, no vertical bar needed
            nt = getfield(getfield(V, p), k)
            if nt isa NamedTuple
                print(io, "\n$maybe_bar1 $s2 ", styled"{success:$k}")
                for (l, m) in enumerate(keys(nt))
                    s3 = l == length(keys(nt)) ? "└" : "├"  # choose ending └ for last variable
                    maybe_bar2 = lastj ? " " : "│"          # if last variable in namespace, no vertical bar needed
                    smry = Base.summary(getfield(nt, m))
                    line = "$maybe_bar1 $maybe_bar2 $s3 " * styled"{magenta:$m}: $smry"
                    line_short = textwidth(line) > 75 ? first(line, 75) * "..." : line
                    print(io, "\n", line_short)
                end
            else
                smry = Base.summary(nt)
                line = "$maybe_bar1 $s2 " * styled"{magenta:$k}: $smry"
                line_short = textwidth(line) > 75 ? first(line, 75) * "..." : line
                print(io, "\n", line_short)
            end
        end
    end
    return nothing
end

# runic: off
"""$(TYPEDSIGNATURES) Allocate all variables for a `model` as defined by its components.
Filters out duplicates and sorts the variables into groups and namespaces."""
function Variables(model::AbstractModel)
    all_vars = all_variables(model)     # one long tuple for all required variables of model and its components
    prognostic        = allocate(filter_variables(all_vars,       PrognosticVariable), model)
    grid              = allocate(filter_variables(all_vars,             GridVariable), model)
    tendencies        = allocate(filter_variables(all_vars,         TendencyVariable), model)
    dynamics          = allocate(filter_variables(all_vars,         DynamicsVariable), model)
    parameterizations = allocate(filter_variables(all_vars, ParameterizationVariable), model)
    particles         = allocate(filter_variables(all_vars,         ParticleVariable), model)
    scratch           = allocate(filter_variables(all_vars,          ScratchVariable), model)
    return Variables(; prognostic, grid, tendencies, dynamics, parameterizations, particles, scratch)
end
# runic: on

"""$(TYPEDSIGNATURES) Allocates all variables within a group given a tuple of variables
expected to be `<: AbstractVariable` definitions. Determines the namespaces,
allocates the arrays with zeros and collects them into NamedTuples."""
function allocate(group, model)
    length(group) == 0 && return NamedTuple()  # return empty NamedTuple if no variables to initialize
    namespaces = filter(k -> k != Symbol(), tuple(keys(group)...))

    # variables without namespace identified by empty symbol Symbol() go directly into the main NamedTuple
    # that way we have variables.prognostic.vorticity skipping the namespace between prognostic and vor
    nt1 = NamedTuple{Tuple(map(v -> v.name, group[Symbol()]))}(Tuple(map(var -> zero(var, model), group[Symbol()])))

    # other variables grouped by namespace
    # e.g. variables.prognostic.ocean.sea_surface_temperature, variables.prognostic.land.soil_moisture, etc.
    nt2 = NamedTuple{namespaces}(
        Tuple(
            map(
                ns ->
                NamedTuple{Tuple(map(v -> v.name, group[ns]))}(Tuple(map(var -> zero(var, model), group[ns])))
                , namespaces
            )
        )
    )

    return merge(nt1, nt2)
end

variables(::Nothing) = ()                                   # to allow for model.component = nothing
variables(::Any) = ()                                       # fallback for any component

# fallbacks in case model components are Nothing
initialize!(::Variables, ::Nothing, ::Any) = nothing        # used to define initial conditions
timestep!(::Variables, ::Nothing, ::Any) = nothing

"""$(TYPEDSIGNATURES)
Define variables needed by a model component. This is used to identify all required variables
by a model and initialize them. The default is to return an empty tuple, but components can
define their variables by extending this function for their type."""
variables

"""$(TYPEDSIGNATURES)
Fallback: component can define `variables(::Component, ::Model)` or simply `variables(::Component)`.
In the former, `model` is available to define required variables based on other model components,
in the latter only the component itself determines which variables are needed."""
variables(component, model) = variables(component)

"""$(TYPEDSIGNATURES)
Extracts all variables from the model by iterating over all components and collecting their variables.
Return a tuple of all variables."""
function all_variables(model::AbstractModel)
    t = variables(model)                        # variables from the model itself
    for component in propertynames(model)       # iterate over all components of the model
        # pass on model as well to allow for cross-component information to be used to define required variables
        vars = variables(getproperty(model, component), model)
        if length(vars) > 0
            t = tuple(t..., vars...)
        end
    end
    return t
end

get_namespaces(variables::AbstractVariable...) = unique([v.namespace for v in variables])

"""$(TYPEDSIGNATURES) Filters a tuple of variables `vars` by `VariableTyple`,
removes duplicates such that the variable path remains unique. Returns a
dictionary of the variable definitions (but the arrays aren't allocated here but in `allocate`)"""
function filter_variables(vars, VariableType)
    vars = filter(v -> v isa VariableType, vars)    # filter by variable type
    vars = unique(v -> identifier(v), vars)         # remove duplicates by identifier (group+namespace+name)

    namespaces = get_namespaces(vars...)            # Split by namespace
    group = Dict{Symbol, Tuple}()                   # assemble the variable group
    for ns in namespaces
        group[ns] = Tuple([v for v in vars if v.namespace == ns])
    end

    return group
end

"""$(TYPEDSIGNATURES)
Helper function to warn when a variable is not defined in the Variables. Reports on
group (prognostic etc.) and namespace therein (e.g. ocean, land, etc.) and lists the defined variables in that group and namespace.
Returns true to allow for short-circuiting with `&& return nothing` to exit a function early."""
function warn_undefvar(vars::Variables, key::Symbol, group::Symbol = :prognostic, namespace::Symbol = Symbol())
    path = namespace == Symbol() ? "$group" : "$group.$namespace"
    defined_vars = namespace == Symbol() ? keys(getfield(vars, group)) : keys(getfield(getfield(vars, group), namespace))
    @warn "Variable $key not defined in variables.$path. Defined are: $defined_vars"
    return true # return true to allow for short-circuiting with && return nothing to skip exit the following code early
end

# TODO move get_step, get_steps to LowerTriangularArrays?

function get_steps(coeffs::LowerTriangularArray{T, 2}) where {T}
    nsteps = size(coeffs, 2)
    return ntuple(i -> lta_view(coeffs, :, i), nsteps)
end

function get_steps(coeffs::LowerTriangularArray{T, 3}) where {T}
    nsteps = size(coeffs, 3)
    return ntuple(i -> lta_view(coeffs, :, :, i), nsteps)
end

export get_step

"""$(TYPEDSIGNATURES)
Get the i-th step of a LowerTriangularArray as a view (wrapped into a LowerTriangularArray).
"step" refers to the last dimension, for prognostic variables used for the leapfrog time step.
This method is for a 2D spectral variable (horizontal only) with steps in the 3rd dimension."""
get_step(coeffs::LowerTriangularArray{T, 2}, i) where {T} = lta_view(coeffs, :, i)

"""$(TYPEDSIGNATURES)
Get the i-th step of a LowerTriangularArray as a view (wrapped into a LowerTriangularArray).
"step" refers to the last dimension, for prognostic variables used for the leapfrog time step.
This method is for a 3D spectral variable (horizontal+vertical) with steps in the 4rd dimension."""
get_step(coeffs::LowerTriangularArray{T, 3}, i) where {T} = lta_view(coeffs, :, :, i)
