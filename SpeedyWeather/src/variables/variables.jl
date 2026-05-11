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
        variable group. Variables sharing the same non-empty `fuse` symbol (within the same
        `namespace`) are allocated as views of a single shared parent buffer, which is exposed
        in the resulting NamedTuple under the `fuse` symbol. This definition is only used to
        return an object to define that a given model component wants to define a variable --
        allocation and being sorted into the `Variables` tree happens elsewhere."""
        @kwdef struct $Var{D, U} <: AbstractVariable{D}
            name::Symbol
            dims::D
            units::U = ""
            desc::String = ""
            namespace::Symbol = Symbol()    # empty symbol is used for atmosphere
            fuse::Symbol = Symbol()         # empty symbol = no fusion; same value = share parent buffer
        end

        nonparametric_type(v::$Var) = nonparametric_type(typeof(v))
        nonparametric_type(::Type{<:$Var}) = $Var
        $Var(name, dims; kwargs...) = $Var(; name, dims, kwargs...)
        $Var(name, dims, units; kwargs...) = $Var(; name, dims, units, kwargs...)
        $Var(name, dims, units, desc; kwargs...) = $Var(; name, dims, units, desc, kwargs...)
        $Var(name, dims, units, desc, namespace; kwargs...) = $Var(; name, dims, units, desc, namespace, kwargs...)
        $Var(name, dims, units, desc, namespace, fuse) = $Var(; name, dims, units, desc, namespace, fuse)
    end
end

function Base.show(io::IO, var::AbstractVariable)
    print(
        io, nonparametric_type(var),
        "(:", var.name, ", ", var.dims, ", units = ", var.units,
        ", desc = ", var.desc, ", namespace = ", var.namespace
    )
    var.fuse === Symbol() || print(io, ", fuse = :", var.fuse)
    print(io, ")")
    return nothing
end

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

# defined e.g. for output filters, as fieldnames(Variables) isn't fully type stable
const ALL_VARIABLE_GROUPS = (:prognostic, :grid, :tendencies, :dynamics, :parameterizations, :particles, :scratch)

Adapt.@adapt_structure Variables

"""$(TYPEDSIGNATURES)
Transfer all arrays in `vars` to `arch` using `on_architecture` for each group."""
function Architectures.on_architecture(arch::AbstractArchitecture, vars::Variables)
    return Variables(;
        (k => on_architecture(arch, getfield(vars, k)) for k in fieldnames(Variables))...
    )
end

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
Filters out duplicates and sorts the variables into groups and namespaces. Variables
sharing a non-empty `fuse` symbol (within the same namespace) share a single parent
buffer; their per-variable entries are views of that parent. Fused parents themselves
live under `vars.scratch.fused.<fuse_symbol>`."""
function Variables(model::AbstractModel)
    all_vars = all_variables(model)     # one long tuple for all required variables of model and its components

    # First pass (cross-type): allocate one parent per (namespace, fuse_symbol) group, spanning
    # all variable types. Returns a Dict keyed by (namespace, fuse) with values (parent, view-by-identifier).
    fuse_parents = build_fuse_parents(all_vars, model)

    prognostic        = allocate(filter_variables(all_vars,       PrognosticVariable), model, fuse_parents)
    grid              = allocate(filter_variables(all_vars,             GridVariable), model, fuse_parents)
    tendencies        = allocate(filter_variables(all_vars,         TendencyVariable), model, fuse_parents)
    dynamics          = allocate(filter_variables(all_vars,         DynamicsVariable), model, fuse_parents)
    parameterizations = allocate(filter_variables(all_vars, ParameterizationVariable), model, fuse_parents)
    particles         = allocate(filter_variables(all_vars,         ParticleVariable), model, fuse_parents)
    scratch           = allocate(filter_variables(all_vars,          ScratchVariable), model, fuse_parents)

    # Install fused parents at their canonical home: vars.scratch.fused.<fuse_symbol>.
    # Validates that fuse symbols are globally unique across namespaces.
    if !isempty(fuse_parents)
        fused = build_fused_namespace(fuse_parents)
        scratch = merge(scratch, (fused = fused,))
    end

    return Variables(; prognostic, grid, tendencies, dynamics, parameterizations, particles, scratch)
end
# runic: on

"""$(TYPEDSIGNATURES) Build the `vars.scratch.fused` NamedTuple from the `fuse_parents` Dict.
Keyed flat by fuse symbol; errors if the same fuse symbol is reused across namespaces."""
function build_fused_namespace(fuse_parents)
    seen = Dict{Symbol, Symbol}()  # fuse_symbol => namespace, for collision diagnostics
    pairs = Pair{Symbol, Any}[]
    for ((ns, fuse_sym), entry) in fuse_parents
        if haskey(seen, fuse_sym)
            error(
                "Fuse symbol `$fuse_sym` is used in multiple namespaces (`$(seen[fuse_sym])` and `$ns`). " *
                "Fuse symbols must be globally unique across namespaces so they can live flat under " *
                "`vars.scratch.fused.<fuse_symbol>`."
            )
        end
        seen[fuse_sym] = ns
        push!(pairs, fuse_sym => entry.parent)
    end
    return (; pairs...)
end

"""$(TYPEDSIGNATURES) For every non-empty `fuse` symbol used across `all_vars`, allocate one parent
buffer covering all members in that (namespace, fuse) group across all variable types.
Returns a `Dict{Tuple{Symbol,Symbol}, NamedTuple}` keyed by `(namespace, fuse)`. Each entry holds
the parent and a `Dict{Symbol, AbstractArray}` mapping each member's `identifier(v)` to its view.

Validates at build time that every constructed view actually aliases the parent at its declared
slot range, that slot ranges are pairwise non-overlapping, and that together they tile the parent's
fused axis exactly with no gaps. This catches offset/order bugs in `_split_views_*` and ensures
the right variable is wired to the right slot."""
function build_fuse_parents(all_vars, model)
    # Group all fused vars by (namespace, fuse), deduped by (namespace, name) Tuple.
    seen = Dict{Tuple{Symbol, Symbol}, AbstractVariable}()   # (namespace, name) => first variable seen
    groups = Dict{Tuple{Symbol, Symbol}, Vector{AbstractVariable}}()
    for v in all_vars
        v.fuse === Symbol() && continue
        key = (v.namespace, v.name)
        if haskey(seen, key)
            prev = seen[key]
            nonparametric_type(prev) === nonparametric_type(v) && continue  # same type, same var declared twice — fine
            error(
                "Fuse group `$(v.fuse)` (namespace `$(v.namespace)`): variable `$(v.name)` is " *
                "declared as both $(nonparametric_type(prev)) and $(nonparametric_type(v)). " *
                "Each (namespace, name) pair may appear at most once across all variable types in a fuse group."
            )
        end
        seen[key] = v
        push!(get!(groups, (v.namespace, v.fuse), AbstractVariable[]), v)
    end

    parents = Dict{Tuple{Symbol, Symbol}, NamedTuple}()
    for ((ns, fuse_sym), vars) in groups
        # Validate: same dim type within a fuse group (allocate_fused dispatches on it).
        D = typeof(first(vars).dims)
        for v in vars
            typeof(v.dims) === D || error(
                "Fuse group `$fuse_sym` (namespace `$ns`) mixes dim types: " *
                "$(typeof(v.dims)) for `$(v.name)` vs $D. All members must share the same dim type."
            )
        end
        # Validate: member names are unique within the group.
        names = Tuple(v.name for v in vars)
        length(unique(names)) == length(names) || error(
            "Fuse group `$fuse_sym` (namespace `$ns`) has duplicate member names: $names. " *
            "Each member of a fuse group must have a unique `name`."
        )
        # Allocate one parent for this fuse group; views & slots aligned with `vars`.
        parent, views, slots = allocate_fused(typed_vars(vars, D), model)
        # Build-time correctness check: every view aliases the parent at its declared slot
        # range, slots are pairwise disjoint, and they tile the parent's fused axis exactly.
        _validate_fuse_layout(fuse_sym, ns, vars, parent, views, slots)
        parents[(ns, fuse_sym)] = (
            parent = parent,
            views = Dict{Tuple{Symbol, Symbol}, Any}((v.namespace, v.name) => views[i] for (i, v) in enumerate(vars)),
        )
    end
    return parents
end

# Verify each view actually points at its declared slot inside the parent, slots don't
# overlap, and together they cover the parent's fused axis exactly. Discards the slot
# information after checking — it's a correctness check, not runtime metadata.
function _validate_fuse_layout(fuse_sym, ns, vars, parent, views, slots)
    parent_size = size(parent.data, 2)
    covered = falses(parent_size)
    for (i, v) in enumerate(vars)
        view_data = views[i].data
        slot = slots[i]
        if view_data isa SubArray
            Base.parent(view_data) === parent.data || error(
                "Fuse group `$fuse_sym` (namespace `$ns`): view for `$(v.name)` is not a view " *
                "of the fused parent buffer (got parent $(typeof(Base.parent(view_data))))."
            )
 
            p_index_last = parentindices(view_data)[end]
            # Normalise scalar indices (Grid2D / Spectral2D uses field_view(parent, :, k)
            # which stores an Int rather than a 1-element range in parentindices).
            p_index_last_range = p_index_last isa Integer ? (p_index_last:p_index_last) : p_index_last
            p_index_last_range == slot || error(
                "Fuse group `$fuse_sym` (namespace `$ns`): view for `$(v.name)` has " *
                "parentindices $(p_index_last) but slot map declares $(slot)."
            )
        else # fallback for GPU array types that are not SubArrays
            size(view_data, ndims(view_data)) == length(slot) || error(
                "Fuse group `$fuse_sym` (namespace `$ns`): view for `$(v.name)` is not a " *
                "SubArray (got $(typeof(view_data))) and its last-axis size " *
                "$(size(view_data, ndims(view_data))) does not match slot length $(length(slot))."
            )
        end
        for k in slot

            # chck the slot is within the parent
            (1 <= k <= parent_size) || error(
                "Fuse group `$fuse_sym` (namespace `$ns`): slot $k for `$(v.name)` " *
                "is outside the parent buffer's axis 1:$parent_size."
            )

            # check it's not already claimed
            covered[k] && error(
                "Fuse group `$fuse_sym` (namespace `$ns`): slot $k is claimed by " *
                "multiple members (offending member: `$(v.name)`)."
            )
            covered[k] = true
        end
    end
    # check everything is covered
    all(covered) || error(
        "Fuse group `$fuse_sym` (namespace `$ns`): slots do not tile the parent's " *
        "fused axis 1:$parent_size (uncovered slots: $(findall(!, covered)))."
    )
    return nothing
end

# Narrow the eltype of a variable vector to the concrete dim type so allocate_fused dispatches.
typed_vars(vars, ::Type{D}) where {D} = AbstractVariable{D}[v for v in vars]

"""$(TYPEDSIGNATURES) Allocates all variables within a group given a tuple of variables
expected to be `<: AbstractVariable` definitions. Determines the namespaces,
allocates the arrays with zeros and collects them into NamedTuples. Members of a
fuse group are returned as views of a shared parent (allocated via `build_fuse_parents`),
and the parent itself is exposed under the fuse symbol in the same NamedTuple."""
function allocate(group, model, fuse_parents = Dict{Tuple{Symbol, Symbol}, NamedTuple}())
    length(group) == 0 && return NamedTuple()  # return empty NamedTuple if no variables to initialize
    namespaces = filter(k -> k != Symbol(), tuple(keys(group)...))

    # variables without namespace identified by empty symbol Symbol() go directly into the main NamedTuple
    # that way we have variables.prognostic.vorticity skipping the namespace between prognostic and vor
    nt1 = haskey(group, Symbol()) ? _allocate_namespace(group[Symbol()], model, fuse_parents) : NamedTuple()

    # other variables grouped by namespace
    # e.g. variables.prognostic.ocean.sea_surface_temperature, variables.prognostic.land.soil_moisture, etc.
    nt2 = NamedTuple{namespaces}(
        Tuple(map(ns -> _allocate_namespace(group[ns], model, fuse_parents), namespaces))
    )

    return merge(nt1, nt2)
end

# Build the NamedTuple for a single namespace: fused members become views into a shared parent;
# standalone members allocate normally via `zero(v, model)`. The parent itself is NOT installed
# here — it lives at the canonical location `vars.scratch.fused.<fuse_symbol>` (built once at the
# end of `Variables(model)`).
function _allocate_namespace(vars, model, fuse_parents)
    pairs = Pair{Symbol, Any}[]
    for v in vars
        if v.fuse === Symbol()
            push!(pairs, v.name => zero(v, model))
        else
            entry = fuse_parents[(v.namespace, v.fuse)]
            push!(pairs, v.name => entry.views[(v.namespace, v.name)])
        end
    end
    return (; pairs...)
end

"""$(TYPEDSIGNATURES)
When model components are named tuples themselves then check for
variables required by the elements of the named tuple and pass those one as key-value pairs."""
variables(nt::NamedTuple, model::AbstractModel) = (variables(pair, model) for pair in pairs(nt)) |> Iterators.flatten |> Tuple
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
variables(component, model::AbstractModel) = variables(component)

"""$(TYPEDSIGNATURES)
Components can define `variables(pair::Pair{<:Symbol, <:AbstractComponent}, ::Model)` when they are part of a NamedTuple of components,
e.g. greenhouse gases. The key of the `pair` can then be used to define variables based on the name of the component.
The fallback defined here just drops the key so that it is optional."""
variables(name_component::Pair, model::AbstractModel) = variables(name_component.second, model)

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
