for Var in (
        :PrognosticVariable,
        :GridVariable,
        :TendencyVariable,
        :DynamicsVariable,
        :ParameterizationVariable,
        :ParticleVariable,
        :ScratchVariable,
    )

    @eval begin
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

"""$(TYPEDSIGNATURES)
Symbol representing a unique identifier for a variable, combining its type, namespace and name and namespace. 
Used to remove duplicates when extracting variables from the model. E.g. `:Variable_ocean_sea_surface_temperature`."""
identifier(v::AbstractVariable) =
    Symbol(nonparametric_type(v), "_", v.namespace, v.namespace == Symbol() ? "" : "_", v.name)

Base.unique(variables::AbstractVariable...) = unique(v -> identifier(v), variables)

"""Struct that holds all variables of a simulation. Variables are split into 
groups corresponding in the fields here $(TYPEDFIELDS) Each group can have
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
            print(io, "\n$maybe_bar1 $s2 ")
            nt = getfield(getfield(V, p), k)
            if nt isa NamedTuple
                print(io, styled"{success:$k}")
                for (l, m) in enumerate(keys(nt))
                    s3 = l == length(keys(nt)) ? "└" : "├"  # choose ending └ for last variable
                    maybe_bar2 = lastj ? " " : "│"          # if last variable in namespace, no vertical bar needed
                    smry = Base.summary(getfield(nt, m))
                    print(io, "\n$maybe_bar1 $maybe_bar2 $s3 ", styled"{magenta:$m}: $smry")
                end
            else
                print(io, styled"{magenta:$k}")
                print(io, ": ", Base.summary(nt))
            end
        end
    end
    return nothing
end

# runic: off
function Variables(model::AbstractModel)
    all_vars = all_variables(model)     # one long tuple for all required variables of model and its components
    prognostic        = initialize_variables(filter_variables(all_vars,       PrognosticVariable), model)
    grid              = initialize_variables(filter_variables(all_vars,             GridVariable), model)
    tendencies        = initialize_variables(filter_variables(all_vars,         TendencyVariable), model)
    dynamics          = initialize_variables(filter_variables(all_vars,         DynamicsVariable), model)
    parameterizations = initialize_variables(filter_variables(all_vars, ParameterizationVariable), model)
    particles         = initialize_variables(filter_variables(all_vars,         ParticleVariable), model)
    scratch           = initialize_variables(filter_variables(all_vars,          ScratchVariable), model)
    return Variables(; prognostic, grid, tendencies, dynamics, parameterizations, particles, scratch)
end
# runic: on

function initialize_variables(variables, model)
    length(variables) == 0 && return NamedTuple()  # return empty NamedTuple if no variables to initialize
    namespaces = filter(k -> k != Symbol(), keys(variables))
    return merge(
        # variables without namespace identified by empty symbol Symbol() go directly into the main NamedTuple
        # that way we have variables.prognostic.vor skipping the namespace between prognostic and vor
        NamedTuple{Tuple(map(v -> v.name, variables[Symbol()]))}(Tuple(map(var -> zero(var, model), variables[Symbol()]))),

        # other variables grouped by namespace
        # e.g. variables.prognostic.ocean.sea_surface_temperature, variables.prognostic.land.soil_moisture, etc.
        NamedTuple{namespaces}(
            Tuple(
                map(
                    ns ->
                    NamedTuple{Tuple(map(v -> v.name, variables[ns]))}(Tuple(map(var -> zero(var, model), variables[ns])))
                    , namespaces
                )
            )
        )
    )
end

variables(::Nothing) = ()                                   # to allow for model.component = nothing
variables(::Any) = ()                                       # fallback for any component
variables(model::AbstractModel) = variables(typeof(model))

"""$(TYPEDSIGNATURES)
Define variables needed by a model component. This is used to extract all variables from the model and to initialize them.
The default is to return an empty tuple, but components can define their variables by extending this function for their type."""
variables

"""$(TYPEDSIGNATURES)
Extracts all variables from the model by iterating over all components and collecting their variables.
Return a tuple of all variables."""
function all_variables(model::AbstractModel)
    t = variables(model)                        # variables from the model itself
    for component in propertynames(model)       # iterate over all components of the model
        vars = variables(getproperty(model, component))
        if length(vars) > 0
            t = tuple(t..., vars...)
        end
    end
    return t
end

get_namespaces(variables::AbstractVariable...) = unique([v.namespace for v in variables])

function filter_variables(variables, VariableType)
    # filter by variable type, remove duplicates
    vars = unique(filter(v -> v isa VariableType, variables))

    # Split by namespace
    namespaces = get_namespaces(vars...)
    namespace_dict = Dict{Symbol, Tuple}()
    for ns in namespaces
        namespace_dict[ns] = Tuple([v for v in vars if v.namespace == ns])
    end

    # Convert to NamedTuple with tuples instead of vectors
    return NamedTuple{Tuple(namespaces)}(Tuple(namespace_dict[ns] for ns in namespaces))
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
