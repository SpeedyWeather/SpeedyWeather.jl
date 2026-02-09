abstract type AbstractVariableDims end                          # Dimensions of variables
abstract type AbstractVariable{D <: AbstractVariableDims} end   # Variable type with dimensions D
abstract type AbstractVariables end                             # Container for all variable arrays in the model

# DIMENSIONS --------------------------
# define a a bunch of variable dimensions used to allocate variables of the right size and type
struct ClockVar <: AbstractVariableDims end
struct Grid2D <: AbstractVariableDims end
struct Grid3D <: AbstractVariableDims end
struct Land3D <: AbstractVariableDims end
struct Ocean3D <: AbstractVariableDims end
struct Spectral2D <: AbstractVariableDims end
struct Spectral3D <: AbstractVariableDims end
struct Latitude1D <: AbstractVariableDims end
struct Vertical1D <: AbstractVariableDims end
@kwdef struct Grid4D <: AbstractVariableDims
    n::Int = 1                                                  # length of 4th dimension
end
@kwdef struct Spectral4D <: AbstractVariableDims
    n::Int = 1
end
@kwdef struct Array1D <: AbstractVariableDims
    n::Int = 1
end
@kwdef struct Array2D <: AbstractVariableDims
    m::Int = 1; n::Int = 1
end

# zero constructor for arrays but matching the grid/spectrum syntax
zeeros(::Type{T}, k...) where {T} = fill!(T(undef, k...), 0)

# runic: off
Base.zero(::AbstractVariable{ClockVar},     model::AbstractModel) = Clock()
Base.zero(::AbstractVariable{Grid2D},       model::AbstractModel) = zeros( model.spectral_grid.GridVariable2D, model.spectral_grid.grid)
Base.zero(::AbstractVariable{Grid3D},       model::AbstractModel) = zeros( model.spectral_grid.GridVariable3D, model.spectral_grid.grid, get_nlayers(model))
Base.zero(::AbstractVariable{Land3D},       model::AbstractModel) = zeros( model.spectral_grid.GridVariable3D, model.spectral_grid.grid, get_nlayers(model.land))
Base.zero(::AbstractVariable{Ocean3D},      model::AbstractModel) = zeros( model.spectral_grid.GridVariable3D, model.spectral_grid.grid, get_nlayers(model.ocean))
Base.zero(::AbstractVariable{Spectral2D},   model::AbstractModel) = zeros( model.spectral_grid.SpectralVariable2D, model.spectral_grid.spectrum)
Base.zero(::AbstractVariable{Spectral3D},   model::AbstractModel) = zeros( model.spectral_grid.SpectralVariable3D, model.spectral_grid.spectrum, get_nlayers(model))
Base.zero(::AbstractVariable{Latitude1D},   model::AbstractModel) = zeeros(model.spectral_grid.VectorType, model.spectral_grid.nlat)
Base.zero(::AbstractVariable{Vertical1D},   model::AbstractModel) = zeeros(model.spectral_grid.VectorType, get_nlayers(model))
Base.zero(v::AbstractVariable{Array1D},     model::AbstractModel) = zeeros(model.spectral_grid.VectorType, v.dims.n)
Base.zero(v::AbstractVariable{Array2D},     model::AbstractModel) = zeeros(model.spectral_grid.MatrixType, v.dims.m, v.dims.n)
Base.zero(v::AbstractVariable{Grid4D},      model::AbstractModel) = zeros( model.spectral_grid.GridVariable3D, model.spectral_grid.grid, v.dims.n)
Base.zero(v::AbstractVariable{Spectral4D},  model::AbstractModel) = zeros( model.spectral_grid.SpectralVariable3D, model.spectral_grid.spectrum, v.dims.n)
# runic: on

# VARIABLE TYPES --------------------------
for Var in (
        :PrognosticVariable,
        :GridVariable,
        :TendencyVariable,
        :DynamicsVariable,
        :ParameterizationVariable,
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

@kwdef struct Variables{P, G, T, D, Pa, S} <: AbstractVariables
    prognostic::P = NamedTuple()
    grid::G = NamedTuple()
    tendencies::T = NamedTuple()
    dynamics::D = NamedTuple()
    parameterizations::Pa = NamedTuple()
    scratch::S = NamedTuple()
end

function Base.show(io::IO, V::Variables)
    print(io, "Variables{@NamedTuple{...}, ...}")
    for (i, p) in enumerate(propertynames(V))
        lasti = i == length(propertynames(V))   # check if last property to choose ending └
        s = lasti ? "└" : "├"                   # choose ending
        print(io, "\n$s $p")
        for (j, k) in enumerate(keys(getfield(V, p)))
            lastj = j == length(keys(getfield(V, p)))   # check if last variable in namespace to choose ending └
            s2 = lastj ? "└" : "├"                      # choose ending └ for last variable
            maybe_bar1 = lasti ? " " : "│"               # if last variable in namespace, no vertical bar needed
            print(io, "\n$maybe_bar1 $s2 $k")
            nt = getfield(getfield(V, p), k)
            if nt isa NamedTuple
                for (l, m) in enumerate(keys(nt))
                    s3 = l == length(keys(nt)) ? "└" : "├"  # choose ending └ for last variable
                    maybe_bar2 = lastj ? " " : "│"           # if last variable in namespace, no vertical bar needed
                    print(io, "\n$maybe_bar1 $maybe_bar2 $s3 $m: ", Base.summary(getfield(nt, m)))
                end
            else
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
    scratch           = initialize_variables(filter_variables(all_vars,          ScratchVariable), model)
    return Variables(; prognostic, grid, tendencies, dynamics, parameterizations, scratch)
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
