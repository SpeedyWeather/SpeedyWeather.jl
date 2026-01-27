abstract type AbstractVariables end
abstract type AbstractVariableDims end
abstract type AbstractVariable{D <: AbstractVariableDims} end

@kwdef struct Variables{P, G, T, D, Ph, S} <: AbstractVariables
    prognostic::P = NamedTuple()
    grid::G = NamedTuple()
    tendencies::T = NamedTuple()
    dynamics::D = NamedTuple()
    physics::Ph = NamedTuple()
    scratch::S = NamedTuple()
end

# runic: off
function Variables(model::AbstractModel)
    progn      = initialize_variables(get_prognostic_variables(model), model)
    grid       = initialize_variables(      get_grid_variables(model), model)
    tendencies = initialize_variables(  get_tendency_variables(model), model)
    dynamics   = initialize_variables(  get_dynamics_variables(model), model)
    physics    = initialize_variables(   get_physics_variables(model), model)
    scratch    = initialize_variables(   get_scratch_variables(model), model)
    return Variables(progn, grid, tendencies, dynamics, physics, scratch)
end
# runic: on

struct Grid2D <: AbstractVariableDims end
struct Grid3D <: AbstractVariableDims end
struct Land3D <: AbstractVariableDims end       # TODO redundant given Grid3D exists?
struct Ocean3D <: AbstractVariableDims end
struct Spectral2D <: AbstractVariableDims end
struct Spectral3D <: AbstractVariableDims end
struct Vertical1D <: AbstractVariableDims end

struct Grid4D <: AbstractVariableDims end
struct Spectral4D <: AbstractVariableDims end

zeeros(::Type{T}, k...) where {T} = fill!(T(undef, k...), 0)

# runic: off
Base.zero(::AbstractVariable{Grid2D},     SG::SpectralGrid, nlayers::Int) = zeros(SG.GridVariable2D, SG.grid)
Base.zero(::AbstractVariable{Grid3D},     SG::SpectralGrid, nlayers::Int) = zeros(SG.GridVariable3D, SG.grid, nlayers)
Base.zero(::AbstractVariable{Land3D},     SG::SpectralGrid, nlayers::Int) = zeros(SG.GridVariable3D, SG.grid, nlayers)
Base.zero(::AbstractVariable{Ocean3D},    SG::SpectralGrid, nlayers::Int) = zeros(SG.GridVariable3D, SG.grid, nlayers)
Base.zero(::AbstractVariable{Spectral2D}, SG::SpectralGrid, nlayers::Int) = zeros(SG.SpectralVariable2D, SG.spectrum)
Base.zero(::AbstractVariable{Spectral3D}, SG::SpectralGrid, nlayers::Int) = zeros(SG.SpectralVariable3D, SG.spectrum, nlayers)
Base.zero(::AbstractVariable{Vertical1D}, SG::SpectralGrid, nlayers::Int) = zeeros(SG.VectorType, nlayers)

Base.zero(::AbstractVariable{Grid4D},     SG::SpectralGrid, k::Int...) = zeros(SG.GridVariable3D, SG.grid, k...)
Base.zero(::AbstractVariable{Spectral4D}, SG::SpectralGrid, k::Int...) = zeros(SG.GridVariable3D, SG.grid, k...)
# runic: on

for Var in (
        :PrognosticVariable,
        :GridVariable,
        :TendencyVariable,
        :DynamicsVariable,
        :PhysicsVariable,
        :ScratchVariable,
    )

    @eval begin
        @kwdef struct $Var{D, U} <: AbstractVariable{D}
            name::Symbol
            dims::D
            units::U = ""
            desc::String = ""
            namespace::Symbol = Symbol()
        end

        nonparametric_type(v::$Var) = nonparametric_type(typeof(v))
        nonparametric_type(::Type{<:$Var}) = $Var
        $Var(name, dims; kwargs...) = $Var(; name, dims, kwargs...)
        $Var(name, dims, units; kwargs...) = $Var(; name, dims, units, kwargs...)
        $Var(name, dims, units, desc; kwargs...) = $Var(; name, dims, units, desc, kwargs...)
        $Var(name, dims, units, desc, namespace) = $Var(; name, dims, units, desc, namespace)
    end
end

identifier(v::AbstractVariable) = Symbol(nonparametric_type(v), "_", v.namespace, v.namespace == Symbol() ? "" : "_", v.name)
get_namespaces(variables::Tuple) = get_namespaces(variables...)
get_namespaces(variables::AbstractVariable...) = unique([v.namespace for v in variables])

function named_tuple(::Type{V}, variables::AbstractVariable...) where {V <: AbstractVariable}
    all_vars = filter(v -> v isa V, variables)
    unique_vars = remove_duplicate_variables(all_vars...)
    namespaces = get_namespaces(unique_vars)

    namespace_dict = Dict{Symbol, Vector{V}}()
    for ns in namespaces
        namespace_dict[ns] = [v for v in unique_vars if v.namespace == ns]
    end

    return NamedTuple{Tuple(namespaces)}(Tuple(Tuple(namespace_dict[ns]) for ns in namespaces))
end

Base.unique(variables::AbstractVariable...) = unique(v -> identifier(v), variables)

remove_duplicates(variables::Tuple) = remove_duplicates(variables...)
function remove_duplicates(variables::AbstractVariable...)
    seen = Set{Tuple{Symbol, Symbol}}()
    unique_vars = []
    for v in variables
        key = (v.name, v.namespace)
        if !(key in seen)
            push!(seen, key)
            push!(unique_vars, v)
        end
    end
    return unique_vars
end

function initialize_variables(SG::SpectralGrid, nlayers::Integer, variables...)
    return NamedTuple{Tuple(map(v -> v.name, variables))}(Tuple(map(var -> zero(var, SG, nlayers), variables)))
end


"""
    $TYPEDSIGNATURES

Takes a tuple of `AbstractVariable`s, filters only the `PrognosticVariable`s,
and splits them into a `NamedTuple` organized by their `namespace` field.

Returns a `NamedTuple` with keys corresponding to the namespaces found
(e.g., `:atmosphere`, `:land`, `:ocean`), where each value is a tuple
of `PrognosticVariable`s belonging to that namespace.

# Example
```julia
vars = (
    PrognosticVariable(name=:u, dims=Grid3D(), namespace=:atmosphere),``
    DiagnosticVariable(name=:temp, dims=Grid3D(), namespace=:atmosphere),
    PrognosticVariable(name=:soil_moisture, dims=Grid3D(), namespace=:land),
    PrognosticVariable(name=:v, dims=Grid3D(), namespace=:atmosphere),
)

result = get_prognostic_variables(vars)
# Returns: (atmosphere = (...), land = (...))
```
"""
function get_prognostic_variables(variables::AbstractVariable...)
    # Filter only PrognosticVariable
    prog_vars = filter(v -> v isa PrognosticVariable, variables)

    # Remove duplicates
    unique_vars = remove_duplicates(prog_vars...)

    # Currently hardcoded (might change in the future)
    namespaces = (:atmosphere, :land, :ocean)

    # Split by namespace
    namespace_dict = Dict{Symbol, Vector{PrognosticVariable}}()
    for ns in namespaces
        namespace_dict[ns] = [v for v in unique_vars if v.namespace == ns]
    end

    # Convert to NamedTuple with tuples instead of vectors
    return NamedTuple{Tuple(namespaces)}(Tuple(Tuple(namespace_dict[ns]) for ns in namespaces))
end

get_prognostic_variables(model::AbstractModel) = get_prognostic_variables(variables(model)...)