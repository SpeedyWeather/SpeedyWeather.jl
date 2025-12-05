export Grid2D, Grid3D, Spectral2D, Spectral3D, PrognosticVariable, DiagnosticVariable

"""
    $TYPEDEF

Indicator type for 2D variables on the spatial grid (x, y).
"""
struct Grid2D <: AbstractVariableDims end

"""
    $TYPEDEF

Indicator type for 3D variables on the spatial grid (x, y, z).
"""
struct Grid3D <: AbstractVariableDims end

"""
    $TYPEDEF

Indicator type for spectral 2D variables (l, m).
"""
struct Spectral2D <: AbstractVariableDims end

"""
    $TYPEDEF

Indicator type for spectral 3D variables (l, m, k).
"""
struct Spectral3D <: AbstractVariableDims end

# intialize variables
Base.zero(::AbstractVariable{Grid2D}, SG::SpectralGrid, nlayers::Int) = zeros(SG.GridVariable2D, SG.grid)
Base.zero(::AbstractVariable{Grid3D}, SG::SpectralGrid, nlayers::Int) = zeros(SG.GridVariable3D, SG.grid, nlayers)
Base.zero(::AbstractVariable{Spectral2D}, SG::SpectralGrid, nlayers::Int) = zeros(SG.SpectralVariable2D, SG.spectrum)
Base.zero(::AbstractVariable{Spectral3D}, SG::SpectralGrid, nlayers::Int) = zeros(SG.SpectralVariable3D, SG.spectrum, nlayers)

"""
    $TYPEDEF

Represents a prognostic variable with the given `name` and `dims`.
"""
@kwdef struct PrognosticVariable{
    VD, # <:VarDims
    UT, # <:Units
} <: AbstractVariable{VD}
    "Name of the prognostic variable"
    name::Symbol

    "Grid dimensions on which the variable is defined"
    dims::VD

    "Physical untis associated with this state variable"
    units::UT = nothing

    "Human-readable description of this state variable"
    desc::String = ""

    "Namespace of the variable, :land, :ocean or :atmosphere (default)"
    namespace::Symbol = :atmosphere
end

"""
    $TYPEDEF

Represents a diagnostic variable with the given `name` and `dims`.
"""
@kwdef struct DiagnosticVariable{
    VD, # <:VarDims
    UT, # <:Units
} <: AbstractVariable{VD}
    "Name of the diagnostic variable"
    name::Symbol

    "Grid dimensions on which the variable is defined"
    dims::VD

    "Physical untis associated with this state variable"
    units::UT = nothing

    "Human-readable description of this state variable"
    desc::String = ""

    "Namespace of the variable, :land, :ocean or :atmosphere (default)"
    namespace::Symbol = :atmosphere
end

"""
    remove_duplicate_variables(variables::AbstractVariable...) 

Remove duplicate variables based on their `name` and `namespace`.
Only the first occurrence of each unique combination is kept.
"""
function remove_duplicate_variables(variables::AbstractVariable...) 
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

"""
    $TYPEDSIGNATURES

Takes a tuple of `AbstractVariable`s, filters only the `DiagnosticVariable`s,
and splits them into a `NamedTuple` organized by their `namespace` field.

Returns a `NamedTuple` with keys corresponding to the namespaces found
(e.g., `:atmosphere`, `:land`, `:ocean`), where each value is a tuple
of `DiagnosticVariable`s belonging to that namespace.

# Example
```julia
vars = (
    DiagnosticVariables(name=:temp, dims=Grid3D(), namespace=:atmosphere),
    PrognosticVariable(name=:u, dims=Grid3D(), namespace=:atmosphere),
    DiagnosticVariables(name=:soil_temp, dims=Grid3D(), namespace=:land),
    DiagnosticVariables(name=:pressure, dims=Grid3D(), namespace=:atmosphere),
)

result = get_diagnostic_variables(vars)
# Returns: (atmosphere = (...), land = (...))
```
"""
function get_diagnostic_variables(variables::AbstractVariable...)
    # Filter only DiagnosticVariables
    diag_vars = filter(v -> v isa DiagnosticVariable, variables)
    
    # Remove duplicates
    unique_vars = remove_duplicate_variables(diag_vars...)
    
    # Currently hardcoded (might change in the future)
    namespaces = (:atmosphere, :land, :ocean)
 
    # Split by namespace
    namespace_dict = Dict{Symbol, Vector{DiagnosticVariable}}()
    for ns in namespaces
        namespace_dict[ns] = [v for v in unique_vars if v.namespace == ns]
    end
    
    # Convert to NamedTuple with tuples instead of vectors
    return NamedTuple{Tuple(namespaces)}(Tuple(Tuple(namespace_dict[ns]) for ns in namespaces))
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
    PrognosticVariable(name=:u, dims=Grid3D(), namespace=:atmosphere),
    DiagnosticVariables(name=:temp, dims=Grid3D(), namespace=:atmosphere),
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
    unique_vars = remove_duplicate_variables(prog_vars...)
    
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
get_diagnostic_variables(model::AbstractModel) = get_diagnostic_variables(variables(model)...)

# TODO: not quite sure yet about the nlayers or where it'll go
# initialize a NamedTuple from variables 
function initialize_variables(SG::SpectralGrid, nlayers::Integer, variables...) 
    return NamedTuple{Tuple(map(v -> v.name, variables))}(Tuple(map(var -> zero(var, SG, nlayers), variables)))
end 