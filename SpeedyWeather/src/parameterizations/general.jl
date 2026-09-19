"""Function that defines the actual parameterization of an `AbstractParameterization`.

Takes in the index of the column `ij`, the `Variables` object `vars`, the
parameterization itself and the model. The model includes - among others - land sea mask,
orography and physical constants.

This function is used within a KernelAbstractions kernel and is therefore expected to work on GPU as well.
Don't use any dynamic dispatches, try to avoid allocations and branches in your code and only use scalar
indexing of arrays."""
parameterization!

"""$(TYPEDSIGNATURES) Fallback when setting `parameterization=nothing` in the model constructor."""
parameterization!(ij, vars::AbstractVariables, parameterization::Nothing, model) = nothing

"""$(TYPEDSIGNATURES) Fallback for parameterizations that don't have a global kernel."""
parameterization!(vars::AbstractVariables, parameterization::Any, model) = nothing

"""$(TYPEDSIGNATURES)
Initialize an `AbstractParameterization`. This is called once at when calling initialize!(model). 
The default behaviour is to return `nothing`."""
initialize!(parameterization::AbstractParameterization, model::AbstractModel) = nothing

"""
    Parameterization(spectral_grid::SpectralGrid, scheme; kwargs...)

Create a SpeedyWeather parameterization from `scheme`, a parameterization defined
in another package. That package does not have to depend on SpeedyWeather, instead
it adds a method to `Parameterization` in its SpeedyWeather extension, e.g.

```julia
# in the main code of ExternalPackage, no SpeedyWeather dependency
struct ExternalLongwave
    ...
end

# in ExternalPackageSpeedyWeatherExt
struct SpeedyExternalLongwave{NF} <: SpeedyWeather.AbstractLongwave
    ...
end

SpeedyWeather.Parameterization(spectral_grid::SpectralGrid, scheme::ExternalLongwave; kwargs...) =
    SpeedyExternalLongwave(spectral_grid, scheme; kwargs...)
```

which is then used as

```julia
using SpeedyWeather, ExternalPackage
spectral_grid = SpectralGrid()
longwave_radiation = Parameterization(spectral_grid, ExternalLongwave())
model = PrimitiveWetModel(spectral_grid; longwave_radiation)
```

Types defined in extensions cannot be exported, so this avoids having to access
them via `Base.get_extension`. Parameterizations that are already
an `AbstractParameterization` are returned unchanged."""
function Parameterization end

Parameterization(::SpectralGrid, parameterization::AbstractParameterization) = parameterization

"""$(TYPEDSIGNATURES)
Extract the parameterizations from the model as NamedTuple.
These are the GPU-compatible components of the model."""
@generated function get_parameterizations(model::ModelType) where {ModelType <: PrimitiveEquation}
    # Extract parameterization symbols from the type
    params_type = fieldtype(ModelType, :params)
    param_names = params_type.parameters[1]  # Extract tuple from Val{tuple}

    # Generate literal field accesses for type stability
    return :(NamedTuple{$param_names}(tuple($([:(model.$name) for name in param_names]...))))
end

@inline get_parameterizations(model::Barotropic) = NamedTuple()
@inline get_parameterizations(model::ShallowWater) = NamedTuple()
