"""$TYPEDEF
Abstract type for all column-based parmaetrizations. Custom parametrizations are expected to 
subtype this and implement the [`variables`](@ref), [`initialize!`](@ref), and [`parameterization!`](@ref) for it. In
order to use the parameterization in a model, add it to the `parameterizations` of the model 
at definition."""
abstract type AbstractParameterization <: AbstractModelComponent end

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
