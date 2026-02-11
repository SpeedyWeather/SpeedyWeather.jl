"""
    $TYPEDEF

Abstract type for all column-based parmaetrizations. Custom parametrizations are expected to 
subtype this and implement the [`variables`](@ref), [`initialize!`](@ref), and [`parameterization!`](@ref) for it. In
order to use the parameterization in a model, add it to the `parameterizations` of the model 
at definition.
"""
abstract type AbstractParameterization <: AbstractModelComponent end

"""Function that defines the actual parameterization of an `AbstractParameterization`. 

Takes in the index of the column `ij`, the `DiagnosticVariables` and `PrognosticVariables`, the 
parameterization itself and a tuple with all relevant model parameters (all fields in `model.parameters`). 
The latter includes - among others - land sea mask, orography and physical constants. 

This function is used within a KernelAbstractions kernel and is therefore expected to work on GPU as well. 
Don't use any dynamic dispatches, try to avoid allocations and branches in your code and only use scalar 
indexing of arrays."""
parameterization!

"""$(TYPEDSIGNATURES) Fallback when setting `parameterization=nothing` in the model constructor."""
parameterization!(ij, diagn, progn, parameterization::Nothing, model_parameters) = nothing

"""
    $(TYPEDSIGNATURES)

Initialize an `AbstractParameterization`. This is called once at when calling initialize!(model). 
The default behaviour is to return `nothing`."""
initialize!(parameterization::AbstractParameterization, model::AbstractModel) = nothing

"""$(TYPEDSIGNATURES)
Extract the parameterizations from the model including land and ocean, to infer variables."""
function get_all_parameterizations(model::AbstractModel)
    return merge(get_parameterizations(model), get_extra_parameterizations(model))
end
