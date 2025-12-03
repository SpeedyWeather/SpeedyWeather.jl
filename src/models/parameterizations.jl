"""
    $TYPEDEF

Abstract type for all column-based parmaetrizations. Custom parametrizations are expected to 
subtype this and implement the [`variables`](@ref), [`initialize!`](@ref), and [`parameterization!`](@ref) for it. In
order to use the parameterization in a model, add it to the `parameterizations` of the model 
at definition.
"""
abstract type AbstractParameterization <: AbstractModelComponent end

"""
    $(TYPEDSIGNATURES)

Function that defines the actual parameterization of an `AbstractParameterization`. 

Takes in the index of the column `ij`, the `DiagnosticVariables` and `PrognosticVariables`, the 
parameterization itself and a tuple with all relevant model parameters (all fields in `model.parameters`). 
The latter includes - among others - land sea mask, orography and physical constants. 

This function is used within a KernelAbstractions kernel and is therefore expected to work on GPU as well. 
Don't use any dynamic dispatches, try to avoid allocations and branches in your code and only use scalar 
indexing of arrays. 
"""
# parameterization!(ij, diagn, progn, parameterization::AbstractParameterization, model_parameters) = nothing
parameterization!(ij, diagn, progn, parameterization::Nothing, model_parameters) = nothing

"""
    $(TYPEDSIGNATURES)

Initialize an `AbstractParameterization`. This is called once at when calling initialize!(model). 
The default behaviour is to return `nothing`
"""
initialize!(parameterization::AbstractParameterization, model::AbstractModel) = nothing

"""
    $(TYPEDSIGNATURES)

Function that defines which variables a `AbstractParameterization` needs to use and allocate,
returns a tuple of [`AbstractVariable`](@ref).
"""
variables(parameterization::AbstractParameterization) = () # fallback for parameterizations that don't need variables

# variables of nothing are nothing
variables(::Nothing) = ()

function variables(model::AbstractModel)
    # Collect all variables from all parameterizations and flatten into a single tuple
    all_vars = Tuple(vcat([collect(SpeedyWeather.variables(component)) for component in SpeedyWeather.get_all_parameterizations(model)]...))
    return all_vars
end

"""$(TYPEDSIGNATURES)
Extract the parameterizations from the model including land and ocean, to infer variables."""
function get_all_parameterizations(model::AbstractModel)
    return merge(get_parameterizations(model), get_extra_parameterizations(model))
end
