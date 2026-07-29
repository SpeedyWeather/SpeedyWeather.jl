"""$TYPEDEF
Abstract type for all column-based parmaetrizations. Custom parametrizations are expected to 
subtype this and implement the [`variables`](@ref), [`initialize!`](@ref), and [`parameterization!`](@ref) for it. In
order to use the parameterization in a model, add it to the `parameterizations` of the model 
at definition."""
abstract type AbstractParameterization <: AbstractModelComponent end
abstract type AbstractOcean <: AbstractModelComponent end
abstract type AbstractSeaIce <: AbstractModelComponent end
