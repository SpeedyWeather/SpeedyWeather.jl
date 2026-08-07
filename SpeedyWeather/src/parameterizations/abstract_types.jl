"""$TYPEDEF
Abstract type for all column-based parmaetrizations. Custom parametrizations are expected to 
subtype this and implement the [`variables`](@ref), [`initialize!`](@ref), and [`parameterization!`](@ref) for it. In
order to use the parameterization in a model, add it to the `parameterizations` of the model 
at definition."""
abstract type AbstractParameterization <: AbstractModelComponent end
abstract type AbstractOcean <: AbstractModelComponent end
abstract type AbstractSeaIce <: AbstractModelComponent end

"""Abstract supertype for all land models. Custom land models should subtype `AbstractWetLand`
or `AbstractDryLand` and implement `initialize!`, `timestep!`, and optionally `variables` and `filter!`."""
abstract type AbstractLand <: AbstractModelComponent end

"""Abstract supertype for sub-components of a land model (soil temperature, soil moisture,
snow, vegetation, rivers). Each component implements `initialize!`, `timestep!`, and optionally
`variables` and `filter!`."""
abstract type AbstractLandComponent <: AbstractModelComponent end

"""Abstract land model type that includes hydrological processes (soil moisture, snow,
vegetation, rivers). Intended for use with `PrimitiveWet`. The default concrete type is `LandModel`."""
abstract type AbstractWetLand <: AbstractLand end

"""Abstract land model type that omits hydrological processes and only includes soil temperature.
Can be used with both `PrimitiveDry` and `PrimitiveWet`. The default concrete type is `DryLandModel`."""
abstract type AbstractDryLand <: AbstractLand end