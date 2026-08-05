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

const DEFAULT_NLAYERS_SOIL = 2

# model class is the abstract supertype
model_class(::Type{<:AbstractWetLand}) = AbstractWetLand
model_class(::Type{<:AbstractDryLand}) = AbstractDryLand
model_class(model::AbstractLand) = model_class(typeof(model))

# model type is the parameter-free type of a model
model_type(::Type{<:AbstractWetLand}) = LandModel
model_type(::Type{<:AbstractDryLand}) = DryLandModel
model_type(model::AbstractLand) = model_type(typeof(model))

@inline get_nlayers(model::AbstractLand) = model.geometry.nlayers
@inline get_nlayers(::Nothing) = 0 # fallback

Base.show(io::IO, M::AbstractLand) = Base.show(io, M, values = false)

export LandModel

"""Full land model with all components: geometry, thermodynamics, soil temperature (bucket scheme),
soil moisture, snow, vegetation, and optional rivers. Intended for use with `PrimitiveWet`. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef mutable struct LandModel{G, TD, T, SM, SN, V, R} <: AbstractWetLand
    spectral_grid::SpectralGrid
    @component geometry::G = LandGeometry(spectral_grid, nlayers = DEFAULT_NLAYERS_SOIL)
    @component thermodynamics::TD = LandThermodynamics(spectral_grid, geometry)
    @component temperature::T = LandBucketTemperature(spectral_grid, geometry)
    @component soil_moisture::SM = LandBucketMoisture(spectral_grid, geometry)
    @component snow::SN = SnowModel(spectral_grid, geometry)
    @component vegetation::V = VegetationClimatology(spectral_grid, geometry)
    @component rivers::R = nothing
end

# also allow spectral grid to be passed on as first an only positional argument to model constructors
# and nlayers to be passed as keyword argument directly to the constructor
(L::Type{<:AbstractLand})(SG::SpectralGrid; nlayers = DEFAULT_NLAYERS_SOIL, kwargs...) =
    L(spectral_grid = SG, geometry = LandGeometry(SG, nlayers = nlayers); kwargs...)

# initializing the land model initializes its components
function initialize!(land::LandModel, model::AbstractModel)
    initialize!(model.land.geometry, model)
    initialize!(model.land.thermodynamics, model)
    initialize!(model.land.temperature, model)
    initialize!(model.land.soil_moisture, model)
    initialize!(model.land.snow, model)
    initialize!(model.land.vegetation, model)
    initialize!(model.land.rivers, model)
    return nothing
end

# allocate variables as defined by land components
variables(land::LandModel) = (
    variables(land.temperature)...,
    variables(land.soil_moisture)...,
    variables(land.snow)...,
    variables(land.vegetation)...,
    variables(land.rivers)...,
)

export DryLandModel

"""Simplified land model with only geometry, thermodynamics, and soil temperature, omitting
hydrological components, i.e. no soil moisture, snow, vegetation, nor rivers. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef mutable struct DryLandModel{G, TD, T} <: AbstractDryLand
    spectral_grid::SpectralGrid
    @component geometry::G = LandGeometry(spectral_grid, nlayers = DEFAULT_NLAYERS_SOIL)
    @component thermodynamics::TD = LandThermodynamics(spectral_grid, geometry)
    @component temperature::T = LandBucketTemperature(spectral_grid, geometry)
end

function initialize!(land::DryLandModel, model::AbstractModel)
    initialize!(model.land.geometry, model)
    initialize!(model.land.thermodynamics, model)
    initialize!(model.land.temperature, model)
    return nothing
end

# initializing the land model initializes its components
variables(land::DryLandModel) = (variables(land.temperature)...,)

# unpack land model and call general timestep! function
land_timestep!(vars::Variables, model::PrimitiveEquation) =
    timestep!(vars, model.land, model)

function timestep!(
        vars::Variables,
        land::AbstractLand,
        model::PrimitiveEquation,
    )
    if model isa PrimitiveWet && land isa AbstractWetLand
        timestep!(vars, land.snow, model)
        timestep!(vars, land.soil_moisture, model)
        timestep!(vars, land.vegetation, model)
        timestep!(vars, land.rivers, model)
    end
    timestep!(vars, land.temperature, model)
    return nothing
end

function initialize!(vars::Variables, land_model::AbstractLand, model::PrimitiveEquation)
    # only initialize soil moisture, vegetation, rivers if atmosphere and land are wet
    if model isa PrimitiveWet && land_model isa AbstractWetLand
        initialize!(vars, land_model.snow, model)
        initialize!(vars, land_model.soil_moisture, model)
        initialize!(vars, land_model.vegetation, model)
        initialize!(vars, land_model.rivers, model)
    end
    initialize!(vars, land_model.temperature, model)
    return nothing
end

# any land model has a nothing fallback for filter! (change variables after time step)
Base.filter!(vars::Variables, land_model::AbstractLand, model::PrimitiveEquation) = nothing

# LandModel passes on filter! to its components
function Base.filter!(vars::Variables, land_model::LandModel, model::PrimitiveEquation)
    # only initialize soil moisture, vegetation, rivers if atmosphere and land are wet
    if model isa PrimitiveWet && land_model isa AbstractWetLand
        filter!(vars, land_model.snow, model)
        filter!(vars, land_model.soil_moisture, model)
        filter!(vars, land_model.vegetation, model)
        filter!(vars, land_model.rivers, model)
    end
    filter!(vars, land_model.temperature, model)
    return nothing
end