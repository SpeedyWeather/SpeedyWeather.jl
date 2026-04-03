abstract type AbstractLand <: AbstractModelComponent end
abstract type AbstractLandComponent <: AbstractModelComponent end
abstract type AbstractWetLand <: AbstractLand end
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

# LandModel defined through its components
export LandModel
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

# also allow spectral grid to be passed on as first an only positional argument to model constructors and nlayers to be passed as keyword argument directly to the constructor
(L::Type{<:AbstractLand})(SG::SpectralGrid; nlayers = DEFAULT_NLAYERS_SOIL, kwargs...) = L(spectral_grid = SG, geometry = LandGeometry(SG, nlayers = nlayers); kwargs...)

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
        # TODO think about the order of these
        timestep!(vars, land.snow, model)
        timestep!(vars, land.soil_moisture, model)
        timestep!(vars, land.vegetation, model)
        timestep!(vars, land.rivers, model)
    end
    timestep!(vars, land.temperature, model)
    return nothing
end

function initialize!(vars::Variables, land_model::AbstractLand, model)
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
