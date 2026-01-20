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

function Base.show(io::IO, M::AbstractLand)
    println(io, "$(model_type(M)) <: $(model_class(M))")
    properties = propertynames(M)
    n = length(properties)
    for (i, key) in enumerate(properties)
        val = getfield(M, key)
        s = i == n ? "└" : "├"  # choose ending └ for last property
        p = i == n ? print : println
        p(io, "$s $key: $(typeof(val))")
    end
    return
end

# LandModel defined through its components
export LandModel
@parameterized @kwdef mutable struct LandModel{G, TD, T, SM, SN, V, R} <: AbstractWetLand
    spectral_grid::SpectralGrid
    @component geometry::G = LandGeometry(spectral_grid)
    @component thermodynamics::TD = LandThermodynamics(spectral_grid)
    @component temperature::T = LandBucketTemperature(spectral_grid)
    @component soil_moisture::SM = LandBucketMoisture(spectral_grid)
    @component snow::SN = SnowModel(spectral_grid)
    @component vegetation::V = VegetationClimatology(spectral_grid)
    @component rivers::R = nothing
end

# also allow spectral grid to be passed on as first an only positional argument to model constructors
(L::Type{<:AbstractLand})(SG::SpectralGrid; kwargs...) = L(spectral_grid = SG; kwargs...)

# initializing the land model initializes its components
function initialize!(
        land::LandModel,
        model::PrimitiveEquation
    )
    initialize!(model.land.geometry, model)
    initialize!(model.land.thermodynamics, model)
    initialize!(model.land.temperature, model)
    initialize!(model.land.soil_moisture, model)
    initialize!(model.land.snow, model)
    initialize!(model.land.vegetation, model)
    return initialize!(model.land.rivers, model)
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
@parameterized @kwdef struct DryLandModel{G, TD, T} <: AbstractDryLand
    spectral_grid::SpectralGrid
    @component geometry::G = LandGeometry(spectral_grid)
    @component thermodynamics::TD = LandThermodynamics(spectral_grid)
    @component temperature::T = LandBucketTemperature(spectral_grid)
end

function initialize!(land::DryLandModel, model::PrimitiveEquation)
    initialize!(model.land.geometry, model)
    initialize!(model.land.thermodynamics, model)
    return initialize!(model.land.temperature, model)
end

# initializing the land model initializes its components
variables(land::DryLandModel) = (variables(land.temperature)...,)

# unpack land model and call general timestep! function
land_timestep!(progn::PrognosticVariables, diagn::DiagnosticVariables, model::PrimitiveEquation) =
    timestep!(progn, diagn, model.land, model)

function timestep!(
        progn::PrognosticVariables,
        diagn::DiagnosticVariables,
        land::AbstractLand,
        model::PrimitiveEquation,
    )
    if model isa PrimitiveWet && land isa AbstractWetLand
        # TODO think about the order of these
        timestep!(progn, diagn, land.snow, model)
        timestep!(progn, diagn, land.soil_moisture, model)
        timestep!(progn, diagn, land.rivers, model)
        timestep!(progn, diagn, land.vegetation, model)
    end
    return timestep!(progn, diagn, land.temperature, model)
end

function initialize!(
        land::PrognosticVariablesLand,  # for dispatch
        progn::PrognosticVariables,
        diagn::DiagnosticVariables,
        land_model::AbstractLand,
        model::PrimitiveEquation,
    ) where {PrognosticVariablesLand}
    initialize!(progn, diagn, land_model.temperature, model)

    # only initialize soil moisture, vegetation, rivers if atmosphere and land are wet
    return if model isa PrimitiveWet && land_model isa AbstractWetLand
        initialize!(progn, diagn, land_model.soil_moisture, model)
        initialize!(progn, diagn, land_model.snow, model)
        initialize!(progn, diagn, land_model.vegetation, model)
        initialize!(progn, diagn, land_model.rivers, model)
    end
end
