abstract type AbstractLand <: AbstractModelComponent end
abstract type AbstractWetLand <: AbstractLand end
abstract type AbstractDryLand <: AbstractLand end

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
end

include("geometry.jl")
include("thermodynamics.jl")
include("temperature.jl")
include("soil_moisture.jl")
include("vegetation.jl")
include("rivers.jl")

# LandModel defined through its components
export LandModel
@kwdef mutable struct LandModel{G, TD, T, SM, V, R} <: AbstractWetLand
    spectral_grid::SpectralGrid
    geometry::G = LandGeometry(spectral_grid)
    thermodynamics::TD = LandThermodynamics(spectral_grid)
    temperature::T = LandBucketTemperature(spectral_grid)
    soil_moisture::SM = LandBucketMoisture(spectral_grid)
    vegetation::V = VegetationClimatology(spectral_grid)
    rivers::R = nothing
end

# also allow spectral grid to be passed on as first an only positional argument to model constructors
(L::Type{<:AbstractLand})(SG::SpectralGrid; kwargs...) = L(spectral_grid=SG; kwargs...)

# initializing the land model initializes its components
function initialize!(   land::LandModel,
                        model::PrimitiveEquation)
    initialize!(model.land.geometry, model)
    initialize!(model.land.thermodynamics, model)
    initialize!(model.land.temperature, model)
    initialize!(model.land.soil_moisture, model)
    initialize!(model.land.vegetation, model)
    initialize!(model.land.rivers, model)
end

export DryLandModel
@kwdef struct DryLandModel{G, TD, T} <: AbstractDryLand
    spectral_grid::SpectralGrid
    geometry::G = LandGeometry(spectral_grid)
    thermodynamics::TD = LandThermodynamics(spectral_grid)
    temperature::T = LandBucketTemperature(spectral_grid)
end

function initialize!(land::DryLandModel, model::PrimitiveEquation)
    initialize!(model.land.geometry, model)
    initialize!(model.land.thermodynamics, model)
    initialize!(model.land.temperature, model)
end

land_timestep!(progn::PrognosticVariables, diagn::DiagnosticVariables, model::PrimitiveEquation) =
    land_timestep!(progn, diagn, model.land, model)

function land_timestep!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    land::AbstractLand,
    model::PrimitiveEquation,
)
    if model isa PrimitiveWet && land isa AbstractWetLand
        # TODO think about the order of these
        timestep!(progn, diagn, land.soil_moisture, model)
        timestep!(progn, diagn, land.rivers, model)
        timestep!(progn, diagn, land.vegetation, model)
    end
    timestep!(progn, diagn, land.temperature, model)
end

function initialize!(
    land::PrognosticVariablesLand,  # for dispatch
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::PrimitiveEquation,
)
    initialize!(progn, diagn, model.land.temperature, model)

    # only initialize soil moisture, vegetation, rivers if atmosphere and land are wet
    if model isa PrimitiveWet && model.land isa AbstractWetLand
        initialize!(progn, diagn, model.land.soil_moisture, model)      
        initialize!(progn, diagn, model.land.vegetation, model)     
        initialize!(progn, diagn, model.land.rivers, model)
    end
end