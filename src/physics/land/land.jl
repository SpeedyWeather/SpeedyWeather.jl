abstract type AbstractLand <: AbstractModelComponent end

include("temperature.jl")
include("soil_moisture.jl")
include("vegetation.jl")
include("rivers.jl")

# LandModel defined through its components
export LandModel
@kwdef struct LandModel{Temperature, SoilMoisture, Vegetation, Rivers} <: AbstractLand
    temperature::Temperature
    soil_moisture::SoilMoisture
    vegetation::Vegetation
    rivers::Rivers
end

# default constructor
function LandModel(SG::SpectralGrid)
    return LandModel(
        SeasonalLandTemperature(SG),
        SeasonalSoilMoisture(SG),
        VegetationClimatology(SG),
        NoRivers(SG),
    )
end

# initializing the land model initializes its components
function initialize!(   land::LandModel,
                        model::PrimitiveEquation)
    initialize!(model.land.temperature, model)
    initialize!(model.land.soil_moisture, model)
    initialize!(model.land.vegetation, model)
    initialize!(model.land.rivers, model)
end

export DryLandModel
@kwdef struct DryLandModel{Temperature} <: AbstractLand
    temperature::Temperature
end

DryLandModel(SG::SpectralGrid) = DryLandModel(SeasonalLandTemperature(SG))
initialize!(land::DryLandModel, model::PrimitiveEquation) = initialize!(land.temperature, model)

land_timestep!(progn::PrognosticVariables, diagn::DiagnosticVariables, model::PrimitiveEquation) = land_timestep!(progn, diagn, model.land, model)

function land_timestep!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    land::LandModel,
    model::PrimitiveWet,
)
    # TODO think about the order of these
    timestep!(progn, diagn, land.soil_moisture, model)
    timestep!(progn, diagn, land.rivers, model)
    timestep!(progn, diagn, land.vegetation, model)
    timestep!(progn, diagn, land.temperature, model)
end

function land_timestep!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    land::DryLandModel,
    model::PrimitiveDry,
)
    timestep!(progn, diagn, land.temperature, model)
end

function initialize!(
    land::PrognosticVariablesLand,  # for dispatch
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::PrimitiveWet,
)
    initialize!(progn, diagn, model.land.temperature, model)
    initialize!(progn, diagn, model.land.soil_moisture, model)
    initialize!(progn, diagn, model.land.vegetation, model)
    initialize!(progn, diagn, model.land.rivers, model)
end

function initialize!(
    land::PrognosticVariablesLand,  # for dispatch
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::PrimitiveDry,
)
    initialize!(progn, diagn, model.land.temperature, model)
end