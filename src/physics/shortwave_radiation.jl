abstract type AbstractRadiation <: AbstractParameterization end
abstract type AbstractShortwave <: AbstractRadiation end

# function barrier for all AbstractShortwave
function shortwave_radiation!(column::ColumnVariables, model::PrimitiveEquation)
    shortwave_radiation!(column, model.shortwave_radiation, model)
end

## NO SHORTWAVE RADIATION
export NoShortwave
struct NoShortwave <: AbstractShortwave end
NoShortwave(SG::SpectralGrid) = NoShortwave()
initialize!(::NoShortwave, ::PrimitiveEquation) = nothing
shortwave_radiation!(::ColumnVariables, ::NoShortwave, ::PrimitiveEquation) = nothing

## SHORTWAVE RADIATION FOR A FULLY TRANSPARENT ATMOSPHERE
export TransparentShortwave
Base.@kwdef struct TransparentShortwave{NF} <: AbstractShortwave
    albedo::NF = 0.8
end

TransparentShortwave(SG::SpectralGrid; kwargs...) = TransparentShortwave{SG.NF}(; kwargs...)
initialize!(::TransparentShortwave, ::PrimitiveEquation) = nothing

function shortwave_radiation!(
    column::ColumnVariables,
    scheme::TransparentShortwave,
    model::PrimitiveEquation,
)
    shortwave_radiation!(column, scheme, model.planet)
end

function shortwave_radiation!(
    column::ColumnVariables,
    scheme::TransparentShortwave,
    planet::AbstractPlanet,
)
    (; cos_zenith) = column
    (; albedo) = scheme
    (; solar_constant) = planet
    column.flux_temp_upward[end] += (1-albedo) * solar_constant * cos_zenith
end