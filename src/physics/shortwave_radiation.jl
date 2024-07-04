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
struct TransparentShortwave <: AbstractShortwave end
TransparentShortwave(SG::SpectralGrid) = TransparentShortwave()
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
    (; cos_zenith, albedo) = column
    (; solar_constant) = planet

    # reduce strength by 1/4 as this currently just hits the air temperature in lowermost
    # layer directly and not through a skin temperature
    column.flux_temp_upward[end] += ((1-albedo) * solar_constant * cos_zenith)/4
end