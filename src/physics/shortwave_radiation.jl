abstract type AbstractRadiation <: AbstractParameterization end
abstract type AbstractShortwave <: AbstractRadiation end

function get_nbands(R::AbstractRadiation)
    hasfield(typeof(R), :nbands) && return R.nbands
    return 0
end

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

    # transparent = optical thickness of zero, no vertical changes in flux
    # this will sum up to zero in every layer (=transparent) but yields
    # a non-zero net flux at the surface 
    column.flux_temp_downward .+= solar_constant * cos_zenith
    column.flux_temp_upward .+= albedo * solar_constant * cos_zenith

    # diagnostic: outgoing shortwave radiation
    column.outgoing_shortwave_radiation = albedo * solar_constant * cos_zenith

    return nothing
end