abstract type AbstractRadiation <: AbstractParameterization end
abstract type AbstractShortwave <: AbstractRadiation end

export NoShortwave
struct NoShortwave <: AbstractShortwave end
NoShortwave(SG::SpectralGrid) = NoShortwave()
initialize!(::NoShortwave, ::PrimitiveEquation) = nothing

# function barrier for all AbstractShortwave
function shortwave_radiation!(column::ColumnVariables, model::PrimitiveEquation)
    shortwave_radiation!(column, model.shortwave_radiation, model)
end

shortwave_radiation!(::ColumnVariables, ::NoShortwave, ::PrimitiveEquation) = nothing

export TransparentShortwave
Base.@kwdef struct TransparentShortwave{NF} <: AbstractShortwave
    albedo::NF = 0.3
    S::Base.RefValue{NF} = Ref(zero(NF))
end

TransparentShortwave(SG::SpectralGrid) = TransparentShortwave{SG.NF}()

function initialize!(scheme::TransparentShortwave, model::PrimitiveEquation)
    (;solar_constant, gravity) = model.planet
    (;pres_ref) = model.atmosphere
    cₚ = model.atmosphere.heat_capacity
    Δσ = model.geometry.σ_levels_thick
    scheme.S[] = (1 - scheme.albedo) * solar_constant * gravity / (Δσ[end] * pres_ref * cₚ)
end

function shortwave_radiation!(
    column::ColumnVariables,
    scheme::TransparentShortwave,
    model::PrimitiveEquation,
)
    shortwave_radiation!(column, scheme, model.solar_zenith)
end

function shortwave_radiation!(
    column::ColumnVariables,
    scheme::TransparentShortwave,
    solar_zenith::AbstractZenith,
)
    cos_zenith = solar_zenith.cos_zenith[column.ij]
    column.temp_tend[end] += cos_zenith*scheme.S[]
end