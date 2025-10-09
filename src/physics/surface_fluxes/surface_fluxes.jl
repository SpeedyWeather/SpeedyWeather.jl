abstract type AbstractSurfaceThermodynamics <: AbstractParameterization end
abstract type AbstractSurfaceWind <: AbstractParameterization end
abstract type AbstractSurfaceHeatFlux <: AbstractParameterization end
abstract type AbstractSurfaceHumidityFlux <: AbstractParameterization end

# skip immediately for sensi_heat_flux or 
surface_heat_flux!(c::ColumnVariables, f::Nothing, ::PrognosticVariables, ::PrimitiveEquation) = nothing
surface_humidity_flux!(c::ColumnVariables, f::Nothing, ::PrognosticVariables, ::PrimitiveEquation) = nothing

# skip the prognostic variables if not defined (needed to read in prescribed fluxes)
surface_heat_flux!(c::ColumnVariables, f::AbstractSurfaceHeatFlux,
    p::PrognosticVariables, m::PrimitiveEquation) = surface_heat_flux!(c, f, m)

surface_humidity_flux!(c::ColumnVariables, f::AbstractSurfaceHumidityFlux,
    p::PrognosticVariables, m::PrimitiveEquation) = surface_humidity_flux!(c, f, m)

# defines the order in which they are called und unpacks to dispatch
function surface_fluxes!(
    column::ColumnVariables,
    progn::PrognosticVariables,
    model::PrimitiveEquation)

    # get temperature, humidity and density at surface
    surface_thermodynamics!(column, model.surface_thermodynamics, model)

    # also calculates surface wind speed necessary for other fluxes too
    surface_wind_stress!(column, model.surface_wind, model)

    # now call other heat (wet and dry) and humidity fluxes (PrimitiveWet only)
    surface_heat_flux!(column, model.surface_heat_flux, progn, model)
    model isa PrimitiveWet && surface_humidity_flux!(column, model.surface_humidity_flux, progn, model)
end

## SURFACE THERMODYNAMICS
export SurfaceThermodynamicsConstant
@kwdef mutable struct SurfaceThermodynamicsConstant{NF} <: AbstractSurfaceThermodynamics
    σ_power_minus_κ::NF = 0
end

SurfaceThermodynamicsConstant(SG::SpectralGrid; kwargs...) = SurfaceThermodynamicsConstant{SG.NF}(; kwargs...)
function initialize!(S::SurfaceThermodynamicsConstant, model::PrimitiveEquation)
    (; κ) = model.atmosphere
    σ = model.geometry.σ_levels_full[end]
    S.σ_power_minus_κ = σ^-κ
end

function surface_thermodynamics!(   column::ColumnVariables,
                                    S::SurfaceThermodynamicsConstant,
                                    model::PrimitiveWet)
    (; R_dry) = model.atmosphere
    σ⁻ᵏ = S.σ_power_minus_κ

    # constant potential temperature
    column.surface_temp = column.temp[end]*σ⁻ᵏ
    column.surface_humid = column.humid[end]    # humidity at surface is the same as in lowermost layer

    (; R_dry) = model.atmosphere
    column.surface_air_density = column.pres[end]/(R_dry*column.surface_temp)
end

function surface_thermodynamics!(   column::ColumnVariables,
                                    S::SurfaceThermodynamicsConstant,
                                    model::PrimitiveDry)
    (; R_dry) = model.atmosphere
    σ⁻ᵏ = S.σ_power_minus_κ

    # constant potential temperature
    column.surface_temp = column.temp[end]*σ⁻ᵏ
    column.surface_air_density = column.pres[end]/(R_dry*column.surface_temp)
end