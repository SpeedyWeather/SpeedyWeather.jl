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
struct SurfaceThermodynamicsConstant <: AbstractSurfaceThermodynamics end
SurfaceThermodynamicsConstant(SG::SpectralGrid) = SurfaceThermodynamicsConstant()
initialize!(::SurfaceThermodynamicsConstant,::PrimitiveEquation) = nothing

function surface_thermodynamics!(   column::ColumnVariables,
                                    ::SurfaceThermodynamicsConstant,
                                    model::PrimitiveWet)

    # surface value is same as lowest model level, use previous time step
    # for numerical stability
    column.surface_temp = column.temp[end]      # todo use constant POTENTIAL temperature
    column.surface_humid = column.humid[end]    # humidity at surface is the same as 

    # surface air density via virtual temperature
    (; R_dry) = model.atmosphere
    Tᵥ = column.temp_virt[column.nlayers]
    column.surface_air_density = column.pres[end]/(R_dry*Tᵥ)
end

function surface_thermodynamics!(   column::ColumnVariables,
                                    ::SurfaceThermodynamicsConstant,
                                    model::PrimitiveDry)
    (; R_dry) = model.atmosphere
    # surface value is same as lowest model level, but use previous
    # time step for numerical stability
    column.surface_temp = column.temp[end]   # todo use constant POTENTIAL temperature
    column.surface_air_density = column.pres[end]/(R_dry*column.surface_temp)
end