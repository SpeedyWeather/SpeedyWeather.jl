"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct SensibleHeatFluxOutput <: AbstractOutputVariable
    name::String = "shf"
    unit::String = "W/m^2"
    long_name::String = "Sensible heat fluxes (positive up)"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::SensibleHeatFluxOutput, simulation) =
    simulation.diagnostic_variables.physics.sensible_heat_flux

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct EvaporativeFluxOutput <: AbstractOutputVariable
    name::String = "ef"
    unit::String = "kg/s/m^2"
    long_name::String = "Surface humidity fluxes (positive up)"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::EvaporativeFluxOutput, simulation) =
    simulation.diagnostic_variables.physics.evaporative_flux

SurfaceFluxesOutput() = (
    SensibleHeatFluxOutput(),
    EvaporativeFluxOutput(),
)