"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct SeaSurfaceTemperatureOutput{F} <: AbstractOutputVariable
    name::String = "sst"
    unit::String = "degC"
    long_name::String = "sea surface temperature"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 10
    transform::F = (x) -> x - 273.15
end

path(::SeaSurfaceTemperatureOutput, simulation) =
    simulation.prognostic_variables.ocean.sea_surface_temperature

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct SeaIceConcentrationOutput <: AbstractOutputVariable
    name::String = "sic"
    unit::String = "m²/m²"
    long_name::String = "sea ice concentration"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 10
end

path(::SeaIceConcentrationOutput, simulation) =
    simulation.prognostic_variables.ocean.sea_ice_concentration

OceanOutput() = (
    SeaSurfaceTemperatureOutput(),
    SeaIceConcentrationOutput(),
)