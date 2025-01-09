"""Defines netCDF output for a specific variables, see `VorticityOutput` for details.
Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct OutgoingLongwaveRadiationOutput <: AbstractOutputVariable
    name::String = "olr"
    unit::String = "W/m^2"
    long_name::String = "Outgoing longwave radiation"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::OutgoingLongwaveRadiationOutput, simulation) =
    simulation.diagnostic_variables.physics.outgoing_longwave_radiation

"""Defines netCDF output for a specific variables, see `VorticityOutput` for details.
Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct OutgoingShortwaveRadiationOutput <: AbstractOutputVariable
    name::String = "osr"
    unit::String = "W/m^2"
    long_name::String = "Outgoing shortwave radiation"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::OutgoingShortwaveRadiationOutput, simulation) =
    simulation.diagnostic_variables.physics.outgoing_shortwave_radiation

RadiationOutput() = (OutgoingLongwaveRadiationOutput(), OutgoingShortwaveRadiationOutput())