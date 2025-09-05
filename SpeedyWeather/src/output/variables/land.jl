"""Defines netCDF output for a specific variable, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct SoilTemperatureOutput{F} <: AbstractOutputVariable
    name::String = "st"
    unit::String = "degC"
    long_name::String = "soil temperature"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)
    is_land::Bool = true    # 3D land variables have another vertical dimension
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 10
    transform::F = (x) -> x - 273.15
end

path(::SoilTemperatureOutput, simulation) =
    simulation.prognostic_variables.land.soil_temperature

"""Defines netCDF output for a specific variable, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct SoilMoistureOutput <: AbstractOutputVariable
    name::String = "sm"
    unit::String = "1"
    long_name::String = "soil moisture"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)
    is_land::Bool = true    # 3D land variables have another vertical dimension
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 10
end

path(::SoilMoistureOutput, simulation) =
    simulation.prognostic_variables.land.soil_moisture

"""Defines netCDF output for a specific variable, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct SoilMoistureAvailabilityOutput <: AbstractOutputVariable
    name::String = "sma"
    unit::String = "1"
    long_name::String = "soil moisture availability"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 10
end

path(::SoilMoistureAvailabilityOutput, simulation) =
    simulation.diagnostic_variables.physics.land.soil_moisture_availability

"""Defines netCDF output for a specific variable, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct RiverRunoffOutput <: AbstractOutputVariable
    name::String = "roff"
    unit::String = "m"
    long_name::String = "river runoff"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 15
end

path(::RiverRunoffOutput, simulation) =
    simulation.diagnostic_variables.physics.land.river_runoff

"""Defines netCDF output for a specific variable, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct LandSeaMaskOutput <: AbstractOutputVariable
    name::String = "lsm"
    unit::String = "1"
    long_name::String = "land-sea mask (1=land, 0=sea)"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, false)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 10
end

path(::LandSeaMaskOutput, simulation) =
    simulation.model.land_sea_mask.mask

LandOutput() = (
    SoilTemperatureOutput(),
    SoilMoistureOutput(),
    SoilMoistureAvailabilityOutput(),
    RiverRunoffOutput(),
    LandSeaMaskOutput(),
)