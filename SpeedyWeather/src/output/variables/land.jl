"""Defines netCDF output for a specific variable, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct SoilTemperature{F} <: AbstractOutputVariable
    name::String = "st"
    unit::String = "degC"
    long_name::String = "soil temperature"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)
    is_land::Bool = true    # 3D land variables have another vertical dimension
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 10
    transform::F = (x) -> x - oftype(x, 273.15)
end

path(::SoilTemperature, simulation) =
    simulation.variables.prognostic.land.soil_temperature

"""Defines netCDF output for a specific variable, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct SoilMoisture <: AbstractOutputVariable
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

path(::SoilMoisture, simulation) =
    simulation.variables.prognostic.land.soil_moisture

"""Defines netCDF output for a specific variable, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct SoilMoistureAvailability <: AbstractOutputVariable
    name::String = "sma"
    unit::String = "1"
    long_name::String = "soil moisture availability"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 10
end

path(::SoilMoistureAvailability, simulation) =
    simulation.variables.parameterizations.land.soil_moisture_availability

"""Defines netCDF output for a specific variable, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct RiverRunoff <: AbstractOutputVariable
    name::String = "roff"
    unit::String = "m"
    long_name::String = "river runoff"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 15
end

path(::RiverRunoff, simulation) =
    simulation.variables.parameterizations.land.river_runoff

"""Defines netCDF output for a specific variable, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct SnowMelt <: AbstractOutputVariable
    name::String = "smelt"
    unit::String = "kg/m^2/s"
    long_name::String = "snow melt+runoff rate"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 12
end

path(::SnowMelt, simulation) =
    simulation.variables.parameterizations.land.snow_melt_rate

"""Defines netCDF output for a specific variable, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct SnowDepth <: AbstractOutputVariable
    name::String = "sd"
    unit::String = "m"
    long_name::String = "snow depth (equivalent liquid water height)"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 10
end

path(::SnowDepth, simulation) =
    simulation.variables.prognostic.land.snow_depth

"""Defines netCDF output for a specific variable, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct LandSeaMask <: AbstractOutputVariable
    name::String = "lsm"
    unit::String = "1"
    long_name::String = "land-sea mask (1=land, 0=sea)"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, false)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 10
end

path(::LandSeaMask, simulation) =
    simulation.model.land_sea_mask.land_fraction

Land() = (
    SoilTemperature(),
    SoilMoisture(),
    SoilMoistureAvailability(),
    RiverRunoff(),
    SnowMelt(),
    SnowDepth(),
    LandSeaMask(),
)

"""$(TYPEDSIGNATURES)
Output variables for the state of a Terrarium land model in coupled
Terrarium-SpeedyWeather simulations.
Call as `OutputVariables.Terrarium(terrarium_model)` to create output variables for all
prognostic and auxiliary variables of the Terrarium model, or as
`OutputVariables.Terrarium(terrarium_model, :variable_name)` for a single one. Output
variables are added to a model's output via `add!(model, OutputVariables.Terrarium(terrarium_model)...)`.
Supported by `NetCDFOutput`, and by `ZarrOutput` once Zarr.jl is loaded."""
function Terrarium end
