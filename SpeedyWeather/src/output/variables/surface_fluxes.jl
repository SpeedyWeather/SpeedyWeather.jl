"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct SurfaceSensibleHeatFluxOutput <: AbstractOutputVariable
    name::String = "shf"
    unit::String = "W/m^2"
    long_name::String = "Surface sensible heat flux (positive up)"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::SurfaceSensibleHeatFluxOutput, simulation) =
    simulation.variables.parameterizations.sensible_heat_flux

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct SurfaceHumidityFluxOutput <: AbstractOutputVariable
    name::String = "shuf"
    unit::String = "kg/s/m^2"
    long_name::String = "Surface humidity flux (positive up)"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::SurfaceHumidityFluxOutput, simulation) =
    simulation.variables.parameterizations.surface_humidity_flux

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct SurfaceLatentHeatFluxOutput <: AbstractOutputVariable
    name::String = "slf"
    unit::String = "W/m^2"
    long_name::String = "Surface latent heat flux (positive up)"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::SurfaceLatentHeatFluxOutput, simulation) =
    simulation.variables.parameterizations.surface_latent_heat_flux

SurfaceFluxesOutput() = (
    SurfaceSensibleHeatFluxOutput(),
    SurfaceHumidityFluxOutput(),
    # SurfaceLatentHeatFluxOutput(),    # don't output by default as it is proportional to humidity flux
)
