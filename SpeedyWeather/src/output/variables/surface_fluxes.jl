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
    simulation.diagnostic_variables.physics.sensible_heat_flux

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
    simulation.diagnostic_variables.physics.surface_humidity_flux

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
    simulation.diagnostic_variables.physics.surface_latent_heat_flux

SurfaceFluxesOutput() = (
    SurfaceSensibleHeatFluxOutput(),
    # SurfaceHumidityFluxOutput(),      # don't output by default as it is proportional to latent heat flux
    SurfaceLatentHeatFluxOutput(),
)