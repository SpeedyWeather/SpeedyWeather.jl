"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct SurfaceSensibleHeatFlux <: AbstractOutputVariable
    name::String = "shf"
    unit::String = "W/m^2"
    long_name::String = "Surface sensible heat flux (positive up)"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::SurfaceSensibleHeatFlux, simulation) =
    simulation.variables.parameterizations.sensible_heat_flux

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct SurfaceHumidityFlux <: AbstractOutputVariable
    name::String = "shuf"
    unit::String = "kg/s/m^2"
    long_name::String = "Surface humidity flux (positive up)"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::SurfaceHumidityFlux, simulation) =
    simulation.variables.parameterizations.surface_humidity_flux

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct SurfaceLatentHeatFlux <: AbstractOutputVariable
    name::String = "slf"
    unit::String = "W/m^2"
    long_name::String = "Surface latent heat flux (positive up)"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::SurfaceLatentHeatFlux, simulation) =
    simulation.variables.parameterizations.surface_latent_heat_flux

SurfaceFluxes() = (
    SurfaceSensibleHeatFlux(),
    SurfaceHumidityFlux(),
    # SurfaceLatentHeatFlux(),    # don't output by default as it is proportional to humidity flux
)
