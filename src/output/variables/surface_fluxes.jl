"""Defines netCDF output for a specific variables, see `VorticityOutput` for details.
Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct SurfaceFluxHeatOutput <: AbstractOutputVariable
    name::String = "surface_flux_heat"
    unit::String = "W/m^2"
    long_name::String = "Surface heat fluxes (positive down)"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

"""$(TYPEDSIGNATURES)
`output!` method for `variable`, see `output!(::NetCDFOutput, ::VorticityOutput, ...)` for details."""
function output!(
    output::NetCDFOutput,
    variable::SurfaceFluxHeatOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    flux = output.grid2D
    (; surface_flux_heat) = diagn.physics
    RingGrids.interpolate!(flux, surface_flux_heat, output.interpolator)

    round!(flux, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, i] = flux
    return nothing
end

"""Defines netCDF output for a specific variables, see `VorticityOutput` for details.
Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct SurfaceFluxHumidOutput <: AbstractOutputVariable
    name::String = "surface_flux_humid"
    unit::String = "kg/s/m^2"
    long_name::String = "Surface humidity fluxes (positive down)"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

"""$(TYPEDSIGNATURES)
`output!` method for `variable`, see `output!(::NetCDFOutput, ::VorticityOutput, ...)` for details."""
function output!(
    output::NetCDFOutput,
    variable::SurfaceFluxHumidOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    flux = output.grid2D
    (; surface_flux_humid) = diagn.physics
    RingGrids.interpolate!(flux, surface_flux_humid, output.interpolator)

    round!(flux, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, i] = flux
    return nothing
end

SurfaceFluxesOutput() = (
    SurfaceFluxHeatOutput(),
    SurfaceFluxHumidOutput(),
)