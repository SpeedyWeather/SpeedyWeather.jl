"""Defines netCDF output for a specific variables, see `VorticityOutput` for details.
Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct SeaSurfaceTemperatureOutput <: AbstractOutputVariable
    name::String = "sst"
    unit::String = "degC"
    long_name::String = "sea surface temperature"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 10
end

"""$(TYPEDSIGNATURES)
`output!` method for `variable`, see `output!(::NetCDFOutput, ::VorticityOutput, ...)` for details."""
function output!(
    output::NetCDFOutput,
    variable::SeaSurfaceTemperatureOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    sst = output.grid2D
    sst_grid = progn.ocean.sea_surface_temperature
    RingGrids.interpolate!(sst, sst_grid, output.interpolator)

    sst .-= 273.15  # convert to ˚C

    round!(sst, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, i] = sst
    return nothing
end