"""Defines netCDF output for a specific variables, see `VorticityOutput` for details.
Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct OrographyOutput <: AbstractOutputVariable
    name::String = "orography"
    unit::String = "m"
    long_name::String = "orography"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, false)
    missing_value::Float64 = NaN
    compression_level::Int = DEFAULT_COMPRESSION_LEVEL
    shuffle::Bool = DEFAULT_SHUFFLE
    keepbits::Int = 15
end

"""$(TYPEDSIGNATURES)
`output!` method for `variable`, see `output!(::NetCDFOutput, ::VorticityOutput, ...)` for details."""
function output!(
    output::NetCDFOutput,
    variable::OrographyOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    # escape immediately when initialization counter > 1 to not write orography again
    output.output_counter > 1 || return nothing

    orog = output.grid2D
    (; orography) = model.orography
    RingGrids.interpolate!(orog, orography, output.interpolator)

    round!(orog, variable.keepbits)
    output.netcdf_file[variable.name][:, :] = orog
    return nothing
end

BoundaryOutput() = (
    OrographyOutput(),
)