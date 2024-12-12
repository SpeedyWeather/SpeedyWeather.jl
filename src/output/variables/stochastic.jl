"""Defines netCDF output for a specific variables, see `VorticityOutput` for details.
Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct RandomPatternOutput <: AbstractOutputVariable
    name::String = "random_pattern"
    unit::String = "1"
    long_name::String = "random pattern"
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
    variable::RandomPatternOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    random_pattern = output.grid2D
    random_pattern_grid = diagn.grid.random_pattern
    RingGrids.interpolate!(random_pattern, random_pattern_grid, output.interpolator)

    round!(random_pattern, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, i] = random_pattern
    return nothing
end