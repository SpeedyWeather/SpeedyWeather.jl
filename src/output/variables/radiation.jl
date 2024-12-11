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

"""$(TYPEDSIGNATURES)
`output!` method for `variable`, see `output!(::NetCDFOutput, ::VorticityOutput, ...)` for details."""
function output!(
    output::NetCDFOutput,
    variable::OutgoingLongwaveRadiationOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    olr = output.grid2D
    (; outgoing_longwave_radiation) = diagn.physics
    RingGrids.interpolate!(olr, outgoing_longwave_radiation, output.interpolator)

    round!(olr, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, i] = olr
    return nothing
end

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

"""$(TYPEDSIGNATURES)
`output!` method for `variable`, see `output!(::NetCDFOutput, ::VorticityOutput, ...)` for details."""
function output!(
    output::NetCDFOutput,
    variable::OutgoingShortwaveRadiationOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    osr = output.grid2D
    (; outgoing_shortwave_radiation) = diagn.physics
    RingGrids.interpolate!(osr, outgoing_shortwave_radiation, output.interpolator)

    round!(osr, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, i] = osr
    return nothing
end

RadiationOutput() = (OutgoingLongwaveRadiationOutput(), OutgoingShortwaveRadiationOutput())