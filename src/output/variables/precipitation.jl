"""Defines netCDF output for a specific variables, see `VorticityOutput` for details.
Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct ConvectivePrecipitationOutput <: AbstractOutputVariable
    name::String = "precip_conv"
    unit::String = "mm"
    long_name::String = "accumulated convective precipitation"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
    rate::AbstractOutputVariable = ConvectivePrecipitationRateOutput()
end

"""$(TYPEDSIGNATURES)
`output!` method for `variable`, see `output!(::NetCDFOutput, ::VorticityOutput, ...)` for details."""
function output!(
    output::NetCDFOutput,
    variable::ConvectivePrecipitationOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    # this is accumualted convective precipitation
    precip = output.grid2D
    (; precip_convection) = diagn.physics
    RingGrids.interpolate!(precip, precip_convection, output.interpolator)
    precip .*= 1000             # convert from [m] to [mm]

    round!(precip, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, i] = precip
    return nothing
end

# at finalize step postprocess the convective precipitation to get the rate
finalize!(output::NetCDFOutput, variable::ConvectivePrecipitationOutput, args...) = output!(output, variable.rate, variable)

abstract type AbstractRainRateOutputVariable <: AbstractOutputVariable end

"""Defines netCDF output for a specific variables, see `VorticityOutput` for details.
Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct ConvectivePrecipitationRateOutput <: AbstractRainRateOutputVariable
    name::String = "precip_conv_rate"
    unit::String = "mm/hr"
    long_name::String = "convective precipitation rate"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

function output!(
    output::NetCDFOutput,
    variable::AbstractRainRateOutputVariable,
    acc_variable::AbstractOutputVariable,
)
    # use .var to prevent Union{Missing, Float32} that NCDatasets uses
    accumulated = output.netcdf_file[acc_variable.name].var[:, :, :]

    # rate is defined as average precip since last output step, so first step is 0
    # convert from accumulated [m] to [mm/hr] rain rate over output time step (e.g. 6hours)
    s = Hour(1)/output.output_dt
    nx, ny = size(accumulated)
    rate = cat(zeros(eltype(accumulated), nx, ny), diff(accumulated, dims=3), dims=3)
    rate .*= s

    # DEFINE NEW NETCDF VARIABLE AND WRITE
    define_variable!(output.netcdf_file, variable, eltype(rate))

    output.netcdf_file[variable.name][:, :, :] = rate
    return nothing
end

"""Defines netCDF output for a specific variables, see `VorticityOutput` for details.
Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct LargeScalePrecipitationOutput <: AbstractOutputVariable
    name::String = "precip_cond"
    unit::String = "mm"
    long_name::String = "accumulated large-scale precipitation"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
    rate::AbstractOutputVariable = LargeScalePrecipitationRateOutput()
end

"""$(TYPEDSIGNATURES)
`output!` method for `variable`, see `output!(::NetCDFOutput, ::VorticityOutput, ...)` for details."""
function output!(
    output::NetCDFOutput,
    variable::LargeScalePrecipitationOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    # this is accumulated large-scale precipitation
    (; precip_large_scale) = diagn.physics
    precip = output.grid2D
    RingGrids.interpolate!(precip, precip_large_scale, output.interpolator)
    precip .*= 1000             # convert from [m] to [mm]

    round!(precip, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, i] = precip
    return nothing
end

# at finalize step postprocess the convective precipitation to get the rate
finalize!(output::NetCDFOutput, variable::LargeScalePrecipitationOutput, args...) = output!(output, variable.rate, variable)

"""Defines netCDF output for a specific variables, see `VorticityOutput` for details.
Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct LargeScalePrecipitationRateOutput <: AbstractRainRateOutputVariable
    name::String = "precip_cond_rate"
    unit::String = "mm/hr"
    long_name::String = "large-scale precipitation rate"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

# no need to redefine output! already defined for AbstractRainRateOutputVariable

## CLOUDS -------------

"""Defines netCDF output for a specific variables, see `VorticityOutput` for details.
Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct CloudTopOutput <: AbstractOutputVariable
    name::String = "cloud_top"
    unit::String = "m"
    long_name::String = "cloud top height"
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
    variable::CloudTopOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    cloud = output.grid2D
    (; cloud_top) = diagn.physics
    RingGrids.interpolate!(cloud, cloud_top, output.interpolator)

    round!(cloud, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, i] = cloud
    return nothing
end

# collect all in one for convenience
PrecipitationOutput() = (
    ConvectivePrecipitationOutput(),
    LargeScalePrecipitationOutput(),
    CloudTopOutput(),
)