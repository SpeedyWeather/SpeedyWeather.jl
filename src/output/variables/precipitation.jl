"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct ConvectivePrecipitationOutput{F, R} <: AbstractOutputVariable
    name::String = "precip_conv"
    unit::String = "mm"
    long_name::String = "accumulated convective precipitation"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 20
    transform::F = (x) -> 1000x                     # [m] to [mm]
    rate::R = ConvectivePrecipitationRateOutput()   # include here to be called at finalize!
end

path(::ConvectivePrecipitationOutput, simulation) =
    simulation.diagnostic_variables.physics.precip_convection

# at finalize step postprocess the convective precipitation to get the rate
finalize!(output::NetCDFOutput, variable::ConvectivePrecipitationOutput, args...) = output!(output, variable.rate, variable)

abstract type AbstractRateOutputVariable <: AbstractOutputVariable end

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct ConvectivePrecipitationRateOutput{F} <: AbstractRateOutputVariable
    name::String = "precip_conv_rate"
    unit::String = "mm/hr"
    long_name::String = "convective precipitation rate"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
    transform::F = (x) -> 1000x     # [m] to [mm]
end

"""$TYPEDSIGNATURES
Post-process the netCDF `output` file to convert accumulated precipitation/snow to
rates."""
function output!(
    output::NetCDFOutput,
    variable::AbstractRateOutputVariable,
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

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct LargeScalePrecipitationOutput{F, R} <: AbstractOutputVariable
    name::String = "precip_cond"
    unit::String = "mm"
    long_name::String = "accumulated large-scale precipitation"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 20
    transform::F = (x) -> 1000x                     # [m] to [mm]
    rate::R = LargeScalePrecipitationRateOutput()   # include here to be called at finalize!
end

path(::LargeScalePrecipitationOutput, simulation) =
    simulation.diagnostic_variables.physics.precip_large_scale

# at finalize step postprocess the convective precipitation to get the rate
finalize!(output::NetCDFOutput, variable::LargeScalePrecipitationOutput, args...) = output!(output, variable.rate, variable)

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct LargeScalePrecipitationRateOutput{F} <: AbstractRateOutputVariable
    name::String = "precip_cond_rate"
    unit::String = "mm/hr"
    long_name::String = "large-scale precipitation rate"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
    transform::F = (x) -> 1000x     # [m] to [mm]
end

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct ConvectiveSnowOutput{F, R} <: AbstractOutputVariable
    name::String = "snow_conv"
    unit::String = "mm"
    long_name::String = "accumulated convective snow"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 20
    transform::F = (x) -> 1000x                     # [m] to [mm]
    rate::R = ConvectiveSnowRateOutput()   # include here to be called at finalize!
end

path(::ConvectiveSnowOutput, simulation) =
    simulation.diagnostic_variables.physics.snow_convection

# at finalize step postprocess the convective snow to get the rate
finalize!(output::NetCDFOutput, variable::ConvectiveSnowOutput, args...) = output!(output, variable.rate, variable)

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct ConvectiveSnowRateOutput{F} <: AbstractRateOutputVariable
    name::String = "snow_conv_rate"
    unit::String = "mm/hr"
    long_name::String = "convective snow rate"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
    transform::F = (x) -> 1000x     # [m] to [mm]
end

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct LargeScaleSnowOutput{F, R} <: AbstractOutputVariable
    name::String = "snow_cond"
    unit::String = "mm"
    long_name::String = "accumulated large-scale snow"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 20
    transform::F = (x) -> 1000x                     # [m] to [mm]
    rate::R = LargeScaleSnowRateOutput()   # include here to be called at finalize!
end

path(::LargeScaleSnowOutput, simulation) =
    simulation.diagnostic_variables.physics.snow_large_scale

# at finalize step postprocess the convective snow precipitation to get the rate
finalize!(output::NetCDFOutput, variable::LargeScaleSnowOutput, args...) = output!(output, variable.rate, variable)

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct LargeScaleSnowRateOutput{F} <: AbstractRateOutputVariable
    name::String = "snow_cond_rate"
    unit::String = "mm/hr"
    long_name::String = "large-scale snow rate"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
    transform::F = (x) -> 1000x     # [m] to [mm]
end

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
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

path(::CloudTopOutput, simulation) =
    simulation.diagnostic_variables.physics.cloud_top

# collect all in one for convenience
PrecipitationOutput() = (
    ConvectivePrecipitationOutput(),
    LargeScalePrecipitationOutput(),
    # ConvectiveSnowOutput(),
    LargeScaleSnowOutput(),
    CloudTopOutput(),
)
