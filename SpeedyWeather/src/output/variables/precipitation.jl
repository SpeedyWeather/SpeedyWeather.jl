"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct ConvectiveRainOutput{F, R} <: AbstractOutputVariable
    name::String = "rain_conv"
    unit::String = "mm"
    long_name::String = "accumulated convective rain"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 20
    transform::F = (x) -> 1000x             # [m] to [mm]
    rate::R = ConvectiveRainRateOutput()    # include here to be called at finalize!
end

path(::ConvectiveRainOutput, simulation) =
    simulation.diagnostic_variables.physics.rain_convection

# at finalize step postprocess the convective rain to get the rate
finalize!(output::NetCDFOutput, variable::ConvectiveRainOutput, args...) = output!(output, variable.rate, variable)

abstract type AbstractRateOutputVariable <: AbstractOutputVariable end

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct ConvectiveRainRateOutput{F} <: AbstractRateOutputVariable
    name::String = "rain_conv_rate"
    unit::String = "mm/hr"
    long_name::String = "convective rain rate"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
    transform::F = (x) -> 1000x     # [m] to [mm]
end

"""$TYPEDSIGNATURES
Post-process the netCDF `output` file to convert accumulated precipitation rain/snow to
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
@kwdef mutable struct LargeScaleRainOutput{F, R} <: AbstractOutputVariable
    name::String = "rain_cond"
    unit::String = "mm"
    long_name::String = "accumulated large-scale rain"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 20
    transform::F = (x) -> 1000x             # [m] to [mm]
    rate::R = LargeScaleRainRateOutput()    # include here to be called at finalize!
end

path(::LargeScaleRainOutput, simulation) =
    simulation.diagnostic_variables.physics.rain_large_scale

# at finalize step postprocess the accumulated rain to get the rate
finalize!(output::NetCDFOutput, variable::LargeScaleRainOutput, args...) = output!(output, variable.rate, variable)

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct LargeScaleRainRateOutput{F} <: AbstractRateOutputVariable
    name::String = "rain_cond_rate"
    unit::String = "mm/hr"
    long_name::String = "large-scale rain rate"
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
    transform::F = (x) -> 1000x             # [m] to [mm]
    rate::R = LargeScaleSnowRateOutput()    # include here to be called at finalize!
end

path(::LargeScaleSnowOutput, simulation) =
    simulation.diagnostic_variables.physics.snow_large_scale

# at finalize step postprocess the convective snow to get the rate
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
@kwdef mutable struct TotalPrecipitationOutput{F, R} <: AbstractOutputVariable
    name::String = "precip_rate"
    unit::String = "kg/m^2/s"
    long_name::String = "total precipitation rate (rain+snow)"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 20
    transform::F = (x) -> x
end

path(::TotalPrecipitationOutput, simulation) =
    simulation.diagnostic_variables.physics.total_precipitation_rate

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
    ConvectiveRainOutput(),
    LargeScaleRainOutput(),
    # ConvectiveSnowOutput(),
    LargeScaleSnowOutput(),
    CloudTopOutput(),
)
