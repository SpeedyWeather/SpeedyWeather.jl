"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct ConvectiveRain{F, R} <: AbstractOutputVariable
    name::String = "rain_conv"
    unit::String = "mm"
    long_name::String = "accumulated convective rain"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 20
    transform::F = (x) -> 1000x             # [m] to [mm]
    rate::R = ConvectiveRainRate()    # include here to be called at finalize!
end

path(::ConvectiveRain, simulation) =
    simulation.variables.parameterizations.rain_convection

# at finalize step postprocess the convective rain to get the rate
finalize!(output::NetCDFOutput, variable::ConvectiveRain, args...) = output!(output, variable.rate, variable)

abstract type AbstractRateOutputVariable <: AbstractOutputVariable end

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct ConvectiveRainRate{F} <: AbstractRateOutputVariable
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
    s = Hour(1) / output.interval
    nx, ny = size(accumulated)
    rate = cat(zeros(eltype(accumulated), nx, ny), diff(accumulated, dims = 3), dims = 3)
    rate .*= s

    # DEFINE NEW NETCDF VARIABLE AND WRITE
    define_variable!(output.netcdf_file, variable, eltype(rate))

    output.netcdf_file[variable.name][:, :, :] = rate
    return nothing
end

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct LargeScaleRain{F, R} <: AbstractOutputVariable
    name::String = "rain_cond"
    unit::String = "mm"
    long_name::String = "accumulated large-scale rain"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 20
    transform::F = (x) -> 1000x             # [m] to [mm]
    rate::R = LargeScaleRainRate()    # include here to be called at finalize!
end

path(::LargeScaleRain, simulation) =
    simulation.variables.parameterizations.rain_large_scale

# at finalize step postprocess the accumulated rain to get the rate
finalize!(output::NetCDFOutput, variable::LargeScaleRain, args...) = output!(output, variable.rate, variable)

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct LargeScaleRainRate{F} <: AbstractRateOutputVariable
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

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct LargeScaleSnow{F, R} <: AbstractOutputVariable
    name::String = "snow_cond"
    unit::String = "mm"
    long_name::String = "accumulated large-scale snow"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 20
    transform::F = (x) -> 1000x             # [m] to [mm]
    rate::R = LargeScaleSnowRate()    # include here to be called at finalize!
end

path(::LargeScaleSnow, simulation) =
    simulation.variables.parameterizations.snow_large_scale

# at finalize step postprocess the convective snow to get the rate
finalize!(output::NetCDFOutput, variable::LargeScaleSnow, args...) = output!(output, variable.rate, variable)

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct LargeScaleSnowRate{F} <: AbstractRateOutputVariable
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

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct CloudTop <: AbstractOutputVariable
    name::String = "cloud_top"
    unit::String = "m"
    long_name::String = "cloud top height"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

# cloud top height [m] as diagnosed by the clouds, not the layer index `cloud_top`
path(::CloudTop, simulation) =
    simulation.variables.parameterizations.cloud_top_height

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct CloudCover <: AbstractOutputVariable
    name::String = "cloud_cover"
    unit::String = "1"
    long_name::String = "cloud cover"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::CloudCover, simulation) =
    simulation.variables.parameterizations.cloud_cover

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct ConvectiveSnow{F, R} <: AbstractOutputVariable
    name::String = "snow_conv"
    unit::String = "mm"
    long_name::String = "accumulated convective snow"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 20
    transform::F = (x) -> 1000x             # [m] to [mm]
    rate::R = ConvectiveSnowRate()    # include here to be called at finalize!
end

path(::ConvectiveSnow, simulation) =
    simulation.variables.parameterizations.snow_convection

# at finalize step postprocess the convective snow to get the rate
finalize!(output::NetCDFOutput, variable::ConvectiveSnow, args...) = output!(output, variable.rate, variable)

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct ConvectiveSnowRate{F} <: AbstractRateOutputVariable
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

# collect all in one for convenience
Precipitation() = (
    ConvectiveRain(),
    ConvectiveSnow(),
    LargeScaleRain(),
    LargeScaleSnow(),
    CloudTop(),
    CloudCover(),
)
