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
    simulation.variables.parameterizations.rain_convection

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
    s = Hour(1) / output.interval
    nx, ny = size(accumulated)
    rate = cat(zeros(eltype(accumulated), nx, ny), diff(accumulated, dims = 3), dims = 3)
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
    simulation.variables.parameterizations.rain_large_scale

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
    simulation.variables.parameterizations.snow_large_scale

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

# cloud top is stored as layer index, converted to height via geopotential in output!
path(::CloudTopOutput, simulation) =
    simulation.variables.parameterizations.cloud_top

"""$(TYPEDSIGNATURES)
Output the cloud top height [m] (above sea level) from the cloud top layer index,
using the geopotential on that layer, `Φ/g`. Columns without clouds (index `nlayers + 1`)
are written as 0 m."""
function output!(
        output::AbstractOutput,
        variable::CloudTopOutput,
        simulation::AbstractSimulation,
    )
    # escape immediately after first call if variable doesn't have a time dimension
    ~hastime(variable) && output.output_counter > 1 && return nothing

    var = path_or_nothing(variable, simulation)
    isnothing(var) && return nothing       # silently escape early if variable is not defined

    # index-based lookup of geopotential per column, do on CPU (output is written from CPU anyway)
    cloud_top = on_architecture(CPU(), var)
    geopotential = on_architecture(CPU(), simulation.variables.dynamics.geopotential)
    g = simulation.model.planet.gravity
    nlayers = size(geopotential, 2)

    cloud_top_height = similar(cloud_top)
    for ij in eachindex(cloud_top, cloud_top_height)
        k = round(Int, cloud_top[ij])
        cloud_top_height[ij] = 1 <= k <= nlayers ? geopotential[ij, k] / g : 0
    end

    # interpolate 2D/3D variables
    cloud_top_height_output = output.field2D
    interpolate_output!(output, cloud_top_height_output, cloud_top_height)

    if hasproperty(variable, :keepbits)     # round mantissabits for compression
        round!(cloud_top_height_output, variable.keepbits)
    end

    write_array!(output, variable, cloud_top_height_output)
    return nothing
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
    transform::F = (x) -> 1000x             # [m] to [mm]
    rate::R = ConvectiveSnowRateOutput()    # include here to be called at finalize!
end

path(::ConvectiveSnowOutput, simulation) =
    simulation.variables.parameterizations.snow_convection

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

# collect all in one for convenience
PrecipitationOutput() = (
    ConvectiveRainOutput(),
    ConvectiveSnowOutput(),
    LargeScaleRainOutput(),
    LargeScaleSnowOutput(),
    CloudTopOutput(),
)
