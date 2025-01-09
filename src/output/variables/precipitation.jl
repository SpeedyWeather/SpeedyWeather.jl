"""Defines netCDF output for a specific variables, see `VorticityOutput` for details.
Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct ConvectivePrecipitationOutput{F} <: AbstractOutputVariable
    name::String = "precip_conv"
    unit::String = "mm"
    long_name::String = "accumulated convective precipitation"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
    transform::F = (x) -> 1000x     # [m] to [mm]
end

path(::ConvectivePrecipitationOutput, simulation) =
    simulation.diagnostic_variables.physics.precip_convection

"""Defines netCDF output for a specific variables, see `VorticityOutput` for details.
Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct ConvectivePrecipitationRateOutput{F} <: AbstractRainRateOutputVariable
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

path(::ConvectivePrecipitationOutput, simulation) =
    simulation.diagnostic_variables.physics.precip_rate_convection

# function output!(
#     output::NetCDFOutput,
#     variable::AbstractRainRateOutputVariable,
#     acc_variable::AbstractOutputVariable,
# )
#     # use .var to prevent Union{Missing, Float32} that NCDatasets uses
#     accumulated = output.netcdf_file[acc_variable.name].var[:, :, :]

#     # rate is defined as average precip since last output step, so first step is 0
#     # convert from accumulated [m] to [mm/hr] rain rate over output time step (e.g. 6hours)
#     s = Hour(1)/output.output_dt
#     nx, ny = size(accumulated)
#     rate = cat(zeros(eltype(accumulated), nx, ny), diff(accumulated, dims=3), dims=3)
#     rate .*= s

#     # DEFINE NEW NETCDF VARIABLE AND WRITE
#     define_variable!(output.netcdf_file, variable, eltype(rate))

#     output.netcdf_file[variable.name][:, :, :] = rate
#     return nothing
# end

"""Defines netCDF output for a specific variables, see `VorticityOutput` for details.
Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct LargeScalePrecipitationOutput{F} <: AbstractOutputVariable
    name::String = "precip_cond"
    unit::String = "mm"
    long_name::String = "accumulated large-scale precipitation"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
    transform::F = (x) -> 1000x     # [m] to [mm]
end

path(::LargeScalePrecipitationOutput, simulation) =
    simulation.diagnostic_variables.physics.precip_large_scale

"""Defines netCDF output for a specific variables, see `VorticityOutput` for details.
Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct LargeScalePrecipitationRateOutput{F} <: AbstractRainRateOutputVariable
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

path(::LargeScalePrecipitationRateOutput, simulation) =
    simulation.diagnostic_variables.physics.precip_rate_large_scale

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

path(::CloudTopOutput, simulation) =
    simulation.diagnostic_variables.physics.cloud_top

# collect all in one for convenience
PrecipitationOutput() = (
    ConvectivePrecipitationOutput(),
    # ConvectivePrecipitationRateOutput(),
    LargeScalePrecipitationOutput(),
    # LargeScalePrecipitationRateOutput(),
    CloudTopOutput(),
)