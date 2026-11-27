"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct ConvectiveRain{F} <: AbstractOutputVariable
    name::String = "rain_conv"
    unit::String = "mm"
    long_name::String = "convective rain accumulated since last output"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 20
    transform::F = (x) -> 1000x             # [m] to [mm]
end

path(::ConvectiveRain, simulation) =
    simulation.variables.output.rain_convection

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct LargeScaleRain{F} <: AbstractOutputVariable
    name::String = "rain_cond"
    unit::String = "mm"
    long_name::String = "large-scale rain accumulated since last output"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 20
    transform::F = (x) -> 1000x             # [m] to [mm]
end

path(::LargeScaleRain, simulation) =
    simulation.variables.output.rain_large_scale

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct LargeScaleSnow{F} <: AbstractOutputVariable
    name::String = "snow_cond"
    unit::String = "mm"
    long_name::String = "large-scale snow accumulated since last output"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 20
    transform::F = (x) -> 1000x             # [m] to [mm]
end

path(::LargeScaleSnow, simulation) =
    simulation.variables.output.snow_large_scale

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
@kwdef mutable struct ConvectiveSnow{F} <: AbstractOutputVariable
    name::String = "snow_conv"
    unit::String = "mm"
    long_name::String = "convective snow accumulated since last output"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 20
    transform::F = (x) -> 1000x             # [m] to [mm]
end

path(::ConvectiveSnow, simulation) =
    simulation.variables.output.snow_convection

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct CAPE <: AbstractOutputVariable
    name::String = "CAPE"
    unit::String = "J/kg"
    long_name::String = "convective available potential energy"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 15
end

path(::CAPE, simulation) =
    simulation.variables.parameterizations.CAPE

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct CAPEPerLayer <: AbstractOutputVariable
    name::String = "CAPE_per_layer"
    unit::String = "J/kg"
    long_name::String = "convective available potential energy per vertical layer"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 10
end

path(::CAPEPerLayer, simulation) =
    simulation.variables.parameterizations.CAPE_per_layer

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct PT <: AbstractOutputVariable
    name::String = "PT"
    unit::String = "mm/s"
    long_name::String = "precipitation rate due to warming"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 15
end

path(::PT, simulation) =
    simulation.variables.parameterizations.PT

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct Pq <: AbstractOutputVariable
    name::String = "Pq"
    unit::String = "mm/s"
    long_name::String = "precipitation rate due to drying"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 15
end

path(::Pq, simulation) =
    simulation.variables.parameterizations.Pq

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct PTPerLayer <: AbstractOutputVariable
    name::String = "PT_per_layer"
    unit::String = "mm/s"
    long_name::String = "precipitation rate due to warming per vertical layer"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 10
end

path(::PTPerLayer, simulation) =
    simulation.variables.parameterizations.PT_per_layer

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct PqPerLayer <: AbstractOutputVariable
    name::String = "Pq_per_layer"
    unit::String = "mm/s"
    long_name::String = "precipitation rate due to drying per vertical layer"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 10
end

path(::PqPerLayer, simulation) =
    simulation.variables.parameterizations.Pq_per_layer

# collect all in one for convenience
Precipitation() = (
    ConvectiveRain(),
    ConvectiveSnow(),
    LargeScaleRain(),
    LargeScaleSnow(),
    CloudTop(),
    CloudCover(),
    PT(),
    PTPerLayer(),
    Pq(),
    PqPerLayer(),
    CAPE(),
    CAPEPerLayer()
)
