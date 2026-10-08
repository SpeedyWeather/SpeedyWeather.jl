"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct CloudCondensateOutput <: AbstractOutputVariable
    name::String = "cloud_condensate"
    unit::String = "kg/kg"
    long_name::String = "cloud condensate (liquid + ice)"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::CloudCondensateOutput, simulation) = simulation.variables.grid.cloud_condensate

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct CloudFractionOutput <: AbstractOutputVariable
    name::String = "cloud_fraction"
    unit::String = "1"
    long_name::String = "cloud fraction"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::CloudFractionOutput, simulation) = simulation.variables.parameterizations.cloud_fraction

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct LiquidWaterPathOutput{F} <: AbstractOutputVariable
    name::String = "lwp"
    unit::String = "g/m^2"
    long_name::String = "liquid water path"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
    transform::F = (x) -> 1000x         # [kg/m²] to [g/m²]
end

path(::LiquidWaterPathOutput, simulation) = simulation.variables.parameterizations.liquid_water_path

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct IceWaterPathOutput{F} <: AbstractOutputVariable
    name::String = "iwp"
    unit::String = "g/m^2"
    long_name::String = "ice water path"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
    transform::F = (x) -> 1000x         # [kg/m²] to [g/m²]
end

path(::IceWaterPathOutput, simulation) = simulation.variables.parameterizations.ice_water_path

# collect all in one for convenience
CloudOutput() = (
    CloudCondensateOutput(),
    CloudFractionOutput(),
    LiquidWaterPathOutput(),
    IceWaterPathOutput(),
    OutgoingShortwaveClearSkyOutput(),      # for cloud radiative effects, with OneBandCloudyShortwave
    OutgoingLongwaveClearSkyOutput(),       # and OneBandCloudyLongwave
)
