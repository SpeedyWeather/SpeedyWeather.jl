"""Defines netCDF output for a specific variables, see `VorticityOutput` for details.
Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct OrographyOutput <: AbstractOutputVariable
    name::String = "orography"
    unit::String = "m"
    long_name::String = "orography"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, false)
    missing_value::Float64 = NaN
    compression_level::Int = DEFAULT_COMPRESSION_LEVEL
    shuffle::Bool = DEFAULT_SHUFFLE
    keepbits::Int = 15
end

path(::OrographyOutput, simulation) = simulation.model.orography.orography

BoundaryOutput() = (
    OrographyOutput(),
)