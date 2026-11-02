"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct Orography <: AbstractOutputVariable
    name::String = "orography"
    unit::String = "m"
    long_name::String = "orography"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, false)
    missing_value::Float64 = NaN
    compression_level::Int = DEFAULT_COMPRESSION_LEVEL
    shuffle::Bool = DEFAULT_SHUFFLE
    keepbits::Int = 15
end

path(::Orography, simulation) = simulation.model.orography.orography

Boundary() = (
    Orography(),
)
