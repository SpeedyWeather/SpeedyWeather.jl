"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct RandomPattern <: AbstractOutputVariable
    name::String = "random_pattern"
    unit::String = "1"
    long_name::String = "random pattern"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::RandomPattern, simulation) =
    simulation.variables.grid.random_pattern
