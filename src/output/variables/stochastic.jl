"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct RandomPatternOutput <: AbstractOutputVariable
    name::String = "random_pattern"
    unit::String = "1"
    long_name::String = "random pattern"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::RandomPatternOutput, simulation) =
    simulation.diagnostic_variables.grid.random_pattern