@kwdef mutable struct TracerOutput{F} <: AbstractOutputVariable
    name::String = "tracer1"
    unit::String = "?"
    long_name::String = "tracer1"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 15
    transform::F = (x) -> x
end

path(tracer::TracerOutput, simulation) = 
    simulation.diagnostic_variables.grid.tracers_grid[Symbol(tracer.name)]