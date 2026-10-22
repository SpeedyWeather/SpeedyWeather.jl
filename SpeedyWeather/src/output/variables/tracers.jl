@kwdef mutable struct Tracer{F} <: AbstractOutputVariable
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

# as tracers are identified by their name::Symbol, allow this too
Tracer(name::Symbol; kwargs...) =
    Tracer(; name = string(name), long_name = string(name), kwargs...)

path(tracer::Tracer, simulation) =
    simulation.variables.grid.tracers[Symbol(tracer.name)]
