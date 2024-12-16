abstract type AbstractTracer <: AbstractModelComponent end

export Tracer
@kwdef mutable struct Tracer <: AbstractTracer
    name::Symbol
    active::Bool = true
end

Tracer(name::Symbol; kwargs...) = Tracer(; name, kwargs...)

const TRACER_DICT = Dict{Symbol, AbstractTracer}

function add!(dict::TRACER_DICT, tracers::Tracer...)
    for tracer in tracers
        dict[tracer.name] = tracer
        @info "Tracer $(tracer.name) added."
    end
end

add!(model::AbstractModel, tracers::Tracer...) = add!(model.tracers, tracers...)
function add!(simulation::AbstractSimulation, tracers::Tracer...)
    add!(simulation.model, tracers...)                                      # first add tracer to model
    tracer_tuple = (tracer for (name, tracer) in simulation.model.tracers)  # get all tracers from there
    add!(simulation.prognostic_variables, tracer_tuple...)                  # then add to ParticleVariables
    add!(simulation.diagnostic_variables, tracer_tuple...)
end

function add!(vars::AbstractVariables, tracer_dict::TRACER_DICT)
    tracers = (tracer for (name, tracer) in tracer_dict)
    add!(vars, tracers...)
end