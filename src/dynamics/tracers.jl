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
    end
    return dict
end

add!(model::AbstractModel, tracers::Tracer...) = add!(model.tracers, tracers...)
function add!(simulation::AbstractSimulation, tracers::Tracer...)
    add!(simulation.model, tracers...)                      # first add tracer to model
    tracers_tuple = values(simulation.model.tracers)        # get all tracers from there
    add!(simulation.prognostic_variables, tracers_tuple...) # then add to ParticleVariables
    add!(simulation.diagnostic_variables, tracers_tuple...) # resets existing tracers to 0 if names match
    simulation.model.tracers
end

add!(vars::AbstractVariables, tracer_dict::TRACER_DICT) = add!(vars, values(tracer_dict)...)

export activate!
activate!(simulation::AbstractSimulation, tracers::Tracer...) = activate!(simulation.model, tracers...)
activate!(model::AbstractModel, tracers::Tracer...) = activate!(model.tracers, tracers...)
activate!(dict::TRACER_DICT, tracers::Tracer...) = _activate!(dict, tracers..., value=true)
function _activate!(dict::TRACER_DICT, tracers::Tracer...; value::Bool=true)
    for tracer in tracers
        dict[tracer.name].active = value
    end
end

export deactivate!
deactivate!(simulation::AbstractSimulation, tracers::Tracer...) = deactivate!(simulation.model, tracers...)
deactivate!(model::AbstractModel, tracers::Tracer...) = deactivate!(model.tracers, tracers...)
deactivate!(dict::TRACER_DICT, tracers::Tracer...) = _activate!(dict, tracers..., value=false)

Base.delete!(model::AbstractModel, tracer::Tracer) = delete!(model.tracers, tracer.name)
function Base.delete!(simulation::AbstractSimulation, tracer::Tracer)
    delete!(simulation.model, tracer)
    delete!(simulation.prognostic_variables, tracer)
    delete!(simulation.diagnostic_variables, tracer)
    simulation.model.tracers
end