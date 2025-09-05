abstract type AbstractTracer <: AbstractModelComponent end

export Tracer

"""Define tracer through a `name::Symbol` (used as key in the tracer dictionaries)
and define it as `active` (=true default). `active=false` will (temporarily)
disable the time evolution of the tracer. Fields are
$(TYPEDFIELDS)"""
@kwdef mutable struct Tracer <: AbstractTracer
    name::Symbol
    active::Bool = true
end

# method to create a tracer from a positional symbol
Tracer(name::Symbol; kwargs...) = Tracer(; name, kwargs...)

const TRACER_DICT = Dict{Symbol, Tracer}

# low-level function to add tracers to the tracer dictionary
function add!(dict::TRACER_DICT, tracers::Tracer...)
    for tracer in tracers
        dict[tracer.name] = tracer
    end
    return dict
end

# calls the above as `model.tracers` is a `TRACER_DICT`
add!(model::AbstractModel, tracers::Tracer...) = add!(model.tracers, tracers...)

# adding a tracer to a simulation adds it to the model and the prognostic and diagnostic variables
function add!(simulation::AbstractSimulation, tracers::Tracer...)
    add!(simulation.model, tracers...)                      # first add tracer to model
    tracers_tuple = values(simulation.model.tracers)        # get all tracers from there
    add!(simulation.prognostic_variables, tracers_tuple...) # then add to ParticleVariables
    add!(simulation.diagnostic_variables, tracers_tuple...) # resets existing tracers to 0 if names match
    simulation.model.tracers
end

add!(vars::AbstractVariables, tracer_dict::TRACER_DICT) = add!(vars, values(tracer_dict)...)

export activate!

"""$(TYPEDSIGNATURES) Activate a deactivated (=frozen) tracers in a simulation, which is a setting in `simulation.model` only."""
activate!(simulation::AbstractSimulation, tracers::Tracer...) = activate!(simulation.model, tracers...)

# propagate activate! from simulation to model to the `model.tracers` dictionary
activate!(model::AbstractModel, tracers::Tracer...) = activate!(model.tracers, tracers...)

# pass on "value" to also allow _activate! to be used for decativation
activate!(dict::TRACER_DICT, tracers::Tracer...) = _activate!(dict, tracers..., value=true)
function _activate!(dict::TRACER_DICT, tracers::Tracer...; value::Bool=true)
    for tracer in tracers
        dict[tracer.name].active = value
    end
end

export deactivate!

"""$(TYPEDSIGNATURES) Deactivate a tracer in a simulation, which is a setting in `simulation.model` only."""
deactivate!(simulation::AbstractSimulation, tracers::Tracer...) = deactivate!(simulation.model, tracers...)
deactivate!(model::AbstractModel, tracers::Tracer...) = deactivate!(model.tracers, tracers...)
deactivate!(dict::TRACER_DICT, tracers::Tracer...) = _activate!(dict, tracers..., value=false)


"""$(TYPEDSIGNATURES) Delete a tracer from a simulation, deleted from the model and the variables."""
function Base.delete!(simulation::AbstractSimulation, tracer::Tracer)
    delete!(simulation.model, tracer)
    delete!(simulation.prognostic_variables, tracer)
    delete!(simulation.diagnostic_variables, tracer)
    simulation.model.tracers
end

# delete from tracer dictionary, identified by its `name::Symbol`
Base.delete!(model::AbstractModel, tracer::Tracer) = delete!(model.tracers, tracer.name)