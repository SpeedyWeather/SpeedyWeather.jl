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

export activate!

"""$(TYPEDSIGNATURES) Activate a deactivated (=frozen) tracers in a simulation, which is a setting in `simulation.model` only."""
activate!(simulation::AbstractSimulation, tracers::Tracer...) = activate!(simulation.model, tracers...)

# propagate activate! from simulation to model to the `model.tracers` dictionary
activate!(model::AbstractModel, tracers::Tracer...) = activate!(model.tracers, tracers...)

# pass on "value" to also allow _activate! to be used for decativation
activate!(dict::TRACER_DICT, tracers::Tracer...) = _activate!(dict, tracers..., value = true)
function _activate!(dict::TRACER_DICT, tracers::Tracer...; value::Bool = true)
    for tracer in tracers
        dict[tracer.name].active = value
    end
    return
end

export deactivate!

"""$(TYPEDSIGNATURES) Deactivate a tracer in a simulation, which is a setting in `simulation.model` only."""
deactivate!(simulation::AbstractSimulation, tracers::Tracer...) = deactivate!(simulation.model, tracers...)
deactivate!(model::AbstractModel, tracers::Tracer...) = deactivate!(model.tracers, tracers...)
deactivate!(dict::TRACER_DICT, tracers::Tracer...) = _activate!(dict, tracers..., value = false)

# delete from tracer dictionary, identified by its `name::Symbol`
Base.delete!(model::AbstractModel, tracer::Tracer) = delete!(model.tracers, tracer.name)

function variables(T::Tracer)
    return (
        PrognosticVariable(T.name, Spectral4D(2), desc = "$(T.name)", namespace = :tracers),    # 2 for 2 leapfrog steps
        GridVariable(T.name, Grid3D(), desc = "$(T.name)", namespace = :tracers),
        GridVariable(Symbol(T.name, :_prev), Grid3D(), desc = "$(T.name)", namespace = :tracers),
        TendencyVariable(T.name, Spectral3D(), desc = "Tendency of $(T.name)", namespace = :tracers),
        TendencyVariable(Symbol(T.name, :_grid), Grid3D(), desc = "Tendency of $(T.name)", namespace = :tracers),
        ScratchVariable(:a_grid, Grid3D(), desc = "Work array for tracer advection"),
        ScratchVariable(:b_grid, Grid3D(), desc = "Work array for tracer advection"),
    )
end

function variables(D::TRACER_DICT)
    return (variables(tracer) for tracer in values(D)) |> Iterators.flatten |> Tuple
end
