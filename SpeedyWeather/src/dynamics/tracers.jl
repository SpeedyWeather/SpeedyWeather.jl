abstract type AbstractTracer <: AbstractModelComponent end

export Tracer

"""Define tracer through a `name::Symbol` (used as key in the tracer dictionaries)
and define it as `active` (=true default). `active=false` will (temporarily)
disable the time evolution of the tracer. Fields are
$(TYPEDFIELDS)"""
@kwdef mutable struct Tracer{B} <: AbstractTracer
    name::Symbol
    active::B = true
end

# method to create a tracer from a positional symbol
Tracer(name::Symbol; kwargs...) = Tracer{Bool}(; name, kwargs...)

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

function variables(T::Tracer, model::AbstractModel)
    nsteps = get_nsteps(model.time_stepping, model)
    pg = nsteps.prognostic_grid
    ps = nsteps.prognostic_spectral
    tg = nsteps.tendency_grid
    ts = nsteps.tendency_spectral
    return (
        PrognosticVariable(T.name, SpectralXYZT(ps), desc = "$(T.name)", namespace = :tracers),
        GridVariable(T.name, GridXYZT(pg), desc = "$(T.name)", namespace = :tracers),
        TendencyVariable(T.name, SpectralXYZT(ts), desc = "Tendency of $(T.name)", namespace = :tracers),
        TendencyVariable(T.name, GridXYZT(tg), desc = "Tendency of $(T.name)", namespace = :grid_tracers),
        ScratchVariable(:a, GridXYZ(), desc = "Scratch array", namespace = :grid),
        ScratchVariable(:b, GridXYZ(), desc = "Scratch array", namespace = :grid),
    )
end

function variables(D::TRACER_DICT, model::AbstractModel)
    return (variables(tracer, model) for tracer in values(D)) |> Iterators.flatten |> Tuple
end

"""$(TYPEDSIGNATURES)
Whether the tracer `name` is in the tracer dictionary `model.tracers` and active. A tracer deleted
after initialization is not, so its variables remain unchanged (frozen). A lookup in the dictionary,
which is not differentiated (inactive for Enzyme)."""
tracer_active(model::AbstractModel, name::Symbol) = tracer_active(model.tracers, name)
tracer_active(tracers::TRACER_DICT, name::Symbol) = (haskey(tracers, name) && tracers[name].active)::Bool

"""$(TYPEDSIGNATURES)
Call `f(Val(name), vars, model, args...)` for every tracer `name` that has variables in `vars` and
is active in `model.tracers` (see `tracer_active`). Unrolled at compile time over the tracer names
in the type of `vars`, so that no dictionary is iterated at runtime and `getfield(vars.….tracers, name)`
inside `f` is type-stable, as needed for Enzyme and Reactant. Does nothing for models without tracers."""
@generated function foreach_active_tracer(
        f::F, vars::Variables{Po, G, T}, model::AbstractModel, args...
    ) where {F, Po, G, T}
    calls = [
        :(tracer_active(model, $(QuoteNode(name))) && f(Val($(QuoteNode(name))), vars, model, args...))
            for name in _namespace_names(T, :tracers)
    ]
    return Expr(:block, calls..., :(return nothing))
end

"""$(TYPEDSIGNATURES)
Throw an `ArgumentError` if `model.tracers` contains tracers without variables in `vars`, i.e.
tracers added after the model was initialized, which the time stepping would otherwise skip."""
function check_tracers(vars::Variables, model::AbstractModel)
    names_with_variables = tracer_tendency_names(vars)
    names_without_variables = sort!([name for name in keys(model.tracers) if !(name in names_with_variables)])
    isempty(names_without_variables) || throw(
        ArgumentError(
            "Tracers $(Tuple(names_without_variables)) have no variables as they were added after the model " *
                "was initialized. Add all tracers before `initialize!(model)` or delete them with " *
                "`delete!(model, Tracer(name))`."
        )
    )
    return nothing
end
