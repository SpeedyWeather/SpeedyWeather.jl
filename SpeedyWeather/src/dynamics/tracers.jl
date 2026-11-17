abstract type AbstractTracer <: AbstractModelComponent end

export Tracer

"""A tracer, advected with the flow. Tracers are defined as a NamedTuple through the `tracers`
keyword of the model constructor, the key being the name of the tracer, e.g.

    model = ShallowWaterModel(spectral_grid; tracers = (; abc = Tracer(), co2 = Tracer()))

`active = false` (default `true`) disables the time evolution of the tracer (it is frozen),
which can be changed at any time with `activate!(model, :abc)` and `deactivate!(model, :abc)`.
Fields are
$(TYPEDFIELDS)"""
@kwdef mutable struct Tracer <: AbstractTracer
    active::Bool = true
end

# `Tracer(:abc)` was the old way to define a tracer, give a migration hint
function Tracer(name::Symbol; kwargs...)
    throw(
        ArgumentError(
            "Tracers are not named any more, define them through the `tracers` keyword when constructing the model, " *
                "e.g. `ShallowWaterModel(spectral_grid; tracers = (; $name = Tracer()))`. (De)activate them with " *
                "`activate!(model, :$name)` or `deactivate!(model, :$name)`."
        )
    )
end

# names of the tracers in `model.tracers`, `tracers = nothing` is treated as no tracers
tracer_names(tracers::NamedTuple) = keys(tracers)
tracer_names(::Nothing) = ()

export activate!

"""$(TYPEDSIGNATURES) Activate the deactivated (=frozen) tracers `names` in a simulation,
which is a setting in `simulation.model` only."""
activate!(simulation::AbstractSimulation, names::Symbol...) = activate!(simulation.model, names...)

"""$(TYPEDSIGNATURES) Activate the deactivated (=frozen) tracers `names` in `model.tracers`."""
activate!(model::AbstractModel, names::Symbol...) = set_tracers_active!(model.tracers, true, names...)

export deactivate!

"""$(TYPEDSIGNATURES) Deactivate (=freeze) the tracers `names` in a simulation,
which is a setting in `simulation.model` only."""
deactivate!(simulation::AbstractSimulation, names::Symbol...) = deactivate!(simulation.model, names...)

"""$(TYPEDSIGNATURES) Deactivate (=freeze) the tracers `names` in `model.tracers`."""
deactivate!(model::AbstractModel, names::Symbol...) = set_tracers_active!(model.tracers, false, names...)

"""$(TYPEDSIGNATURES) Set the `active` flag of the tracers `names` in `tracers` to `value`,
throw an `ArgumentError` for names that are not tracers."""
function set_tracers_active!(tracers::NamedTuple, value::Bool, names::Symbol...)
    for name in names
        haskey(tracers, name) ||
            throw(ArgumentError("No tracer :$name in model.tracers, tracers are $(keys(tracers))"))
        tracers[name].active = value
    end
    return tracers
end

set_tracers_active!(::Nothing, value::Bool, names::Symbol...) = set_tracers_active!((;), value, names...)

# tracers are a NamedTuple, the key is the name of the tracer, as for greenhouse gases
variables(name_tracer::Pair{Symbol, Tracer}, model::AbstractModel) =
    variables(name_tracer.second, name_tracer.first, model)

function variables(::Tracer, name::Symbol, model::AbstractModel)
    nsteps = get_nsteps(model.time_stepping, model)
    pg = nsteps.prognostic_grid
    ps = nsteps.prognostic_spectral
    tg = nsteps.tendency_grid
    ts = nsteps.tendency_spectral
    return (
        PrognosticVariable(name, SpectralXYZT(ps), desc = "$name", namespace = :tracers),
        GridVariable(name, GridXYZT(pg), desc = "$name", namespace = :tracers),
        TendencyVariable(name, SpectralXYZT(ts), desc = "Tendency of $name", namespace = :tracers),
        TendencyVariable(name, GridXYZT(tg), desc = "Tendency of $name", namespace = :grid_tracers),
        ScratchVariable(:a, GridXYZ(), desc = "Scratch array", namespace = :grid),
        ScratchVariable(:b, GridXYZ(), desc = "Scratch array", namespace = :grid),
    )
end

"""$(TYPEDSIGNATURES)
Whether `name` is a tracer in `model.tracers` and active. With `name` a compile-time constant
(as in `foreach_active_tracer`) this is a field load, the lookup in the NamedTuple constant-folds.
Not differentiated (inactive for Enzyme)."""
tracer_active(model::AbstractModel, name::Symbol) = tracer_active(model.tracers, name)
tracer_active(tracers::NamedTuple, name::Symbol) = haskey(tracers, name) && getfield(tracers, name).active
tracer_active(::Nothing, name::Symbol) = false

"""$(TYPEDSIGNATURES)
Call `f(Val(name), vars, model, args...)` for every tracer `name` that has variables in `vars` and
is active in `model.tracers` (see `tracer_active`). Unrolled at compile time over the tracer names
in the type of `vars`, so that `getfield(vars.….tracers, name)` inside `f` is type-stable,
as needed for Enzyme and Reactant. Does nothing for models without tracers."""
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
Throw an `ArgumentError` if `model.tracers` is not a NamedTuple of `Tracer`s (or `nothing`), e.g. a
`Tuple` of tracers that would be silently ignored, or if its tracers are not those with variables in
`vars`, e.g. for a simulation assembled from the variables of a model with other tracers."""
function check_tracers(vars::Variables, model::AbstractModel)
    (; tracers) = model
    tracers isa Union{Nothing, NamedTuple{names, <:Tuple{Vararg{Tracer}}} where {names}} || throw(
        ArgumentError(
            "Tracers have to be a NamedTuple of `Tracer`s, e.g. `tracers = (; abc = Tracer(), co2 = Tracer())`, " *
                "but `model.tracers` is a $(typeof(tracers))."
        )
    )
    names_model = tracer_names(tracers)
    names_variables = tracer_tendency_names(vars)
    issetequal(names_model, names_variables) || throw(
        ArgumentError(
            "The tracers $names_model of the model are not the tracers $names_variables in the variables. " *
                "Create the simulation from the model with `initialize!(model)`."
        )
    )
    return nothing
end
