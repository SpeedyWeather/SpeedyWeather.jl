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

# All tracers are allocated in fused parents, stacked along the layer axis in the order of `model.tracers`:
# their spectral and grid variables in `vars.fused.tracers_spectral`, `vars.fused.tracers_grid`, and their spectral
# and grid tendencies followed by the fluxes u*A, v*A of their horizontal advection in
# `vars.fused.tracer_tendencies_spectral`, `vars.fused.tracer_tendencies_grid`, i.e. [tendencies; u fluxes; v fluxes]
# of N tracers with L layers each in layers 1:NL, NL+1:2NL, 2NL+1:3NL. So that all tracers can be processed
# together, see `fused_tracers`, with one transform of the tendencies and both fluxes of all tracers.
# The per-tracer variables (views of the parents) are used when they are processed one by one.
function variables(tracers::NamedTuple{names, <:Tuple{Vararg{Tracer}}}, model::AbstractModel) where {names}
    per_tracer = Tuple(Iterators.flatten(variables(tracer, name, model) for (name, tracer) in pairs(tracers)))
    u_fluxes = Tuple(Iterators.flatten(tracer_flux_variables(name, :u, model) for name in names))
    v_fluxes = Tuple(Iterators.flatten(tracer_flux_variables(name, :v, model) for name in names))
    return (per_tracer..., u_fluxes..., v_fluxes...)
end

function variables(::Tracer, name::Symbol, model::AbstractModel)
    nsteps = get_nsteps(model.time_stepping, model)
    pg = nsteps.prognostic_grid
    ps = nsteps.prognostic_spectral
    tg = nsteps.tendency_grid
    ts = nsteps.tendency_spectral
    return (
        PrognosticVariable(name, SpectralXYZT(ps), desc = "$name", namespace = :tracers, fuse = :tracers_spectral),
        GridVariable(name, GridXYZT(pg), desc = "$name", namespace = :tracers, fuse = :tracers_grid),
        TendencyVariable(name, SpectralXYZT(ts), desc = "Tendency of $name", namespace = :tracers, fuse = :tracer_tendencies_spectral),
        TendencyVariable(name, GridXYZT(tg), desc = "Tendency of $name", namespace = :grid_tracers, fuse = :tracer_tendencies_grid),
    )
end

# name of the flux `component` (:u or :v) of the horizontal advection of tracer `name`, e.g. `:abc_u_flux`
tracer_flux_name(name::Symbol, component::Symbol) = Symbol(name, :_, component, :_flux)

# the flux `component` (:u or :v, as u*A/coslat, v*A/coslat) of the horizontal advection of tracer `name` (=A),
# on the grid and in spectral space, members of the tracer tendencies' fused parents (with the same step dimension)
function tracer_flux_variables(name::Symbol, component::Symbol, model::AbstractModel)
    nsteps = get_nsteps(model.time_stepping, model)
    flux_name = tracer_flux_name(name, component)
    desc = "$component*$name/coslat flux of the horizontal advection of $name"
    return (
        DynamicsVariable(flux_name, SpectralXYZT(nsteps.tendency_spectral); desc, namespace = :tracers, fuse = :tracer_tendencies_spectral),
        DynamicsVariable(flux_name, GridXYZT(nsteps.tendency_grid); desc, namespace = :grid_tracers, fuse = :tracer_tendencies_grid),
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
Whether all tracers in `tracers` are active, `false` without tracers."""
all_tracers_active(tracers::NamedTuple) = !isempty(tracers) && all(tracer -> tracer.active, values(tracers))
all_tracers_active(::Nothing) = false

"""$(TYPEDSIGNATURES)
Default of the model option `fuse_tracers`: process all tracers together (fused) on GPU, where that is faster as
the transforms and kernels of a single tracer are too small to use the GPU, but one by one on CPU, where processing
them together is not faster (see `docs/dev/2026-11/2026-11-17-fused-tracers.md`)."""
default_fuse_tracers(architecture::AbstractArchitecture) = false
default_fuse_tracers(::GPU) = true

"""$(TYPEDSIGNATURES)
Whether all tracers are processed together (fused): one batched transform and one kernel per operation for all
tracers, stacked along the layer axis of the fused parents in `vars.fused` (see `variables(::NamedTuple{names,
<:Tuple{Vararg{Tracer}}}, ::AbstractModel)`). Requires `model.fuse_tracers` and all tracers to be active (a deactivated
tracer has to stay frozen); otherwise the tracers are processed one by one, see `foreach_active_tracer`.
Not differentiated (inactive for Enzyme)."""
fused_tracers(vars::Variables, model::AbstractModel) =
    haskey(vars.fused, :tracers_spectral) && model.fuse_tracers && all_tracers_active(model.tracers)

# Layers per call of a transform (or a spectral gradient) of slabs with `nlayers` layers each stacked along the layer
# axis (e.g. one slab per tracer): as many whole slabs as fit into the transform's largest batch `S.nlayers` (see
# `transform_batch` of the `SpectralGrid`), at least one
stacked_layers_per_call(nlayers::Integer, S::AbstractSpectralTransform) = max(1, S.nlayers ÷ nlayers) * nlayers

# Layer ranges of the chunks of `nlayers_total` stacked layers with `layers_per_call` layers each (the last one shorter)
stacked_chunks(nlayers_total::Integer, layers_per_call::Integer) =
    (layer:min(layer + layers_per_call - 1, nlayers_total) for layer in 1:layers_per_call:nlayers_total)

"""$(TYPEDSIGNATURES)
`transform!` of `input` to `output`, both holding slabs of `nlayers` layers stacked along the layer axis (e.g. one
slab per tracer), in chunks of whole slabs that fit into the transform's largest batch `S.nlayers` (see
`transform_batch` of the `SpectralGrid`), with one call if all slabs fit."""
function transform_stacked!(output, input, nlayers::Integer, scratch_memory, S::AbstractSpectralTransform)
    layers_per_call = stacked_layers_per_call(nlayers, S)
    size(input, 2) <= layers_per_call && return transform!(output, input, scratch_memory, S)
    for layers in stacked_chunks(size(input, 2), layers_per_call)
        transform!(wrapped_view(output, :, layers), wrapped_view(input, :, layers), scratch_memory, S)
    end
    return output
end

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
