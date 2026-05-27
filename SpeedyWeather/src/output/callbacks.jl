abstract type AbstractCallback <: AbstractModelComponent end
const CALLBACK_DICT = Dict{Symbol, AbstractCallback}
const RANDSTRING_LENGTH = 4

export CallbackDict
function CallbackDict(callbacks::AbstractCallback...)
    callback_pairs = (
        Pair(Symbol("callback_" * randstring(RANDSTRING_LENGTH)), callback)
            for callback in callbacks
    )
    return CALLBACK_DICT(callback_pairs...)
end

"""$(TYPEDSIGNATURES)
Empty Callback dictionary generator."""
CallbackDict() = CALLBACK_DICT()

"""$(TYPEDSIGNATURES)
Create Callback dictionary like normal dictionaries."""
CallbackDict(pairs::Pair{Symbol, <:AbstractCallback}...) = CALLBACK_DICT(pairs...)

# dummy callback
export NoCallback

"""Dummy callback that doesn't do anything."""
struct NoCallback <: AbstractCallback end
initialize!(::NoCallback, args...) = nothing     # executed once before the main time loop
callback!(::NoCallback, args...) = nothing       # executed after every time step
finalize!(::NoCallback, args...) = nothing       # executed after main time loop finishes

# simply loop over dict of callbacks, passing the simulation through to each callback
for func in (:initialize!, :callback!, :finalize!)
    @eval begin
        function $func(callbacks::CALLBACK_DICT, simulation::AbstractSimulation)
            for key in keys(callbacks)
                $func(callbacks[key], simulation)
            end
            return
        end
    end
end

export add!
"""
$(TYPEDSIGNATURES)
Add a or several callbacks to a Dict{String, AbstractCallback} dictionary. To be used like

    add!(simulation.callbacks, :my_callback => callback)
    add!(simulation.callbacks, :my_callback1 => callback, :my_callback2 => other_callback)
"""
function add!(D::CALLBACK_DICT, key_callbacks::Pair{Symbol, <:AbstractCallback}...)
    for key_callback in key_callbacks
        key = key_callback.first
        callback = key_callback.second
        D[key] = callback
    end
    return D
end

"""
$(TYPEDSIGNATURES)
Add a or several callbacks to a `simulation::AbstractSimulation`."""
add!(sim::AbstractSimulation, key_callbacks::Pair{Symbol, <:AbstractCallback}...) =
    add!(sim.callbacks, key_callbacks...)
add!(D::CALLBACK_DICT, key::Symbol, callback::AbstractCallback) = add!(D, Pair(key, callback))
add!(sim::AbstractSimulation, key::Symbol, callback::AbstractCallback) =
    add!(sim.callbacks, Pair(key, callback))


# also with string but flag conversion
function add!(D::CALLBACK_DICT, key::String, callback::AbstractCallback)
    key_symbol = Symbol(key)
    @warn "Callback keys are Symbols. String \"$key\" converted to Symbol :$key_symbol."
    return add!(D, key_symbol, callback)
end

"""
$(TYPEDSIGNATURES)
Add a or several callbacks to a Dict{Symbol, AbstractCallback} dictionary without specifying the
key which is randomly created like callback_????. To be used like

    add!(simulation.callbacks, callback)
    add!(simulation.callbacks, callback1, callback2)."""
function add!(D::CALLBACK_DICT, callbacks::AbstractCallback...; verbose = true)
    for callback in callbacks
        key = Symbol("callback_" * Random.randstring(4))
        verbose && @info "$(typeof(callback)) callback added with key $key"
        add!(D, key => callback)
    end
    return D
end

"""
$(TYPEDSIGNATURES)
Add a or several callbacks to a `simulation::AbstractSimulation` without specifying the
key which is randomly created like callback_????. To be used like

    add!(simulation, callback)
    add!(simulation, callback1, callback2)."""
add!(sim::AbstractSimulation, callbacks::AbstractCallback...) =
    add!(sim.callbacks, callbacks..., verbose = sim.feedback.verbose)

# adding via tuple splats the tuple
add!(sim::AbstractSimulation, tuple::Tuple) = add!(sim, tuple...)

# delete!(dict, key) already defined in Base

export GlobalSurfaceTemperatureCallback

"""
Callback that records the global mean surface temperature on every time step.
$(TYPEDFIELDS)"""
Base.@kwdef mutable struct GlobalSurfaceTemperatureCallback{NF} <: AbstractCallback
    timestep_counter::Int = 0
    temperature::Vector{NF} = zeros(DEFAULT_NF, 0)
end

GlobalSurfaceTemperatureCallback(SG::SpectralGrid) = GlobalSurfaceTemperatureCallback{SG.NF}()

"""
$(TYPEDSIGNATURES)
Initializes callback.temperature vector that records the global mean surface temperature on every time step.
Allocates vector of correct length (number of elements = total time steps plus one) and stores the
global surface temperature of the initial conditions"""
function initialize!(
        callback::GlobalSurfaceTemperatureCallback{NF},
        simulation::AbstractSimulation,
    ) where {NF}
    vars, model = simulation.variables, simulation.model
    callback.temperature = Vector{NF}(undef, vars.prognostic.clock.n_timesteps + 1)    # replace with vector of correct length
    nlayers = model.geometry.nlayers
    callback.temperature[1] = vars.grid.temp_average[nlayers]       # set initial conditions
    callback.timestep_counter = 1                                   # (re)set counter to 1
    return nothing
end

"""
$(TYPEDSIGNATURES)
Pulls the average temperature from the lowermost layer and stores it in the next
element of the callback.temperature vector."""
function callback!(
        callback::GlobalSurfaceTemperatureCallback,
        simulation::AbstractSimulation,
    )
    vars, model = simulation.variables, simulation.model
    callback.timestep_counter += 1
    i = callback.timestep_counter
    nlayers = model.geometry.nlayers
    callback.temperature[i] = vars.grid.temp_average[nlayers]
    return nothing
end

# nothing to finalize
finalize!(::GlobalSurfaceTemperatureCallback, args...) = nothing
