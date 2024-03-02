abstract type AbstractCallback end
const CALLBACK_DICT = Dict{Symbol,AbstractCallback}
const RANDSTRING_LENGTH = 4

export CallbackDict
function CallbackDict(callbacks::AbstractCallback...)
    callback_pairs = (Pair(Symbol("callback_"*randstring(RANDSTRING_LENGTH)),callback)
                        for callback in callbacks)
    CALLBACK_DICT(callback_pairs...)
end

CallbackDict() = CALLBACK_DICT()
CallbackDict(pairs::Pair{Symbol,<:AbstractCallback}...) = CALLBACK_DICT(pairs...)

function Base.show(io::IO,C::AbstractCallback)
    println(io,"$(typeof(C)) <: AbstractCallback")
    keys = propertynames(C)
    print_fields(io,C,keys)
end

# dummy callback
export NoCallback

"""Dummy callback that doesn't do anything."""
struct NoCallback <: AbstractCallback end
initialize!(::NoCallback,args...) = nothing     # executed once before the main time loop
callback!(::NoCallback,args...) = nothing       # executed after every time step
finish!(::NoCallback,args...) = nothing         # executed after main time loop finishes

# simply loop over dict of callbacks
for func in (:initialize!, :callback!, :finish!)
    @eval begin
        function $func(callbacks::CALLBACK_DICT,args...)
            for key in keys(callbacks)
                $func(callbacks[key],args...)
            end
        end
    end
end

export add!

"""
$(TYPEDSIGNATURES)
Add a callback to a Dict{String,AbstractCallback} dictionary. To be used like

    add!(model.callbacks,"mycallback",callback)

"""
function add!(D::CALLBACK_DICT,key::Symbol,callback::AbstractCallback)
    D[key] = callback
    return nothing
end

function add!(D::CALLBACK_DICT,key::String,callback::AbstractCallback)
    key_symbol = Symbol(key)
    @info "Callback keys are Symbols. String \"$key\" converted to Symbol :$key_symbol."
    add!(D,key_symbol,callback)
end

"""
$(TYPEDSIGNATURES)
Add a callback to a Dict{Symbol,AbstractCallback} dictionary without specifying the
key which is randomly created like callback_????. To be used like

    add!(model.callbacks,callback)

"""
function add!(D::CALLBACK_DICT,callback::AbstractCallback)
    key = Symbol("callback_"*Random.randstring(4))
    @info "$(typeof(callback)) callback added with key $key"
    D[key] = callback
    return nothing
end

# delete!(dict,key) already defined in Base

export GlobalSurfaceTemperatureCallback

"""
Callback that records the global mean surface temperature on every time step
$(TYPEDFIELDS)."""
Base.@kwdef mutable struct GlobalSurfaceTemperatureCallback{NF} <: AbstractCallback
    timestep_counter::Int = 0
    temp::Vector{NF} = zeros(DEFAULT_NF,0)
end

GlobalSurfaceTemperatureCallback(SG::SpectralGrid) = GlobalSurfaceTemperatureCallback{SG.NF}()

"""
$(TYPEDSIGNATURES)
Initializes callback.temp vector that records the global mean surface temperature on every time step.
Allocates vector of correct length (number of elements = total time steps plus one) and stores the
global surface temperature of the initial conditions"""
function initialize!(
    callback::GlobalSurfaceTemperatureCallback{NF},
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::ModelSetup,
) where NF
    callback.temp = Vector{NF}(undef,progn.clock.n_timesteps+1) # replace with vector of correct length
    callback.temp[1] = diagn.layers[diagn.nlev].temp_average[]  # set initial conditions
    callback.timestep_counter = 1                               # (re)set counter to 1
end

"""
$(TYPEDSIGNATURES)
Pulls the average temperature from the lowermost layer and stores it in the next
element of the callback.temp vector."""
function callback!(
    callback::GlobalSurfaceTemperatureCallback,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::ModelSetup,
)
    callback.timestep_counter += 1  
    i = callback.timestep_counter
    callback.temp[i] = diagn.layers[diagn.nlev].temp_average[]
end

# nothing to finish
finish!(::GlobalSurfaceTemperatureCallback,args...) = nothing