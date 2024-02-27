abstract type AbstractCallback end

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

# simply loop over vector of callbacks
for func in (:initialize!, :callback!, :finish!)
    @eval begin
        function $func(callbacks::Vector{<:AbstractCallback},args...)
            for callback in callbacks
                $func(callback,args...)
            end
        end
    end
end

# define to make append!(::AbstractVector{<:AbstractCallback},::AbstractCallback) possible
Base.length(::SpeedyWeather.AbstractCallback) = 1
Base.iterate(c::SpeedyWeather.AbstractCallback) = (c,nothing)

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