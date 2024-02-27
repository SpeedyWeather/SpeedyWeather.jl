# Callbacks

SpeedyWeather.jl implements a callback system to let users include a flexible piece of code
into the time stepping. You can think about the main time loop *calling back* to check whether
anything else should be done before continuing with the next time step. The callback system
here is called *after* the time step only (excluding calls at with `initialize!` and `finish!`),
we currently do not implement other callsites.

Callbacks are mainly introduced for diagnostic purposes, meaning that they do not influence
the simulation, and access the prognostic variables and the model components in a read-only
fashion. However, you may define a callback that changes the orography during the simulation.
In general, one has to keep the general order of executions during a time step in mind
(valid for all models)

1. set tendencies to zero
2. compute parameterizations, forcing, or drag terms. Accumulate tendencies.
3. compute dynamics, accumulate tendencies.
4. time stepping
5. output
6. callbacks

This means that, at the current callsite, a callback can read the tendencies but writing
into it would be overwritten by the zeroing of the tendencies in 1. anyway. At the moment,
if a callback wants to implement an additional tendency then it currently should be
implemented as a parameterization, forcing or drag term. 

## Defining a callback

You can (and are encouraged!) to write your own callbacks to diagnose SpeedyWeather simulations.
Let us implement a `StormChaser` callback, recording the highest surface wind speed
on every time step, that we want to use to illustrate how a callback needs
to be defined.

Every custom callback needs to be defined as a (`mutable`) `struct`, subtype of `AbstractCallback`,
i.e. `struct` or `mutable struct CustomCallback <: SpeedyWeather.AbstractCallback`. In our case, this is

```@example callback
using SpeedyWeather

Base.@kwdef mutable struct StormChaser{NF} <: SpeedyWeather.AbstractCallback
    timestep_counter::Int = 0
    maximum_surface_wind_speed::Vector{NF} = [0]
end

# Generator function
StormChaser(SG::SpectralGrid) = StormChaser{SG.NF}()
```

We decide to have a field `timestep_counter` in the callback that allows us to track
the number of times the callback was called to create a time series of our highest
surface wind speeds. The actual `maximum_surface_wind_speed` is then a vector
of a given type `NF` (= number format), which is where we'll write into. Both
are initialised with zeros. We also add a generator function, similar as to
many other components in SpeedyWeather that just pulls the number format from
the `SpectralGrid` object.

Now every callback needs to extend three methods

1. `initialize!`, called once before the main time loop starts
2. `callback!`, called after every time step
3. `finish!`, called once after the last time step

And we'll go through them one by one.

```@example callbacks
function SpeedyWeather.initialize!(
    callback::StormChaser,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::ModelSetup,
)
    # allocate recorder: number of time steps (incl initial conditions) in simulation  
    callback.maximum_surface_wind_speed = zeros(progn.clock.n_timesteps + 1)
    
    # where surface (=lowermost model layer) u,v on the grid are stored
    (;u_grid, v_grid) = diagn.layers[diagn.nlev].grid_variables
    
    # maximum wind speed of initial conditions
    callback.maximum_surface_wind_speed[1] = max_2norm(u_grid,v_grid)
    
    # (re)set counter to 1
    callback.timestep_counter = 1
end
```
The `initialize!` function has to be extended for the new callback `::StormChaser` as first
argument, then followed by prognostic and diagnostic variables and model. For correct
multiple dispatch it is important to restrict the first argument to the new `StormChaser` type
(to not call another callback instead), but the other type declarations are for clarity only.
`initialize!(::AbstractCallback,args...)` is called once just before the main time loop,
meaning after the initial conditions are set and after all other components are initialized.
We replace the vector inside our storm chaser with a vector of the correct length so that
we have a "recorder" allocated, a vector that can store the maximum surface wind speed on
every time step. We then also compute that maximum for the initial conditions and set the
time step counter to 1. We define the `max_2norm` function as follows

```@example callbacks
"""Maximum of the 2-norm for two arrays."""
function max_2norm(u::AbstractArray{T},v::AbstractArray{T}) where T
    max_norm = zero(T)      # = u² + v²
    for ij in eachindex(u, v)
        # find largest wind speed squared
        max_norm = max(max_norm, u[ij]^2 + v[ij]^2)
    end
    return sqrt(max_norm)   # take sqrt only once
end
```

Then we need to extend the `callback!` function as follows

```@example callbacks
function SpeedyWeather.callback!(
    callback::StormChaser,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::ModelSetup,
)

    # increase counter
    callback.timestep_counter += 1  
    i = callback.timestep_counter

    # where surface (=lowermost model layer) u,v on the grid are stored
    (;u_grid, v_grid) = diagn.layers[diagn.nlev].grid_variables

    # maximum wind speed at current time step
    callback.maximum_surface_wind_speed[i] = max_2norm(u_grid,v_grid)
end
```
The function signature for `callback!` is the same as for `initialize!`. You may
access anything from `progn`, `diagn` or `model`, although for a purely diagnostic
callback this  be read-only. While you could change other model components like the
land sea mask in `model.land_sea_mask` or orography etc. then you interfere with the
simulation which is more advanced and will not be discussed here.

Lastly, we extend the `finish!` function which is called once after the last time step.
This could be used, for example, to save the `maximum_surface_wind_speed` vector to
file or in case you want to find the highest wind speed across all time steps.
But in many cases you may not need to do anything, in which case you just just let
it return `nothing`.

```@example callbacks
SpeedyWeather.finish!(::StormChaser,args...) = nothing
```

!!! note "Always extend `initialize!`, `callback!` and `finish!`"
    For a custom callback you need to extend all three, `initialize!`, `callback!` and `finish!`,
    even if your callback doesn't need it. Just return `nothing` in that case. Otherwise a
    `MethodError` will occur. While we could have defined all callbacks by default to do nothing
    on each of these, this may give you the false impression that your callback is already defined
    correctly, although it's not.

## Adding a callback

Every model has a field `callbacks::AbstractVector{<:AbstractCallback}` such that the `callbacks`
keyword can be used to create a model with a vector of callbacks

```@example callbacks
spectral_grid = SpectralGrid()
dummy_callback = [NoCallback()]    # 1-element vector with dummy NoCallback only
model = PrimitiveWetModel(;spectral_grid, callbacks=dummy_callback)
model.callbacks
```

but, maybe more conveniently, a callback can be added after model construction too

```@example callbacks
storm_chaser = StormChaser(spectral_grid)
record_surface_temperature = GlobalSurfaceTemperatureCallback(spectral_grid)
append!(model.callbacks, storm_chaser)
append!(model.callbacks, record_surface_temperature)
```

which means that now in the calls to `callback!` first the dummy `NoCallback` is called
and then our storm chaser callback and then the `GlobalSurfaceTemperatureCallback` which
records the global mean surface temperature on every time step. From normal [NetCDF output](@ref)
the information these callbacks analyse would not be available,
only at the frequency of the model output, which for every time step would create way more data
and considerably slow down the simulation. Let's run the simulation and check the callbacks

```@example callbacks
simulation = initialize!(model)
run!(simulation,period=Day(3))
v = model.callbacks[2].maximum_surface_wind_speed
maximum(v)      # highest surface wind speeds in simulation [m/s]
```
The second callback is our `storm_chaser::StormChaser` (remember the first callback was a
dummy `NoCallback`), the third is the `GlobalSurfaceTemperatureCallback` with
the field `temp` is a vector of the global mean surface temperature on every
time step while the model ran for 3 days.

```@example callbacks
model.callbacks[3].temp
```
