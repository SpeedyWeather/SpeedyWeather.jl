# Callbacks

SpeedyWeather.jl implements a callback system to let users include a flexible piece of code
into the time stepping. You can think about the main time loop *calling back* to check whether
anything else should be done before continuing with the next time step. The callback system
here is called *after* the time step only (plus one call at `initialize!` and one at `finish!`),
we currently do not implement other callsites.

Callbacks are mainly introduced for diagnostic purposes, meaning that they do not influence
the simulation, and access the prognostic variables and the model components in a read-only
fashion. However, a callback is not strictly prevented from changing prognostic or diagnostic
variables or the model. For example, you may define a callback that changes the orography
during the simulation. In general, one has to keep the general order of executions during a
time step in mind (valid for all models)

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

```@example callbacks
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
"""Maximum of the 2-norm of elements across two arrays."""
function max_2norm(u::AbstractArray{T},v::AbstractArray{T}) where T
    max_norm = zero(T)      # = u² + v²
    for ij in eachindex(u, v)
        # find largest wind speed squared
        max_norm = max(max_norm, u[ij]^2 + v[ij]^2)
    end
    return sqrt(max_norm)   # take sqrt only once
end
```

Note that this function is defined in the scope `Main` and not inside SpeedyWeather, this is absolutely
possible due to Julia's scope of variables which will use `max_2norm` from `Main` scope if it doesn't
exist in the global scope inside the `SpeedyWeather` module scope.
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
callback this should be read-only. While you could change other model components like the
land sea mask in `model.land_sea_mask` or orography etc. then you interfere with the
simulation which is more advanced and will be discussed in Intrusive callbacks below.

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

## Adding a callback

Every model has a field `callbacks::Dict{Symbol,AbstractCallback}` such that the `callbacks`
keyword can be used to create a model with a dictionary of callbacks. Callbacks are identified
with a `Symbol` key inside such a dictionary. We have a convenient `CallbackDict` generator function
which can be used like `Dict` but the key-value pairs have to be of type `Symbol`-`AbstractCallback`.
Let us illustrate this with the dummy callback `NoCallback` (which is a callback that returns `nothing`
on `initialize!`, `callback!` and `finish!`)

```@example callbacks
callbacks = CallbackDict()                                  # empty dictionary
callbacks = CallbackDict(:my_callback => NoCallback())      # key => callback
```
If you don't provide a key a random key will be assigned
```@example callbacks
callbacks = CallbackDict(NoCallback())
```
and you can add (or delete) additional callbacks
```@example callbacks
add!(callbacks, NoCallback())                   # this will also pick a random key
add!(callbacks, :my_callback => NoCallback())   # use key :my_callback
delete!(callbacks, :my_callback)                # remove by key
callbacks
```
And you can chain them too
```@example callbacks
add!(callbacks, NoCallback(), NoCallback())                     # random keys
add!(callbacks, :key1 => NoCallback(), :key2 => NoCallback())   # keys provided
```
Meaning that callbacks can be added before and after model construction


```@example callbacks
spectral_grid = SpectralGrid()
callbacks = CallbackDict(:callback_added_before => NoCallback())
model = PrimitiveWetModel(;spectral_grid, callbacks)
add!(model.callbacks,:callback_added_afterwards => NoCallback())
```
Let us add two more meaningful callbacks

```@example callbacks
storm_chaser = StormChaser(spectral_grid)
record_surface_temperature = GlobalSurfaceTemperatureCallback(spectral_grid)
add!(model.callbacks, :storm_chaser => storm_chaser)
add!(model.callbacks, :temperature => record_surface_temperature)
```

which means that now in the calls to `callback!` first the two dummy `NoCallback`s are called
and then our storm chaser callback and then the `GlobalSurfaceTemperatureCallback` which
records the global mean surface temperature on every time step. From normal [NetCDF output](@ref)
the information these callbacks analyse would not be available,
only at the frequency of the model output, which for every time step would create way more data
and considerably slow down the simulation. Let's run the simulation and check the callbacks

```@example callbacks
simulation = initialize!(model)
run!(simulation, period=Day(3))
v = model.callbacks[:storm_chaser].maximum_surface_wind_speed
maximum(v)      # highest surface wind speeds in simulation [m/s]
```
Cool, our `StormChaser` callback with the key `:storm_chaser` has been recording maximum
surface wind speeds in [m/s]. And the `:temperature` callback a time series of
global mean surface temperatures in Kelvin on every time step while the model ran for 3 days.

```@example callbacks
model.callbacks[:temperature].temp
```

## Intrusive callbacks

In the sections above, callbacks were introduced as a tool to define custom
diagnostics or simulation output. This is the simpler and recommended way of using 
them but nothing stops you from defining a callback that is *intrusive*, meaning
that it can alter the prognostic or diagnostic variables or the model.

Changing any components of the model, e.g. boundary conditions like orography
or the land-sea mask through a callback is possible although one should notice
that this only comes into effect on the next time step given the execution
order mentioned above. One could for example run a simulation for a certain
period and then start moving continents around. Note that for physical consistency
this should be reflected in the orography, land-sea mask, as well as in the available
sea and land-surface temperatures, but one is free to do this only partially too.
Another example would be to switch on/off certain model components over time.
If these components are implemented as *mutable* struct then one could define
a callback that weakens their respective strength parameter over time.

Changing the diagnostic variables, however, will not have any effect. All of
them are treated as work arrays, meaning that their state is completely
overwritten on every time step.  Changing the prognostic variables in spectral space
directly is not advised though possible because this can easily lead to stability issues.
It is generally easier to implement something like this as a parameterization, forcing or
drag term (which can also be made time-dependent).

Overall, callbacks give the user a wide range of possibilities to diagnose 
the simulation while running or to interfere with a simulation. We therefore
encourage users to use callbacks as widely as possible, but if you run
into any issues please open an issue in the repository and explain what
you'd like to achieve and which errors you are facing. We are happy to help.
