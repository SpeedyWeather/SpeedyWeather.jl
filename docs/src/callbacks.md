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
    model::AbstractModel,
)
    # allocate recorder: number of time steps (incl initial conditions) in simulation  
    callback.maximum_surface_wind_speed = zeros(progn.clock.n_timesteps + 1)
    
    # where surface (=lowermost model layer) u, v on the grid are stored
    u_grid = diagn.grid.u_grid[:, diagn.nlayers]
    v_grid = diagn.grid.u_grid[:, diagn.nlayers]

    # maximum wind speed of initial conditions
    callback.maximum_surface_wind_speed[1] = max_2norm(u_grid, v_grid)
    
    # (re)set counter to 1
    callback.timestep_counter = 1
end
```
The `initialize!` function has to be extended for the new callback `::StormChaser` as first
argument, then followed by prognostic and diagnostic variables and model. For correct
multiple dispatch it is important to restrict the first argument to the new `StormChaser` type
(to not call another callback instead), but the other type declarations are for clarity only.
`initialize!(::AbstractCallback, args...)` is called once just before the main time loop,
meaning after the initial conditions are set and after all other components are initialized.
We replace the vector inside our storm chaser with a vector of the correct length so that
we have a "recorder" allocated, a vector that can store the maximum surface wind speed on
every time step. We then also compute that maximum for the initial conditions and set the
time step counter to 1. We define the `max_2norm` function as follows

```@example callbacks
"""Maximum of the 2-norm of elements across two arrays."""
function max_2norm(u::AbstractArray{T}, v::AbstractArray{T}) where T
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
    model::AbstractModel,
)

    # increase counter
    callback.timestep_counter += 1  
    i = callback.timestep_counter

    # where surface (=lowermost model layer) u, v on the grid are stored
    u_grid = diagn.grid.u_grid[:, diagn.nlayers]
    v_grid = diagn.grid.u_grid[:, diagn.nlayers]

    # maximum wind speed at current time step
    callback.maximum_surface_wind_speed[i] = max_2norm(u_grid, v_grid)
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
SpeedyWeather.finish!(::StormChaser, args...) = nothing
```

!!! note "Always extend `initialize!`, `callback!` and `finish!`"
    For a custom callback you need to extend all three, `initialize!`, `callback!` and `finish!`,
    even if your callback doesn't need it. Just return `nothing` in that case. Otherwise a
    `MethodError` will occur. While we could have defined all callbacks by default to do nothing
    on each of these, this may give you the false impression that your callback is already defined
    correctly, although it's not.

## Adding a callback

Every model has a field `callbacks::Dict{Symbol, AbstractCallback}` such that the `callbacks`
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
model = PrimitiveWetModel(; spectral_grid, callbacks)
add!(model.callbacks, :callback_added_afterwards => NoCallback())
add!(model, :callback_added_afterwards2 => NoCallback())
```

Note how the first argument can be `model.callbacks` as outlined in the sections above
because this is the callbacks dictionary, but also simply
`model`, which will add the callback to `model.callbacks`. It's equivalent.
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

## [Intrusive callbacks](@id intrusive_callbacks)

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

As an example of a callback that changes the model components see

- Millenium flood: [Time-dependent land-sea mask](@ref)

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

## Schedules

For convenience, SpeedyWeather.jl implements a [`Schedule`](@ref) which helps to
schedule when callbacks are called. Because in many situations you don't want to call
them on every time step but only periodically, say once a day, or only on specific
dates and times, e.g. Jan 1 at noon. Several examples how to create schedules

```@example schedule
using SpeedyWeather

# execute on timestep at or after Jan 2 2000
event_schedule = Schedule(DateTime(2000,1,2))   

# several events scheduled
events = (DateTime(2000,1,3), DateTime(2000,1,5,12))
several_events_schedule = Schedule(events...)

# provided as Vector{DateTime} with times= keyword
always_at_noon = [DateTime(2000,1,i,12) for i in 1:10]
noon_schedule = Schedule(times=always_at_noon)

# or using every= for periodic execution, here once a day
periodic_schedule = Schedule(every=Day(1))
```

A `Schedule` has 5 fields, see [`Schedule`](@ref). `every` is an option
to create a periodic schedule to execute every time that indicated
period has passed. `steps` and `counter` will let you know how many
callback execution steps there are and count them up. `times` is a 
`Vector{DateTime}` containing scheduled events. `schedule` is the
actual schedule inside a `Schedule`, implemented as `BitVector`
indicating whether to execute on a given time step (`true`) or not
(`false`).

Let's show how to use a `Schedule` inside a callback

```@example schedule
struct MyScheduledCallback <: SpeedyWeather.AbstractCallback
    schedule::Schedule
    # add other fields here that you need
end

function SpeedyWeather.initialize!(
    callback::MyScheduledCallback,
    progn::PrognosticVariables,
    args...
)
    # when initializing a scheduled callback also initialize its schedule!
    initialize!(callback.schedule, progn.clock)

    # initialize other things in your callback here
end

function SpeedyWeather.callback!(
    callback::MyScheduledCallback,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    # scheduled callbacks start with this line to execute only when scheduled!
    # else escape immediately
    isscheduled(callback.schedule, progn.clock) || return nothing

    # Just print the North Pole surface temperature to screen
    (;time) = progn.clock
    temp_at_north_pole = diagn.grid.temp_grid[1,end]

    @info "North pole has a temperature of $temp_at_north_pole on $time."
end

# nothing needs to be done when finishing
SpeedyWeather.finish!(::MyScheduledCallback, args...) = nothing
```

So in summary
- add a field `schedule::Schedule` to your callback
- add the line `initialize!(callback.schedule, progn.clock)` when initializing your callback
- start your `callback!` method with `isscheduled(callback.schedule, progn.clock) || return nothing` to execute only when scheduled

A `Schedule` is a field inside a callback as this allows you the set the callbacks
desired schedule when creating it. In the example above we can create our callback
that is supposed to print the North Pole's temperature like so
```@example schedule
north_pole_temp_at_noon_jan9 = MyScheduledCallback(Schedule(DateTime(2000,1,9,12)))
```
The default for `every` is `typemax(Int)` indicating "never". This just means that there
is no periodically reoccuring schedule, only `schedule.times` would include some times
for events that are scheduled. Now let's create a primitive equation model with that callback

```@example schedule
spectral_grid = SpectralGrid(trunc=31, nlayers=5)
model = PrimitiveWetModel(;spectral_grid)
add!(model.callbacks, north_pole_temp_at_noon_jan9)

# start simulation 7 days earlier
simulation = initialize!(model, time = DateTime(2000,1,2,12))
run!(simulation, period=Day(10))
nothing # hide
```

So the callback gives us the temperature at the North Pole exactly when scheduled.
We could have also stored this temperature, or conditionally changed parameters inside
the model. There are plenty of ways how to use the scheduling, either by event, or in contrast,
we could also schedule for once a day. As illustrated in the following

```@example schedule
north_pole_temp_daily = MyScheduledCallback(Schedule(every=Day(1)))
add!(model.callbacks, north_pole_temp_daily)

# resume simulation, start time is now 2000-1-12 noon
run!(simulation, period=Day(5))
nothing # hide
```

Note that the previous callback is still part of the model, we haven't
deleted it with `delete!`. But because it's scheduled for a specific
time that's in the past now that we resume the simulation it's schedule
is empty (which is thrown as a warning). However, our new callback,
scheduled daily, is active and prints daily at noon, because the
simulation start time was noon.

### Scheduling logic

An event `Schedule` (created with `DateTime` object(s)) for callbacks, executes 
on or after the specified times.
For two consecutive time steps ``i``, ``i+1``, an event is scheduled at ``i+1``
when it occurs in ``(i,i+1]``. So a simulation with timestep `i` on Jan-1 at 1am,
and ``i+1`` at 2am, will execute a callback scheduled for 1am at ``i`` but scheduled
for 1am and 1s (=01:00:01 on a 24H clock) at 2am. Because callbacks are
always executed _after_ a timestep this also means that a simulation starting
at midnight with a callback scheduled for midnight will not execute this
callback as it is outside of the ``(i, i+1]`` range. You'd need to include
this execution into the initialization. If several events inside the
`Schedule` fall into the same time step (in the example above, 1am and 1s and
1am 30min) the execution will not happen twice. Think of a scheduled callback
as a binary "should the callback be executed now or not?". Which is in fact
how it's implemented, as a `BitVector` of the length of the number of time steps.
If the bit at a given timestep is true, execute, otherwise not.

A periodic `Schedule` (created with `every = Hour(2)` or similar) will execute
on the timestep _after_ that period (here 2 hours) has passed. If a simulation
starts at midnight with one hour time steps then execution would take place
after the timestep from 1am to 2am because that's when the clock switches to
2am which is 2 hours after the start of the simulation. Note that therefore
the initial timestep is not included, however, the last time step would be
if the period is a multiple of the scheduling period. If the first timestep
should be included (e.g. you want to do something with the initial conditions)
then you'll need to include that into the initialization of the callback.

Periodic schedules which do not match the simulation time step will be adjusted
by rounding. Example, if you want a schedule which executes every hour
but your simulation time step is 25min then it will be adjusted to execute
every 2nd time step, meaning every 50min and not 1 hour. However, an info
will be thrown if that is the case

```@example schedule
odd_schedule = MyScheduledCallback(Schedule(every = Minute(70)))
add!(model.callbacks, odd_schedule)

# resume simulation for 4 hours
run!(simulation, period=Hour(4))
nothing # hide
```

Now we get two empty schedules, one from callback that's supposed to
execute on Jan 9 noon (this time has passed in our simulation) and
one from the daily callback (we're not simulating for a day).
You could just `delete!` those callbacks. You can see that while we
wanted our `odd_schedule` to execute every 70min, it has to adjust it
to every 60min to match the simulation time step of 30min.

After the model initialization you can always check the simulation time step
from `model.time_stepping` 

```@example schedule
model.time_stepping
```

Or converted into minutes (the time step internally is at millisecond accuracy)

```@example schedule
Minute(model.time_stepping.Δt_millisec)
```

which illustrates why the adjustment of our callback frequency was necessary.