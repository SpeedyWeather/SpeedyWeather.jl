abstract type AbstractSchedule <: AbstractModelComponent end

# otherwise one needs to write SpeedyWeather.isscheduled inside custom callbacks or components
export isscheduled

export Schedule

"""
A schedule for callbacks, to execute them periodically after a given period has
passed (first timestep excluded) or on/after specific times (initial time excluded).
For two consecutive time steps i, i+1, an event is scheduled at i+1 when it occurs
in (i,i+1]. Similarly, a periodic schedule of period p will be executed at
start+p, but p is rounded to match the multiple of the model timestep.
Periodic schedules and event schedules can be combined, executing at both.
A Schedule is supposed to be added into callbacks as fields

    @kwdef struct MyCallback
        schedule::Schedule = Schedule(every=Day(1))
        other_fields
    end

see also `initialize!(::Schedule, ::Clock)` and `isscheduled(::Schedule, ::Clock)`.

Fields
$(TYPEDFIELDS)"""
@kwdef mutable struct Schedule <: AbstractSchedule
    "[OPTION] Execute every time period, first timestep excluded. Default=never."
    every::Second = Second(typemax(Int))

    "[OPTION] Events scheduled at times"
    times::Vector{DateTime} = zeros(DateTime, 0)

    "Actual schedule, true=execute this timestep, false=don't"
    schedule::BitVector = BitVector(undef, 0)

    "Number of scheduled executions"
    steps::Int = length(schedule)

    "Count up the number of executions"
    counter::Int = 0
end

"""$(TYPEDSIGNATURES)
A Schedule based on DateTime arguments, For two consecutive time steps i, i+1, an event is
scheduled at i+1 when it occurs in (i,i+1]."""
Schedule(times::DateTime...) = Schedule(times = DateTime[times...])

"""$(TYPEDSIGNATURES)
Initialize a Schedule with a Clock (which is assumed to have been initialized).
Takes both scheduler.every and scheduler.times into account, such that
both periodic and events can be scheduled simultaneously. But execution will
happen only once if they coincide on a given time step."""
function initialize!(scheduler::Schedule, clock::Clock)
    clock.Δt == Millisecond(0) && error("Clock needs to be initialized with a non-zero time step before initializing schedule.")
    schedule = falses(clock.n_time_steps)    # initialize schedule as BitVector

    # PERIODIC SCHEDULE, always AFTER scheduler.every time period has passed
    if scheduler.every.value < typemax(Int)
        every_n_time_steps = max(1, round(Int, scheduler.every / clock.Δt))
        schedule[every_n_time_steps:every_n_time_steps:end] .= true

        prev_every = readable_secs(scheduler.every.value)
        scheduler.every = Second(second(every_n_time_steps * clock.Δt))
        now_every = readable_secs(scheduler.every.value)
        s = "Scheduler adjusted from every $prev_every to every $now_every to match timestep."
        now_every != prev_every && @info s
    end

    # ADD EVENT SCHEDULE, on first timestep ON or AFTER the scheduled time
    # event on clock.start will not be executed, ()
    for event in scheduler.times
        i = ceil(Int, (event - clock.start) / clock.Δt)   #
        if 0 < i <= clock.n_time_steps  # event needs to take place in (start, end]
            schedule[i] = true          # add to schedule,
        end
    end

    scheduler.schedule = schedule   # set schedule
    scheduler.steps = sum(schedule) # number of executions scheduled periodic + events
    scheduler.steps == 0 && @warn "Empty schedule."
    scheduler.counter = 0           # reset counter

    return scheduler
end

"""
$(TYPEDSIGNATURES)
Evaluate whether (e.g. a callback) should be scheduled at the timestep given
in clock. Returns true for scheduled executions, false for no execution on
this time step."""
function isscheduled(S::Schedule, clock::Clock)
    clock.time_step_counter == 0 && return false
    is_scheduled = S.schedule[clock.time_step_counter]
    S.counter += is_scheduled
    return is_scheduled
end

export PeriodicSchedule

"""
A periodic schedule, to execute something every period `every` (first time step excluded),
rounded to a multiple of the time step. Unlike `Schedule` it does not need the clock (start
time or length of the simulation) and can therefore be initialized with the time stepper only,
e.g. in `initialize!(model)`, converting `every` into `every_n_steps`. Use as

    @kwdef struct MyComponent
        schedule::PeriodicSchedule = PeriodicSchedule(every = Day(1))
    end

with `isscheduled(component.schedule, clock) || return nothing`.
Fields
$(TYPEDFIELDS)"""
@kwdef mutable struct PeriodicSchedule <: AbstractSchedule
    "[OPTION] Execute every time period, rounded to a multiple of the time step, `Second(0)` = never"
    every::Second = Day(1)

    "[DERIVED] `every` in time steps (at least 1), 0 = never"
    every_n_steps::Int = 0

    "[DERIVED] Count up the number of executions"
    counter::Int = 0
end

"""$(TYPEDSIGNATURES)
A `PeriodicSchedule` executing `every` period."""
PeriodicSchedule(every::Period) = PeriodicSchedule(; every = Second(every))

"""$(TYPEDSIGNATURES)
Initialize a `PeriodicSchedule` with the time step of `time_stepping`: precomputes `every_n_steps`
from `schedule.every` (at least 1, 0 = never for `every = Second(0)`) and resets the counter."""
function initialize!(schedule::PeriodicSchedule, time_stepping::AbstractTimeStepper)
    (; every) = schedule
    Δt = time_stepping.Δt_millisec
    schedule.every_n_steps = every > Second(0) ? max(1, round(Int, Millisecond(every) / Δt)) : 0
    schedule.counter = 0
    return schedule
end

"""$(TYPEDSIGNATURES)
Whether `schedule` executes on the time step of `clock`, i.e. every `every_n_steps` time steps
(spin-up steps excluded). Counts the executions."""
function isscheduled(schedule::PeriodicSchedule, clock::Clock)
    is_scheduled = isscheduled(schedule.every_n_steps, clock)
    count_execution!(schedule, is_scheduled)
    return is_scheduled
end

# count the executions of a schedule, see SpeedyWeatherReactantExt for a traced `is_scheduled`
count_execution!(schedule::AbstractSchedule, is_scheduled::Bool) = (schedule.counter += is_scheduled; nothing)

export StepSchedule

"""
An immutable schedule executing every `every_n_steps` time steps (first time step and spin-up
steps excluded), 0 = never, 1 = every time step. Needs no initialization. Use as
`isscheduled(component.schedule, clock) || return nothing`.
Fields
$(TYPEDFIELDS)"""
@kwdef struct StepSchedule <: AbstractSchedule
    "[OPTION] Execute every n time steps, 0 = never"
    every_n_steps::Int = 1
end

"""$(TYPEDSIGNATURES)
A `StepSchedule` executing every `every_n_steps` time steps."""
StepSchedule(every_n_steps::Integer) = StepSchedule(; every_n_steps)

"""$(TYPEDSIGNATURES)
Nothing to initialize for a `StepSchedule`."""
initialize!(schedule::StepSchedule, ::AbstractTimeStepper) = schedule

"""$(TYPEDSIGNATURES)
Whether `schedule` executes on the time step of `clock`, i.e. every `every_n_steps` time steps
(spin-up steps excluded)."""
isscheduled(schedule::StepSchedule, clock::Clock) = isscheduled(schedule.every_n_steps, clock)

# every n time steps of size Δt (spin-up steps excluded), 0 = never
function isscheduled(every_n_steps::Integer, clock::Clock)
    every_n_steps > 0 || return false
    n = clock.time_step_counter
    return (n > 0) & (n % every_n_steps == 0)
end
