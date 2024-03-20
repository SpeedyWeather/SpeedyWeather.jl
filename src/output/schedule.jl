abstract type AbstractSchedule <: AbstractModelComponent end

export Schedule

"""
A schedule for callbacks, to execute them periodically after a given period has
passed (first timestep excluded) or on/after specific times (initial time excluded).
For two consecutive time steps i, i+1, an event is scheduled at i+1 when it occurs
in (i,i+1]. Similarly, a periodic schedule of period p will be executed at
start+p, but p is rounded to match the multiple of the model timestep.
Periodic schedules and event schedules can be combined, executing at both.
A Schedule is supposed to be added into callbacks as fields

    Base.@kwdef struct MyCallback
        schedule::Schedule = Schedule(every=Day(1))
        other_fields
    end

see also initialize!(::Schedule,::Clock) and isscheduled(::Schedule,::Clock).
Fields
$(TYPEDFIELDS)"""
Base.@kwdef mutable struct Schedule <: AbstractSchedule
    "[OPTION] Execute every time period, first timestep excluded. Default=never."
    every::Second = Second(typemax(Int))
    
    "[OPTION] Events scheduled at times"
    times::Vector{DateTime} = zeros(DateTime,0)

    "Actual schedule, true=execute this timestep, false=don't"
    schedule::BitVector = BitVector(undef,0)

    "Number of scheduled executions"
    steps::Int = length(schedule)

    "Count up the number of executions"
    counter::Int = 0
end

"""
$(TYPEDSIGNATURES)
A Schedule based on DateTime arguments, For two consecutive time steps i, i+1, an event is
scheduled at i+1 when it occurs in (i,i+1]."""
Schedule(times::DateTime...) = Schedule(times=DateTime[times...])

"""
$(TYPEDSIGNATURES)
Initialize a Schedule with a Clock (which is assumed to have been initialized).
Takes both scheduler.every and scheduler.times into account, such that
both periodic and events can be scheduled simulataneously. But execution will
happen only once if they coincide on a given time step."""
function initialize!(scheduler::Schedule, clock::Clock)
    schedule = falses(clock.n_timesteps)    # initialize schedule as BitVector

    # PERIODIC SCHEDULE, always AFTER scheduler.every time period has passed
    if scheduler.every.value < typemax(Int)
        every_n_timesteps = max(1,round(Int, scheduler.every/clock.Δt))
        schedule[every_n_timesteps:every_n_timesteps:end] .= true

        prev_every = readable_secs(scheduler.every.value)
        scheduler.every = Second(Dates.second(every_n_timesteps*clock.Δt))
        now_every = readable_secs(scheduler.every.value)
        s = "Scheduler adjusted from every $prev_every to every $now_every to match timestep."
        now_every != prev_every && @info s
    end

    # ADD EVENT SCHEDULE, on first timestep ON or AFTER the scheduled time
    # event on clock.start will not be executed, ()
    for event in scheduler.times
        i = ceil(Int, (event - clock.start)/clock.Δt)   #
        if 0 < i <= clock.n_timesteps   # event needs to take place in (start, end] 
            schedule[i] = true          # add to schedule, 
        end
    end

    scheduler.schedule = schedule   # set schedule
    scheduler.steps = sum(schedule) # number of executions scheduled periodic + events
    scheduler.steps == 0 && @warn "Empty schedule."
    scheduler.counter = 0           # reset counter

    return scheduler
end

# otherwise one needs to write SpeedyWeather.isscheduled inside custom callbacks
export isscheduled  

"""
$(TYPEDSIGNATURES)
Evaluate whether (e.g. a callback) should be scheduled at the timestep given
in clock. Returns true for scheduled executions, false for no execution on
this time step."""
function isscheduled(S::Schedule, clock::Clock)
    is_scheduled = S.schedule[clock.timestep_counter]
    S.counter += is_scheduled
    return is_scheduled
end
