abstract type AbstractSchedule <: AbstractModelComponent end

export Schedule
Base.@kwdef mutable struct Schedule <: AbstractSchedule
    "Execute every time period, first timestep included. Default=never."
    every::Second = Second(typemax(Int))
    
    "Events scheduled at times"
    times::Vector{DateTime} = zeros(DateTime,0)

    "Actual schedule, true=execute this timestep, false=don't"
    schedule::BitVector = BitVector(undef,0)

    "Number of scheduled executions"
    steps::Int = length(schedule)

    "Count up the number of executions"
    counter::Int = 0
end

Schedule(times::DateTime...) = Schedule(times=DateTime[times...])

function initialize!(scheduler::Schedule, clock::Clock)
    schedule = falses(clock.n_timesteps)    # initialize schedule as BitVector

    # PERIODIC SCHEDULE, always AFTER scheduler.every time period has passed
    if scheduler.every.value < typemax(Int)
        every_n_timesteps = max(1,round(Int, scheduler.every/clock.Δt))
        schedule[every_n_timesteps:every_n_timesteps:end] .= true

        prev_every = scheduler.every
        now_every = convert(typeof(prev_every), every_n_timesteps*clock.Δt)
        scheduler.every = now_every
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

function isscheduled(S::Schedule, clock::Clock)
    is_scheduled = S.schedule[clock.timestep_counter]
    S.counter += is_scheduled
    return is_scheduled
end
