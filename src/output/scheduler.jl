abstract type AbstractScheduler <: AbstractModelComponent end

export Scheduler
Base.@kwdef mutable struct Scheduler <: AbstractScheduler
    every::Second = Second(typemax(Int))
    times::Vector{DateTime} = zeros(DateTime,0)
    schedule::BitVector = BitVector(undef,0)
    steps::Int = length(schedule)
end

export PeriodicScheduler, EventScheduler
PeriodicScheduler(every::Dates.TimePeriod) = Scheduler(; every)
EventScheduler(times::DateTime...) = Scheduler(times=DateTime[times...])
EventScheduler(times::Vector{DateTime}) = Scheduler(; times)

function initialize!(scheduler::Scheduler, clock::Clock)
    schedule = falses(clock.n_timesteps)    # initialize schedule as BitVector

    # periodic schedule
    every_n_timesteps = max(1,round(Int, Millisecond(scheduler.every).value/clock.Δt_millisec.value))
    schedule[1:every_n_timesteps:end] .= true

    # event schedule
    time = clock.start
    for event in scheduler.at_times
        i::Int = 1
        event_scheduled::Bool = false
        while ~event_scheduled && i <= clock.n_timesteps
            event_scheduled = time >= event >= clock.start
            schedule[i] |= event_scheduled
            time += clock.Δt_millisec
            i += 1
        end
    end

    scheduler.schedule = schedule   # set schedule
    scheduler.steps = sum(schedule)
    scheduler.steps == 0 && @warn "Empty schedule for Scheduler."

    return scheduler
end

isscheduled(S::Scheduler, clock::Clock) = S.schedule[clock.timestep_counter]