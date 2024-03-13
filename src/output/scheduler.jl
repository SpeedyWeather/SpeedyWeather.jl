abstract type AbstractScheduler <: AbstractModelComponent end

export Scheduler
Base.@kwdef mutable struct Scheduler <: AbstractScheduler
    every::Second
    at_times::Vector{DateTime} = zeros(DateTime,0)
    schedule::BitVector = BitVector(undef,0)
    steps::Int = length(schedule)
end

export PeriodicScheduler, EventScheduler
PeriodicScheduler(every::Dates.TimePeriod = Day(1)) = Scheduler(; every)
EventScheduler(times::DateTime...) = Scheduler(every=Second(-1), at_times=DateTime[times...])
EventScheduler(times::Vector{DateTime}) = Scheduler(every=Second(-1), at_times=times)

function initialize!(scheduler::Scheduler, clock::Clock)
    is_periodic = scheduler.every.value > 0
    is_event = length(scheduler.at_times) > 0

    schedule = falses(clock.n_timesteps)    # initialize schedule as BitVector

    if is_periodic
        every_n_timesteps = max(1,round(Int, Millisecond(scheduler.every).value/clock.Δt_millisec.value))
        schedule[1:every_n_timesteps:end] .= true
    elseif is_event
        time = clock.start
        for event in scheduler.at_times
            i::Int = 1
            event_scheduled::Bool = false
            while ~event_scheduled && i <= clock.n_timesteps
                event_scheduled = time >= event >= clock.start
                schedule[i] = event_scheduled
                time += clock.Δt_millisec
                i += 1
            end
        end
    else
        @warn "Empty schedule for Scheduler."
    end

    scheduler.schedule = schedule   # set schedule
    scheduler.steps = sum(schedule)
    return scheduler
end

isscheduled(S::Scheduler, clock::Clock) = S.schedule[clock.timestep_counter]