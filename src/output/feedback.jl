abstract type AbstractFeedback end

export Feedback

"""
Feedback struct that contains options and object for command-line feedback
like the progress meter.
$(TYPEDFIELDS)"""
@kwdef mutable struct Feedback <: AbstractFeedback
    "print feedback to REPL?, default is isinteractive(), true in interactive REPL mode"
    verbose::Bool = isinteractive()

    "check for NaRs in the prognostic variables"
    debug::Bool = true
    
    "write a progress.txt file? State synced with NetCDFOutput.output"
    output::Bool = false

    "identification of run, taken from ::OutputWriter"
    id::String = ""

    "path to run folder, taken from ::OutputWriter"
    run_path::String = ""

    # PROGRESS
    "struct containing everything progress related"
    progress_meter::ProgressMeter.Progress =
        ProgressMeter.Progress(1, enabled=isinteractive, showspeed=true; desc="Weather is speedy: ")

    "txt is a Nothing in case of no output"
    progress_txt::Union{IOStream, Nothing} = nothing

    "did Infs/NaNs occur in the simulation?"
    nars_detected::Bool = false
end

function Base.show(io::IO, F::AbstractFeedback)
    println(io, "$(typeof(F)) <: AbstractFeedback")
    keys = propertynames(F)
    print_fields(io, F, keys)
end

function Base.show(io::IO, P::ProgressMeter.Progress)
    println(io, "$(typeof(P)) <: ProgressMeter.AbstractProgress")
    keys = propertynames(P)
    print_fields(io, P, keys)
end

"""
$(TYPEDSIGNATURES)
Initializes the a `Feedback` struct."""
function initialize!(feedback::Feedback, clock::Clock, model::AbstractModel)

    # set to false to recheck for NaRs
    feedback.nars_detected = false

    # hack: redefine element in global constant dt_in_sec
    # used to pass on the time step to ProgressMeter.speedstring
    DT_IN_SEC[] = model.time_stepping.Δt_sec

    if feedback.output   # with netcdf output write progress.txt
        (; run_path, id) = feedback
        SG = model.spectral_grid
        L = model.time_stepping
        days = Second(clock.period).value/(3600*24)
        
        # create progress.txt file in run_????/
        progress_txt = open(joinpath(run_path, "progress.txt"), "w")
        s = "Starting SpeedyWeather.jl run $id on "*
                Dates.format(Dates.now(), Dates.RFC1123Format)
        write(progress_txt, s*"\n")
        write(progress_txt, "Integrating:\n")
        write(progress_txt, "$SG\n")
        write(progress_txt, "Time: $days days at Δt = $(L.Δt_sec)s\n")
        write(progress_txt, "\nAll data will be stored in $run_path\n")
        feedback.progress_txt = progress_txt
    end

    # reinitalize progress meter, minus one to exclude first_timesteps! which contain compilation
    # only do now for benchmark accuracy
    (; showspeed, desc) = feedback.progress_meter
    (; verbose) = feedback
    feedback.progress_meter = ProgressMeter.Progress(clock.n_timesteps-1; enabled=verbose, showspeed, desc)
end

"""
$(TYPEDSIGNATURES)
Calls the progress meter and writes every 5% progress increase to txt."""
function progress!(feedback::Feedback)

    # update progress meter and unpack counter after update
    ProgressMeter.next!(feedback.progress_meter)    
    (; counter, n ) = feedback.progress_meter

    # write progress to txt file too
    if (counter/n*100 % 1) > ((counter+1)/n*100 % 1)  
        percent = round(Int, (counter+1)/n*100)      # % of time steps completed
        if feedback.output && (percent % 5 == 0)    # write every 5% step in txt 
            write(feedback.progress_txt, @sprintf("\n%3d%%", percent))
            r = remaining_time(feedback.progress_meter)
            write(feedback.progress_txt, ", ETA: $r")

            time_elapsed = feedback.progress_meter.tlast - feedback.progress_meter.tinit
            s = speedstring(time_elapsed/counter, DT_IN_SEC[])
            write(feedback.progress_txt, ", $s")

            feedback.nars_detected && write(feedback.progress_txt, ", NaRs detected.")
            flush(feedback.progress_txt)
        end
    end
end

function progress!( feedback::Feedback,
                    progn::PrognosticVariables)
    progress!(feedback)
    feedback.debug && nar_detection!(feedback, progn)
end

"""
$(TYPEDSIGNATURES)
Finalises the progress meter and the progress txt file."""
function finish!(F::Feedback)
    ProgressMeter.finish!(F.progress_meter)
    
    if F.output     # write final progress to txt file
        time_elapsed = F.progress_meter.tlast - F.progress_meter.tinit
        s = "\nIntegration done in $(readable_secs(time_elapsed))."
        write(F.progress_txt, "\n$s\n")
        flush(F.progress_txt)
    end
end

"""
$(TYPEDSIGNATURES)
Detect NaR (Not-a-Real) in the prognostic variables."""
function nar_detection!(feedback::Feedback, progn::PrognosticVariables)

    feedback.nars_detected && return nothing    # escape immediately if nans already detected
    i = feedback.progress_meter.counter         # time step
    vor0 = progn.vor[2][1, end]                 # only check 0-0 mode of surface vorticity

    # just check first harmonic, spectral transform propagates NaRs globally anyway
    nars_detected_here = ~isfinite(vor0)
    nars_detected_here && @warn "NaN or Inf detected at time step $i"
    feedback.nars_detected = nars_detected_here
end

"""
$(TYPEDSIGNATURES)
Estimates the remaining time from a `ProgresssMeter.Progress`.
Adapted from ProgressMeter.jl"""
function remaining_time(p::ProgressMeter.Progress)
    elapsed_time = time() - p.tinit
    est_total_time = elapsed_time * (p.n - p.start) / (p.counter - p.start)
    if 0 <= est_total_time <= typemax(Int)
        eta_sec = round(Int, est_total_time - elapsed_time )
        eta = ProgressMeter.durationstring(eta_sec)
    else
        eta = "N/A"
    end
    return eta
end

"""
$(TYPEDSIGNATURES)
Define a ProgressMeter.speedstring method that also takes a time step
`dt_in_sec` to translate sec/iteration to days/days-like speeds."""
function speedstring(sec_per_iter, dt_in_sec)
    if sec_per_iter == Inf
        return "  N/A  days/day"
    end

    sim_time_per_time = dt_in_sec/sec_per_iter

    for (divideby, unit) in (   (365*1_000, "millenia"),
                                (365, "years"),
                                (1, "days"),
                                (1/24, "hours"))    
        if (sim_time_per_time / divideby) > 2
            return @sprintf "%5.2f %2s/day" (sim_time_per_time / divideby) unit
        end
    end
    return " <2 hours/days"
end

# hack: define global constant whose element will be changed in initialize_feedback
# used to pass on the time step to ProgressMeter.speedstring via calling this
# constant from the ProgressMeter module
const DT_IN_SEC = Ref(1800.0)

# "extend" the speedstring function from ProgressMeter by defining it for ::AbstractFloat
# not just ::Any to effectively overwrite it
function ProgressMeter.speedstring(sec_per_iter::AbstractFloat)
    dt_in_sec = SpeedyWeather.DT_IN_SEC[]   # pull global "constant"
    speedstring(sec_per_iter, dt_in_sec)
end

"""
$(TYPEDSIGNATURES)
Returns `Dates.CompoundPeriod` rounding to either (days, hours), (hours, minutes), (minutes,
seconds), or seconds with 1 decimal place accuracy for >10s and two for less.
E.g.
```@example
julia> using SpeedyWeather: readable_secs

julia> readable_secs(12345)
```
"""
function readable_secs(secs::Real)
    millisecs = Dates.Millisecond(round(secs * 1000))
    if millisecs >= Dates.Day(1)
        return Dates.canonicalize(round(millisecs, Dates.Hour))
    elseif millisecs >= Dates.Hour(1)
        return Dates.canonicalize(round(millisecs, Dates.Minute))
    elseif millisecs >= Dates.Minute(1)
        return Dates.canonicalize(round(millisecs, Dates.Second))
    elseif millisecs >= Dates.Second(10)
        return Dates.canonicalize(round(millisecs, Dates.Millisecond(100)))
    end
    return Dates.canonicalize(round(millisecs, Dates.Millisecond(10)))
end