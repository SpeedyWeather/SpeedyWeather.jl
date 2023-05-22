mutable struct Feedback <: AbstractFeedback
    verbose::Bool                           # print feedback to REPL?
    debug::Bool                             # run nan_detection code?
    
    output::Bool                            # write a progress.txt?
    id::Union{String,Int}                   # identification of run, taken from ::OutputWriter
    run_path::String                        # path to run folder, taken from ::OutputWriter

    # PROGRESS
    progress_meter::ProgressMeter.Progress  # struct containing everything progress related
    progress_txt::Union{IOStream,Nothing}   # txt is a Nothing in case of no output

    # NaRS (Not-a-Real) AND OTHER MODEL STATE FEEDBACK
    nars_detected::Bool                     # did Infs/NaNs occur in the simulation?
end

function Feedback(  outputter::OutputWriter,
                    time_stepping::TimeStepper,
                    verbose::Bool=true,
                    debug::Bool=true)
    
    (;output) = outputter
    id = ""             # do this in initialize!(::OutputWriter, ...) to avoid folder-race conditions
    run_path = ""

    # PROGRESSMETER
    (;Δt_sec, n_timesteps) = time_stepping
    DT_IN_SEC[] = Δt_sec        # hack: redefine element in global constant dt_in_sec
                                # used to pass on the time step to ProgressMeter.speedstring
    desc = "Weather is speedy: "
    
    # show progress meter via `enabled` through verbose parameter
    # minus one to exclude first_timesteps! which contain compilation
    progress_meter = ProgressMeter.Progress(n_timesteps-1, enabled=verbose, showspeed=true; desc)
    progress_txt = nothing          # initialize with nothing, initialize in initialize!(::Feedback,...)

    nars_detected = false
    return Feedback(verbose,debug,
                    output,id,run_path,
                    progress_meter,progress_txt,
                    nars_detected)
end

"""Initialises the progress txt file."""
function initialize!(feedback::Feedback,model::ModelSetup)

    # reinitalize progress meter
    (;n, enabled, showspeed, desc) = feedback.progress_meter
    feedback.progress_meter = ProgressMeter.Progress(n;enabled,showspeed, desc)

    if feedback.output   # with netcdf output write progress.txt
        (; run_path, id) = feedback
        SG = model.spectral_grid
        L = model.time_stepping

        # create progress.txt file in run_????/
        progress_txt = open(joinpath(run_path,"progress.txt"),"w")
        s = "Starting SpeedyWeather.jl run $id on "*
                Dates.format(Dates.now(),Dates.RFC1123Format)
        write(progress_txt,s*"\n")
        write(progress_txt,"Integrating:\n")
        write(progress_txt,"$SG\n")
        write(progress_txt,"Time: $(L.n_days) days at Δt = $(L.Δt_sec)s\n")
        write(progress_txt,"\nAll data will be stored in $run_path\n")
        feedback.progress_txt = progress_txt
    end
end

function initialize!(pm::ProgressMeter.Progress)
    pm.counter = 0
    pm.tinit = time()
    pm.tlast = pm.tinit
end

"""Calls the progress meter and writes every 5% progress increase to txt."""
function progress!(feedback::Feedback)

    # update progress meter and unpack counter after update
    ProgressMeter.next!(feedback.progress_meter)    
    (; counter, n ) = feedback.progress_meter

    # write progress to txt file too
    if (counter/n*100 % 1) > ((counter+1)/n*100 % 1)  
        percent = round(Int,(counter+1)/n*100)      # % of time steps completed
        if feedback.output && (percent % 5 == 0)    # write every 5% step in txt 
            write(feedback.progress_txt,@sprintf("\n%3d%%",percent))
            r = remaining_time(feedback.progress_meter)
            write(feedback.progress_txt,", ETA: $r")

            time_elapsed = feedback.progress_meter.tlast - feedback.progress_meter.tinit
            s = speedstring(time_elapsed/counter,DT_IN_SEC[])
            write(feedback.progress_txt,", $s")

            feedback.nars_detected && write(feedback.progress_txt,", NaRs detected.")
            flush(feedback.progress_txt)
        end
    end
end

function progress!( feedback::Feedback,
                    progn::PrognosticVariables)
    progress!(feedback)
    feedback.debug && nar_detection!(feedback,progn)
end

"""Finalises the progress meter and the progress txt file."""
function progress_finish!(F::Feedback)
    ProgressMeter.finish!(F.progress_meter)
    
    if F.output     # write final progress to txt file
        time_elapsed = F.progress_meter.tlast - F.progress_meter.tinit
        s = "\nIntegration done in $(readable_secs(time_elapsed))."
        write(F.progress_txt,"\n$s\n")
        flush(F.progress_txt)
    end
end

"""Detect NaR (Not-a-Real) in the prognostic variables."""
function nar_detection!(feedback::Feedback,progn::PrognosticVariables)

    feedback.nars_detected && return nothing    # escape immediately if nans already detected
    i = feedback.progress_meter.counter         # time step
    nars_detected_here = false
    (; vor ) = progn.layers[end].timesteps[2] # only check for surface vorticity

    if ~nars_detected_here
        nars_vor = ~isfinite(vor[1])    # just check first mode
                                        # spectral transform propagates NaRs globally anyway
        nars_vor && @warn "NaR detected at time step $i"
        nars_detected_here |= nars_vor
    end

    feedback.nars_detected |= nars_detected_here
end

# adapted from ProgressMeter.jl
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

# adapted from ProgressMeter.jl
function speedstring(sec_per_iter,dt_in_sec)
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
const DT_IN_SEC = Ref(1800)

function ProgressMeter.speedstring(sec_per_iter,dt_in_sec=SpeedyWeather.DT_IN_SEC)
    speedstring(sec_per_iter,dt_in_sec[])
end

"""
    readable_secs(secs::Real) -> Dates.CompoundPeriod

Returns `Dates.CompoundPeriod` rounding to either (days, hours), (hours, minutes), (minutes,
seconds), or seconds with 1 decimal place accuracy for >10s and two for less.
E.g.
```julia
julia> readable_secs(12345)
3 hours, 26 minutes
```
"""
function readable_secs(secs::Real)
    millisecs = Dates.Millisecond(round(secs * 10 ^ 3))
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