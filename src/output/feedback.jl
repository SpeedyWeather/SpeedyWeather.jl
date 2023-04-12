mutable struct Feedback <: AbstractFeedback
    verbose::Bool                           # print feedback to REPL?
    debug::Bool                             # run nan_detection code?
    output::Bool                            # store output (here: write to parameters and progress.txt?)

    # PROGRESS
    progress_meter::ProgressMeter.Progress  # struct containing everything progress related
    progress_txt::Union{IOStream,Nothing}   # txt is a Nothing in case of no output

    # NaRS (Not-a-Real) AND OTHER MODEL STATE FEEDBACK
    nars_detected::Bool                     # did Infs/NaNs occur in the simulation?
end

"""Initialises the progress txt file."""
function initialize_feedback(outputter::Output,M::ModelSetup)
    @unpack output, write_restart = outputter
    @unpack run_id, run_path = outputter

    if output   # with netcdf output write parameters.txt and progress.txt
        @unpack NF, n_days, trunc = M.parameters
        @unpack Grid, npoints = M.geometry

        # create progress.txt file in run????/
        progress_txt = open(joinpath(run_path,"progress.txt"),"w")
        s = "Starting SpeedyWeather.jl run $run_id on "*
                Dates.format(Dates.now(),Dates.RFC1123Format)
        write(progress_txt,s*"\n")      # and in file

        # add some information on resolution and number format
        write(progress_txt,"Integrating $(n_days) days at a spectral resolution of "*
                            "T$trunc on a $Grid with $npoints grid points.\n")
        write(progress_txt,"Number format is "*string(NF)*".\n")
        write(progress_txt,"All data will be stored in $run_path.\n")

        # also export parameters into run????/parameters.txt
        parameters_txt = open(joinpath(run_path,"parameters.txt"),"w")
        print(parameters_txt,M.parameters)
        close(parameters_txt)

    else        # no netcdf output
        progress_txt = nothing      # for no ouput, allocate dummies for Feedback struct
    end

    nans_detected = false           # don't check again if true

    # PROGRESSMETER
    @unpack verbose, debug = M.parameters
    @unpack n_timesteps = M.constants
    DT_IN_SEC[] = M.constants.Δt_sec        # hack: redefine element in global constant dt_in_sec
                                            # used to pass on the time step to ProgressMeter.speedstring
    desc = "Weather is speedy$(output ? " run $run_id: " : ": ")"
    
    # show progress meter via `enabled` through verbose parameter
    # one more time steps for the first Euler time step in first_timesteps!
    progress_meter = ProgressMeter.Progress(n_timesteps+1, enabled=verbose, showspeed=true; desc)

    return Feedback(verbose,debug,output,progress_meter,progress_txt,nans_detected)
end

"""Calls the progress meter and writes every 5% progress increase to txt."""
function progress!(feedback::Feedback)
    ProgressMeter.next!(feedback.progress_meter)    # update progress meter
    (;counter, n) = feedback.progress_meter         # unpack counter after update

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
    (;vor) = progn.layers[end].timesteps[2]     # only check for surface vorticity

    if ~nars_detected_here
        nars_vor = ~isfinite(vor[1])    # just check first mode
                                        # spectral transform propagates NaRs globally anyway
        nars_vor && @warn "NaR detected at time step $i"
        nars_detected_here |= nars_vor
    end

    feedback.nars_detected |= nars_detected_here
end

