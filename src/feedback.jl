@with_kw mutable struct Feedback
    
    # DURATION ESTIMATE
    i::Int=0                                # time step
    i_start_durest::Int=2                   # first time steps to ignore for duration estimate
    n_steps_durest::Int=5                   # use that many time steps to extrapolate total duration
    t_start::Float64=time()                 # start time (sec since UNIX epoch 1/1/1970)
    t_start_durest::Float64=t_start         # start time for duration estimation (to be updated)
    t_end::Float64=t_start                  # end time (to be updated)
    progress_txt::Union{IOStream,Nothing}   # txt is a Nothing in case of no output
    
    # OUTPUT
    output::Bool                            # output to netCDF?
    n_timesteps::Int                        # number of time steps
    n_outputsteps::Int                      # number of time steps with output
    run_id::Int                             # run identification number
    run_path::String                        # output path plus run????/

    # NANS AND OTHER MODEL STATE FEEDBACK
    nans_detected::Bool=false               # did NaNs occur in the simulation?
end

"""Estimates the total time the model integration will take."""
function duration_estimate(feedback::Feedback)

    @unpack t_start_durest,i,n_timesteps,output,progress_txt = feedback

    # time estimates
    time_per_step = (time()-t_start_durest) / (i-1)
    time_total = round(Int,time_per_step*n_timesteps)
    time_to_go = round(Int,time_per_step*(n_timesteps-i))

    s1 = "Model integration will take approximately "*readable_secs(time_total)*","
    s2 = "and is hopefully done on "*
            Dates.format(Dates.now() + Dates.Second(time_to_go),Dates.RFC1123Format)

    # print time information inline
    print("\r\u1b[K")   
    println(s1)
    println(s2)

    if output           # print in txt file
        write(progress_txt,"\n"*s1*"\n")
        write(progress_txt,s2*"\n")
        flush(progress_txt)
    end
end

"""Initialises the progress txt file."""
function initialize_feedback(M::ModelSetup)
    @unpack output = M.parameters
    @unpack n_timesteps, n_outputsteps = M.constants

    if output   # with netcdf output
        @unpack NF,n_days,trunc,nlon,nlat = M.parameters

        run_id,run_path = get_run_id_path(M.parameters)     # create output folder and get its id and path
        
        # create progress.txt file in run????/
        progress_txt = open(joinpath(run_path,"progress.txt"),"w")
        s = "Starting SpeedyWeather.jl run $run_id on "*
                Dates.format(Dates.now(),Dates.RFC1123Format)
        println(s)                      # also write to REPL
        write(progress_txt,s*"\n")      # and in file

        # add some information on resolution and number format
        write(progress_txt,"Integrating $(n_days) days at a spectral resolution of"*
                            "T$trunc with $(nlon)x$(nlat) grid points.\n")
        write(progress_txt,"Number format is "*string(NF)*".\n")
        write(progress_txt,"All data will be stored in $run_path.\n")

        # also export parameters into run????/parameters.txt
        parameters_txt = open(joinpath(run_path,"parameters.txt"),"w")
        print(parameters_txt,M.parameters)
        close(parameters_txt)

    else        # no netcdf output
        println("Starting SpeedyWeather.jl on "*
                    Dates.format(Dates.now(),Dates.RFC1123Format)*" without output.")
        progress_txt = nothing      # for no ouput, allocate dummies for Feedback struct
        run_id = -1                 # dummy
        run_path = ""               # dummy
    end

    return Feedback(;progress_txt,output,n_timesteps,n_outputsteps,run_id,run_path)
end

"""Feedback function that calls duration estimate, nan_detection and progress."""
function feedback!( feedback::Feedback,     # Feedback struct
                    i::Int)                 # time integration time step index

    @unpack output, i_start_durest, n_steps_durest = feedback
    feedback.i = i      # update time step counter

    if output
        if i == i_start_durest                      # start the clock after ith time step
            feedback.t_start_durest = time()        # start time counter
        elseif i == i_start_durest+n_steps_durest   # after some time steps
            duration_estimate(feedback)             # estimate total model duration
        end
    end

    progress!(feedback)                             # feedback on progress in %
end

"""Finalises the progress txt file."""
function feedback_end!(feedback::Feedback)
    @unpack output,t_start,progress_txt = feedback
    
    t_end = time()          # time when the simulation ends
    feedback.t_end = t_end  # store in struct

    s = " Integration done in "*readable_secs(t_end-t_start)*"."
    println(s)
    if output
        write(progress_txt,"\n"*s[2:end]*"\n")  # close txt file with last output
        flush(progress_txt)
    end
end

"""Converts time step into percent for feedback."""
function progress!(feedback::Feedback)

    @unpack i,n_timesteps,progress_txt,output = feedback

    # update every 1% steps in REPL
    if (i/n_timesteps*100 % 1) > ((i+1)/n_timesteps*100 % 1)  
        percent = round(Int,(i+1)/n_timesteps*100)  # % of time steps completed
        print("\r\u1b[K")                           # remove previous p% in the REPL
        print("$percent%")                          # print new
        if output && (percent % 5 == 0)             # write every 5% step in txt 
            write(progress_txt,"\n$percent%")
            flush(progress_txt)
        end
    end
end