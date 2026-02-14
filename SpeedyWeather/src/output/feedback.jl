abstract type AbstractFeedback end

export Feedback

"""
Feedback struct that contains options and object for command-line feedback
like the progress meter.
$(TYPEDFIELDS)"""
@kwdef mutable struct Feedback <: AbstractFeedback
    "[OPTION] print feedback to REPL?, default is isinteractive(), true in interactive REPL mode"
    verbose::Bool = isinteractive()

    "[OPTION] check for NaNs in the prognostic variables"
    debug::Bool = true

    "[OPTION] Progress description"
    description::String = ""

    "[OPTION] show speed in progress meter?"
    showspeed::Bool = true

    "[DERIVED] struct containing everything progress related"
    progress_meter::ProgressMeter.Progress =
        ProgressMeter.Progress(1, enabled = verbose)

    "[DERIVED] did NaNs occur in the simulation?"
    nans_detected::Bool = false
end

function Base.show(io::IO, F::AbstractFeedback)
    println(io, "$(typeof(F)) <: AbstractFeedback")
    keys = propertynames(F)
    return print_fields(io, F, keys)
end

function Base.show(io::IO, P::ProgressMeter.Progress)
    println(io, "$(typeof(P)) <: ProgressMeter.AbstractProgress")
    keys = propertynames(P)
    return print_fields(io, P, keys)
end

"""
$(TYPEDSIGNATURES)
Initializes the a `Feedback` struct."""
function initialize!(feedback::Feedback, clock::Clock, model::AbstractModel)
    # set to false to recheck for NaNs
    feedback.nans_detected = false

    # hack: redefine element in global constant dt_in_sec
    # used to pass on the time step to ProgressMeter.speedstring
    DT_IN_SEC[] = model.time_stepping.Δt_sec

    # reinitalize progress meter, minus one to exclude first_timesteps! which contain compilation
    # only do now for benchmark accuracy
    (; showspeed, description, verbose) = feedback
    desc = description * (model.output.active ? " $(model.output.run_folder) " : " ")
    feedback.progress_meter = ProgressMeter.Progress(clock.n_timesteps - 1; 
        enabled = verbose,
        showspeed,
        desc,
        color = :blue,
        barlen = 10,
        barglyphs = ProgressMeter.BarGlyphs(" ━━  "),
        )
    return nothing
end

progress!(feedback::Feedback) = ProgressMeter.next!(feedback.progress_meter)

function progress!(feedback::Feedback, progn::PrognosticVariables, diagn::DiagnosticVariables)
    mod(feedback.progress_meter.core.counter, 100) == 0 && max_speed(diagn)
    progress!(feedback)
    feedback.debug && nan_detection!(feedback, progn)
    return nothing
end

"""
$(TYPEDSIGNATURES)
Finalises the progress meter and the progress txt file."""
finalize!(F::Feedback) = ProgressMeter.finish!(F.progress_meter)

"""$(TYPEDSIGNATURES)
Detect NaN (Not-a-Number, or Inf) in the prognostic variables."""
function nan_detection!(feedback::Feedback, progn::PrognosticVariables)
    feedback.nans_detected && return nothing            # escape immediately if nans already detected
    i = feedback.progress_meter.counter                 # time step
    GPUArrays.@allowscalar vor0 = progn.vor[2, end, 2]  # only check 1-0 mode of surface vorticity

    # just check first harmonic, spectral transform propagates NaNs globally anyway
    nans_detected_here = ~isfinite(vor0)
    nans_detected_here && @warn "NaN or Inf detected at time step $i"
    return feedback.nans_detected = nans_detected_here
end

"""
$(TYPEDSIGNATURES)
Define a ProgressMeter.speedstring method that also takes a time step
`dt_in_sec` to translate sec/iteration to days/days-like speeds."""
function speedstring(sec_per_iter, dt_in_sec)
    if sec_per_iter == Inf
        return "  N/A  days/day"
    end

    sim_time_per_time = dt_in_sec / sec_per_iter

    for (divideby, unit) in (
            (365 * 1_000, "millenia"),
            (365, "years"),
            (1, "days"),
            (1 / 24, "hours"),
        )
        if (sim_time_per_time / divideby) > 2
            return @sprintf "%5.2f %2s/day" (sim_time_per_time / divideby) unit
        end
    end
    return " <2 hours/days"
end

# hack: define global constant whose element will be changed in initialize_feedback
# used to pass on the time step to ProgressMeter.speedstring via calling this
# constant from the ProgressMeter module
const DT_IN_SEC = Ref(1.0)
const FEEDBACK_UMAX = Ref(0f0)

# "extend" the speedstring function from ProgressMeter by defining it for ::AbstractFloat
# not just ::Any to effectively overwrite it
function ProgressMeter.speedstring(sec_per_iter::AbstractFloat)
    dt_in_sec = SpeedyWeather.DT_IN_SEC[]   # pull global "constant"
    U = SpeedyWeather.FEEDBACK_UMAX[]
    return progress_string(sec_per_iter, dt_in_sec, U)
end

function progress_string(sec_per_iter, dt_in_sec, U)
    speed = speedstring(sec_per_iter, dt_in_sec)
    umax = @sprintf "U = %2i m/s" U
    return speed * ", " * umax
end

function max_speed(diagn::DiagnosticVariables)
    umin, umax = extrema(diagn.grid.u_grid)
    FEEDBACK_UMAX[] = max(abs(umin), abs(umax))
end

export ParametersTxt

"""ParametersTxt callback. Writes a parameters.txt file with all model parameters.
Options are $(TYPEDFIELDS)"""
@kwdef mutable struct ParametersTxt <: AbstractCallback
    "[OPTION] Path for parameters.txt file, uses model.output.run_path if not specified"
    path::String = ""

    "[OPTION] File name for parameters.txt file"
    filename::String = "parameters.txt"

    "[OPTION] Only write with model.output.active = true?"
    write_only_with_output::Bool = true
end

"""$(TYPEDSIGNATURES)
Initialize ParametersTxt by writing the model parameters (via defined show of model components) into a text file."""
function initialize!(parameters_txt::ParametersTxt, progn, diagn, model)

    # escape in case of no output
    parameters_txt.write_only_with_output && (model.output.active || return nothing)

    (; filename) = parameters_txt
    path = parameters_txt.path == "" ? model.output.run_path : parameters_txt.path
    mkpath(path)

    # also export parameters into run????/parameters.txt
    file = open(joinpath(path, filename), "w")
    for property in propertynames(model)
        println(file, "model.$property")
        println(file, getfield(model, property), "\n")
    end
    close(file)

    model.output.active || @info "Parameter summary written to $(joinpath(path, filename)) although output=false"
    return nothing
end

callback!(::ParametersTxt, args...) = nothing
finalize!(::ParametersTxt, arg...) = nothing


export ProgressTxt

"""ProgressTxt callback. Writes a progress.txt file with time stepping progress.
Options are $(TYPEDFIELDS)"""
@kwdef mutable struct ProgressTxt <: AbstractCallback
    "[OPTION] Path for progress.txt file, uses model.output.run_path if not specified"
    path::String = ""

    "[OPTION] File name for progress.txt file"
    filename::String = "progress.txt"

    "[OPTION] Only write with model.output.active = true?"
    write_only_with_output::Bool = true

    "[OPTION] Every n% of time steps write to progress.txt, default is 5%"
    every_n_percent::Int = 5

    "[DERIVED] IOStream for progress.txt file"
    file::IOStream = IOStream("")
end

"""$(TYPEDSIGNATURES)
Initializes the ProgressTxt callback by creating a progress.txt file and writing some initial information to it."""
function initialize!(progress_txt::ProgressTxt, progn, diagn, model)
    # escape in case of no output
    progress_txt.write_only_with_output && (model.output.active || return nothing)

    (; filename) = progress_txt
    path = progress_txt.path == "" ? model.output.run_path : progress_txt.path
    mkpath(path)

    (; run_folder, run_path) = model.output
    SG = model.spectral_grid
    L = model.time_stepping
    days = Second(progn.clock.period).value / (3600 * 24)

    # create progress.txt file in run_????/
    file = open(joinpath(path, filename), "w")
    s = "Starting SpeedyWeather.jl $run_folder on " *
        Dates.format(Dates.now(), Dates.RFC1123Format)
    write(file, s * "\n")
    write(file, "Integrating:\n")
    write(file, "$SG\n")
    write(file, "Time: $days days at Δt = $(L.Δt_sec)s\n")
    model.output.active && write(file, "\nAll data will be stored in $run_path\n")
    model.output.active || write(file, "\nNo output will be written (output=false)\n")
    progress_txt.file = file

    model.output.active || @info "Progress is being written to $(joinpath(path, filename)) although output=false"
    return nothing
end

"""$(TYPEDSIGNATURES)
Writes the time stepping progress to the progress.txt file every `every_n_percent` % of time steps."""
function callback!(progress_txt::ProgressTxt, progn, diagn, model)
    # escape in case of no output
    progress_txt.write_only_with_output && (model.output.active || return nothing)

    (; progress_meter, nans_detected) = model.feedback
    (; counter, n) = progress_meter
    (; file, every_n_percent) = progress_txt

    # occasionally write progress to txt file
    return if (counter / n * 100 % 1) > ((counter + 1) / n * 100 % 1)
        percent = round(Int, (counter + 1) / n * 100)             # % of time steps completed
        if (percent % every_n_percent == 0)                 # write every p% step in txt
            write(file, @sprintf("\n%3d%%", percent))
            r = remaining_time(progress_meter)
            write(file, ", ETA: $r")

            time_elapsed = progress_meter.tlast - progress_meter.tinit
            s = speedstring(time_elapsed / counter, DT_IN_SEC[])
            write(file, ", $s")

            nans_detected && write(file, ", NaN/Inf detected.")
            flush(file)
        end
    end
end

"""$(TYPEDSIGNATURES)
Finalizes the ProgressTxt callback by writing the total time taken to the progress.txt file and closing it."""
function finalize!(progress_txt::ProgressTxt, progn, diagn, model)
    # escape in case of no output
    progress_txt.write_only_with_output && (model.output.active || return nothing)

    (; file) = progress_txt
    (; progress_meter) = model.feedback

    time_elapsed = progress_meter.tlast - progress_meter.tinit
    s = "\nIntegration done in $(readable_secs(time_elapsed))."
    write(file, "\n$s\n")
    flush(file)
    close(file)
    return nothing
end

"""$(TYPEDSIGNATURES)
Estimates the remaining time from a `ProgresssMeter.Progress`. Adapted from ProgressMeter.jl"""
function remaining_time(p::ProgressMeter.Progress)
    elapsed_time = time() - p.tinit
    est_total_time = elapsed_time * (p.n - p.start) / (p.counter - p.start)
    if 0 <= est_total_time <= typemax(Int)
        eta_sec = round(Int, est_total_time - elapsed_time)
        eta = ProgressMeter.durationstring(eta_sec)
    else
        eta = "N/A"
    end
    return eta
end
