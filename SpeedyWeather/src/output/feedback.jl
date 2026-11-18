abstract type AbstractFeedback <: AbstractModelComponent end

export Feedback, ProgressElements

"""
Feedback of a running simulation: the progress lines printed to the REPL (built from `elements`, see
`ProgressElements`) and the detection of NaNs (`nan_detection`, see `NaNDetector`). Options are
$(TYPEDFIELDS)"""
@kwdef mutable struct Feedback{E <: Tuple, N <: Union{AbstractNaNDetector, Nothing}} <: AbstractFeedback
    "[OPTION] print the progress lines to the REPL?, default is isinteractive(), true in interactive REPL mode"
    verbose::Bool = isinteractive()

    "[OPTION] Description in front of the progress lines"
    description::String = ""

    "[OPTION] Minimum wallclock time [s] between redraws of the progress lines"
    feedback_dt::Float32 = 0.1

    """[OPTION] Elements of the progress lines, see `ProgressElements`, they are shown with
    `print_element(element, progress, vars, model)` on every redraw. E.g. a single line with a bar of
    fixed length `Feedback(elements = (ProgressElements.Percentage(), ProgressElements.Bar(length = 20), ProgressElements.ETA()))`"""
    elements::E = ProgressElements.default_elements()

    "[OPTION] Separator between the elements of the progress lines"
    separator::String = ", "

    "[OPTION] Detection of NaNs in the simulation, which stops it, `nothing` to disable"
    nan_detection::N = NaNDetector()

    "[DERIVED] progress meter, created in `initialize!`"
    progress_meter::ProgressMeters.Progress = ProgressMeters.Progress(1, enabled = false)
end

# the default outer constructor of a parametric struct doesn't convert its arguments (e.g. feedback_dt = 0.1)
Feedback(verbose, description, feedback_dt, elements::E, separator, nan_detection::N, progress_meter) where {E, N} =
    Feedback{E, N}(verbose, description, feedback_dt, elements, separator, nan_detection, progress_meter)

# variables of the NaN detection and the elements, e.g. the scratch vector of `VerticalCourantNumber`
variables(feedback::Feedback) = (
    variables(feedback.nan_detection)...,
    Iterators.flatten(map(variables, feedback.elements))...,
)

"""
$(TYPEDSIGNATURES)
Initializes the `Feedback`: a new progress meter for the time steps of this run and the NaN detection."""
function initialize!(feedback::Feedback, vars::Variables, model::AbstractModel)
    (; clock) = vars.prognostic
    (; description, verbose, feedback_dt, elements, separator) = feedback
    model.output.active && (description = strip(description * " " * model.output.run_folder))
    feedback.progress_meter = ProgressMeters.Progress(
        Int(clock.n_steps);     # use time stepper steps (regardless Δt) not time steps of size Δt, Int for Reactant's numbers
        enabled = verbose,
        description,
        elements,
        separator,
        color = :blue,
        dt = feedback_dt,
    )
    isnothing(feedback.nan_detection) || initialize!(feedback.nan_detection, vars, model)
    return nothing
end

"""$(TYPEDSIGNATURES)
Feedback after `nsteps` time steps (1 normally, several after a compiled program of several time steps
with Reactant): advance the progress meter, which redraws the progress lines every `feedback_dt` seconds
with `print_element(element, progress, vars, model)`, and detect NaNs. After several time steps only the
last one can be checked for NaNs, so that's done immediately, as on the last time step."""
function progress!(feedback::Feedback, vars::Variables, model::AbstractModel, nsteps::Integer = 1)
    progress_meter = feedback.progress_meter
    ProgressMeters.next!(progress_meter, vars, model; step = nsteps)
    check = progress_meter.counter >= progress_meter.n || nsteps > 1
    detect_nans!(feedback.nan_detection, vars, model, progress_meter.counter; check)
    return nothing
end

# fallback for feedback = nothing
progress!(::Nothing, vars::Variables, model::AbstractModel, nsteps::Integer = 1) = nothing

"""$(TYPEDSIGNATURES)
Were NaNs detected in the simulation? `false` without feedback or NaN detection."""
nans_detected(feedback::Feedback) = nans_detected(feedback.nan_detection)

"""
$(TYPEDSIGNATURES)
Finishes the progress meter, drawing the final progress lines."""
finalize!(feedback::Feedback, vars::Variables, model::AbstractModel) =
    ProgressMeters.finish!(feedback.progress_meter, vars, model)

# fallback if feedback is set to nothing
finalize!(::Nothing, vars::Variables, model::AbstractModel) = nothing

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
function initialize!(parameters_txt::ParametersTxt, vars, model)

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
function initialize!(progress_txt::ProgressTxt, vars, model)
    # escape in case of no output
    progress_txt.write_only_with_output && (model.output.active || return nothing)

    (; filename) = progress_txt
    path = progress_txt.path == "" ? model.output.run_path : progress_txt.path
    mkpath(path)

    (; run_folder, run_path) = model.output
    SG = model.spectral_grid
    L = model.time_stepping
    days = Second(vars.prognostic.clock.period).value / (3600 * 24)

    # create progress.txt file in run_????/
    file = open(joinpath(path, filename), "w")
    s = "Starting SpeedyWeather.jl $run_folder on " *
        Dates.format(Dates.now(), Dates.RFC1123Format)
    write(file, s * "\n")
    write(file, "Integrating:\n")
    write(file, "$SG\n")
    write(file, "Time: $days days at Δt = $(L.Δt)s\n")
    model.output.active && write(file, "\nAll data will be stored in $run_path\n")
    model.output.active || write(file, "\nNo output will be written (output=false)\n")
    progress_txt.file = file

    model.output.active || @info "Progress is being written to $(joinpath(path, filename)) although output=false"
    return nothing
end

"""$(TYPEDSIGNATURES)
Writes the time stepping progress to the progress.txt file every `every_n_percent` % of time steps."""
function callback!(progress_txt::ProgressTxt, vars, model)
    # escape in case of no output
    progress_txt.write_only_with_output && (model.output.active || return nothing)
    isnothing(model.feedback) && return nothing

    (; progress_meter) = model.feedback
    (; counter, n) = progress_meter
    (; file, every_n_percent) = progress_txt

    # occasionally write progress to txt file
    if (counter / n * 100 % 1) > ((counter + 1) / n * 100 % 1)
        percent = round(Int, (counter + 1) / n * 100)       # % of time steps completed
        if (percent % every_n_percent == 0)                 # write every p% step in txt
            write(file, @sprintf("\n%3d%%", percent))
            t = time()
            r = ProgressMeters.durationstring(ProgressMeters.eta(progress_meter, t))
            write(file, ", ETA: $r")

            time_elapsed = t - progress_meter.tinit
            s = ProgressElements.speedstring(time_elapsed / max(1, counter), model.time_stepping.Δt)
            write(file, ", $s")

            nans_detected(model.feedback) && write(file, ", NaN/Inf detected.")
            flush(file)
        end
    end
    return nothing
end

"""$(TYPEDSIGNATURES)
Finalizes the ProgressTxt callback by writing the total time taken to the progress.txt file and closing it."""
function finalize!(progress_txt::ProgressTxt, vars, model)
    # escape in case of no output
    progress_txt.write_only_with_output && (model.output.active || return nothing)
    isnothing(model.feedback) && return nothing

    (; file) = progress_txt
    (; progress_meter) = model.feedback

    time_elapsed = time() - progress_meter.tinit
    s = "\nIntegration done in $(readable_secs(time_elapsed))."
    write(file, "\n$s\n")
    flush(file)
    close(file)
    return nothing
end
