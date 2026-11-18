abstract type AbstractFeedback <: AbstractModelComponent end

export Feedback, ProgressElements

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

    "[OPTION] Minimum wallclock time between feedback updates"
    feedback_dt::Float32 = 0.1

    """[OPTION] Interval in time steps between NaN checks. NaNs are recorded on every time step (without
    synchronizing a GPU) but only checked every `interval` steps, the warning then reports the time step
    at which they occurred"""
    interval::Int = 50

    """[OPTION] Elements of the progress line, see `ProgressElements`, `nothing` for
    `ProgressElements.default_elements(architecture)`. Elements specific to SpeedyWeather are bound to
    the simulation in `initialize!`, e.g. for a single line with a bar of fixed length
    `Feedback(elements = (ProgressElements.Percentage(), ProgressElements.Bar(length = 20), ProgressElements.ETA()))`"""
    elements::Union{Tuple, Nothing} = nothing

    "[OPTION] Separator between the elements of the progress line"
    separator::String = ", "

    "[DERIVED] progress meter, created in `initialize!`"
    progress_meter::ProgressMeters.Progress = ProgressMeters.Progress(1, enabled = false)

    "[DERIVED] did NaNs occur in the simulation?"
    nans_detected::Bool = false

    "[DERIVED] time step at which NaNs were first detected, 0 if none"
    nans_detected_step::Int = 0

    "[DERIVED] simulation time at which NaNs were first detected"
    nans_detected_time::DateTime = DateTime(2000)

    "[DERIVED] time steps recorded in `vars.dynamics.feedback.nan_record` (on the CPU)"
    nan_record_steps::Vector{Int} = zeros(Int, 1)

    "[DERIVED] simulation times of the time steps recorded in `nan_record`"
    nan_record_times::Vector{DateTime} = fill(DateTime(2000), 1)

    "[DERIVED] number of time steps recorded since the last check"
    nan_record_count::Int = 0
end

# variables needed by the elements of the progress line, e.g. the scratch vector of `VerticalCourantNumber`,
# and the record of the last `interval` time steps for the NaN checks
variables(feedback::Feedback) = (
    DynamicsVariable(
        :nan_record, VectorDim(max(1, feedback.interval)), namespace = :feedback,
        desc = "1 if the 1-0 mode of the surface vorticity was NaN/Inf in the recorded time step, 0 otherwise", units = "1"
    ),
    Iterators.flatten(map(variables, something(feedback.elements, ProgressElements.default_elements())))...,
)

"""
$(TYPEDSIGNATURES)
Initializes the a `Feedback` struct."""
function initialize!(feedback::Feedback, variables::Variables, model::AbstractModel)
    (; clock) = variables.prognostic

    # set to false to recheck for NaNs, records of the last `interval` time steps for the NaN checks
    feedback.nans_detected = false
    feedback.nans_detected_step = 0
    interval = length(variables.dynamics.feedback.nan_record)
    feedback.nan_record_steps = zeros(Int, interval)
    feedback.nan_record_times = fill(clock.time, interval)
    feedback.nan_record_count = 0

    # reinitalize progress meter, minus one to exclude first_timesteps! which contain compilation
    # only do now for benchmark accuracy
    (; description, verbose, feedback_dt, separator) = feedback
    elements = something(feedback.elements, ProgressElements.default_elements(model.architecture))
    elements = ProgressElements.bind_elements(elements, variables, model)
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

    return nothing
end

progress!(feedback::Feedback) = ProgressMeters.next!(feedback.progress_meter)

function progress!(feedback::Feedback, vars::Variables, model::AbstractModel)
    progress!(feedback)
    feedback.debug || return nothing

    # record the 1-0 mode on every time step, but check it for NaNs only every `interval` steps
    # (and on the last) as that requires copying it to the CPU, synchronizing a GPU
    record_nans!(feedback, vars, model)
    (; counter, n) = feedback.progress_meter
    (feedback.nan_record_count == length(feedback.nan_record_steps) || counter >= n) &&
        nan_detection!(feedback, vars)
    return nothing
end

# fallback for feedback = nothing
progress!(::Nothing, vars::Variables, model::AbstractModel) = nothing

"""$(TYPEDSIGNATURES)
Feedback after `nsteps` time steps at once, e.g. after a compiled program of several time steps
with Reactant: advance the progress meter by `nsteps` and check for NaNs (if `feedback.debug`)."""
function progress!(feedback::Feedback, vars::Variables, model::AbstractModel, nsteps::Integer)
    ProgressMeters.next!(feedback.progress_meter; step = nsteps)
    (feedback.debug && !feedback.nans_detected) || return nothing

    # only the last of the `nsteps` time steps can be checked, directly (no record, as it's not traced)
    vor = get_prognostic_step(vars.prognostic.vorticity, model.time_stepping, feedback)
    GPUArrays.@allowscalar vor_10 = vor[2, end]     # 1-0 mode of the surface vorticity
    isfinite(vor_10) || nans_detected!(feedback, feedback.progress_meter.counter, vars.prognostic.clock.time)
    return nothing
end

progress!(::Nothing, vars::Variables, model::AbstractModel, nsteps::Integer) = nothing

"""
$(TYPEDSIGNATURES)
Finalises the progress meter and the progress txt file."""
finalize!(F::Feedback) = ProgressMeters.finish!(F.progress_meter)

# fallback if feedback is set to nothing
finalize!(::Nothing) = nothing

"""$(TYPEDSIGNATURES)
Record the 1-0 mode of the surface vorticity (on the device, without synchronizing a GPU) together
with the time step and simulation time, for [`nan_detection!`](@ref)."""
function record_nans!(feedback::Feedback, vars::Variables, model::AbstractModel)
    feedback.nans_detected && return nothing        # stop recording once NaNs were detected
    k = feedback.nan_record_count += 1
    vor = get_prognostic_step(vars.prognostic.vorticity, model.time_stepping, feedback)
    # record whether the 1-0 mode of the surface (last layer) vorticity is finite, the spectral
    # transform propagates NaNs globally, so checking one coefficient is enough
    record_nans!(vars.dynamics.feedback.nan_record, k, vor.data, size(vor.data, 2))
    feedback.nan_record_steps[k] = feedback.progress_meter.counter
    feedback.nan_record_times[k] = vars.prognostic.clock.time
    return nothing
end

# nan_record[k] = !isfinite(vor[2, layer]) on the device, without synchronizing a GPU
function record_nans!(nan_record::AbstractVector, k::Integer, vor::AbstractArray, layer::Integer)
    view(nan_record, k:k) .= .!isfinite.(view(vor, 2, layer:layer))
    return nothing
end

"""$(TYPEDSIGNATURES)
Detect NaN (Not-a-Number, or Inf) in the time steps recorded since the last check (see
[`record_nans!`](@ref)) and warn with the time step and simulation time at which they first occurred."""
function nan_detection!(feedback::Feedback, vars::Variables)
    k = feedback.nan_record_count
    feedback.nan_record_count = 0
    (feedback.nans_detected || k == 0) && return nothing    # escape immediately if NaNs already detected
    record = Array(view(vars.dynamics.feedback.nan_record, 1:k))    # copy to the CPU (synchronizes a GPU)
    first_nan = findfirst(!iszero, record)
    isnothing(first_nan) && return nothing

    nans_detected!(feedback, feedback.nan_record_steps[first_nan], feedback.nan_record_times[first_nan])
    return nothing
end

# store and report time `step` and simulation `time` at which NaNs were detected
function nans_detected!(feedback::Feedback, step::Integer, time::DateTime)
    feedback.nans_detected = true
    feedback.nans_detected_step = step
    feedback.nans_detected_time = time
    @warn "NaN or Inf detected at time step $step ($time)"
    return nothing
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

    (; progress_meter, nans_detected) = model.feedback
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

            nans_detected && write(file, ", NaN/Inf detected.")
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
