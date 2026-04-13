export OutputWriterCore

abstract type AbstractOutput <: AbstractModelComponent end
abstract type AbstractOutputVariable <: AbstractModelComponent end

"""Shared runtime state for output writers: derived file paths, output frequency,
and step/output counters. Embed as `core::OutputWriterCore` in every output writer.
Option fields (path, id, overwrite, output_dt, …) stay as flat fields on the writer itself.

Call `initialize!(core, output, model)` inside the writer's `initialize!` method to
populate all derived fields from the writer's option fields and the model time step.
Call `output!(core, output)` each time step to check if output should be written.
Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct OutputWriterCore
    "[DERIVED] folder name where data is stored, determined at initialize!"
    run_folder::String = ""

    "[DERIVED] full path to folder where data is stored, determined at initialize!"
    run_path::String = ""

    "[DERIVED] output frequency in time steps, computed at initialize!"
    output_every_n_steps::Int = 0

    "[DERIVED] time step counter, incremented every time step"
    timestep_counter::Int = 0

    "[DERIVED] output step counter, incremented every time output is written"
    output_counter::Int = 0
end

const OUTPUT_WRITER_CORE_FIELDS = fieldnames(OutputWriterCore)

"""$(TYPEDSIGNATURES)
Forward `OutputWriterCore` field reads transparently from any `AbstractOutput` that
embeds a `core::OutputWriterCore` field, so callers can write `output.run_path` etc."""
function Base.getproperty(output::AbstractOutput, s::Symbol)
    s === :core && return getfield(output, :core)
    s in OUTPUT_WRITER_CORE_FIELDS && return getproperty(getfield(output, :core), s)
    return getfield(output, s)
end

"""$(TYPEDSIGNATURES)
Forward `OutputWriterCore` field writes transparently from any `AbstractOutput` that
embeds a `core::OutputWriterCore` field, so callers can write `output.run_path = ...` etc."""
function Base.setproperty!(output::AbstractOutput, s::Symbol, v)
    s === :core && return setfield!(output, :core, v)
    s in OUTPUT_WRITER_CORE_FIELDS && return setproperty!(getfield(output, :core), s, v)
    return setfield!(output, s, v)
end

"""$(TYPEDSIGNATURES)
Initialize the `OutputWriterCore` from the option fields of `output` and `model`:
determine the run folder (skipped for in-memory writers when `create_folder=false`),
compute the output frequency from `model.time_stepping`, reset counters, and add
standard callbacks (parameters txt, progress txt, restart file) if enabled.
Returns `true` if output is active, `false` otherwise."""
function initialize!(
        core::OutputWriterCore,
        output::AbstractOutput,
        model::AbstractModel;
        create_folder::Bool = true,
    )
    output.active || return false

    # GET RUN ID, CREATE FOLDER (only for outputs that write to disk)
    if create_folder
        determine_run_folder!(output)
        create_run_folder!(output)
    end

    # OUTPUT FREQUENCY
    core.output_every_n_steps = max(
        1, round(
            Int,
            Millisecond(output.output_dt).value / model.time_stepping.Δt_millisec.value
        )
    )
    output.output_dt = Second(round(Int, core.output_every_n_steps * model.time_stepping.Δt_sec))

    # RESET COUNTERS
    core.timestep_counter = 0
    core.output_counter = 0

    # CALLBACKS
    output.write_parameters_txt && add!(model.callbacks, :parameters_txt => ParametersTxt())
    output.write_progress_txt && add!(model.callbacks, :progress_txt => ProgressTxt())
    output.write_restart && add!(model.callbacks, :restart_file => RestartFile())

    return true
end

"""$(TYPEDSIGNATURES)
Increment the timestep counter and return `true` if output should be written on
this step (i.e. `active` is true and the counter is a multiple of `output_every_n_steps`)."""
function output!(core::OutputWriterCore, output::AbstractOutput)
    core.timestep_counter += 1
    output.active || return false
    core.timestep_counter % core.output_every_n_steps == 0 || return false
    return true
end
