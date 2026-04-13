export JLD2Output, ArrayOutput

"""Abstract supertype for output writers that save `Variables` (or subsets thereof)
at regular intervals. Subtypes must have fields: `active`, `path`, `run_prefix`, `id`,
`run_number`, `run_digits`, `overwrite`, `write_restart`, `write_parameters_txt`,
`write_progress_txt`, `output_dt`, `groups`, and `core::OutputWriterCore`."""
abstract type AbstractVariablesOutput <: AbstractOutput end

# GROUPS FILTERING
"""$(TYPEDSIGNATURES)
Filter `vars::Variables` to only include the groups specified in `output.groups`.
Returns the `Variables` struct directly when `groups` contains `:all`,
otherwise a NamedTuple with the selected groups."""
function filter_groups(vars::Variables, output::AbstractVariablesOutput)
    groups = output.groups
    :all in groups && return vars
    selected = Tuple(g for g in ALL_VARIABLE_GROUPS if g in groups)
    return NamedTuple{selected}(Tuple(getfield(vars, g) for g in selected))
end

"""Output writer for a JLD2 file that saves the Variables struct directly to a JLD2 file.
All internal scalings and units are still applied to these outputs. Fields are
$(TYPEDFIELDS)"""
@kwdef mutable struct JLD2Output <: AbstractVariablesOutput

    # FILE OPTIONS
    active::Bool = false

    "[OPTION] path to output folder, run_???? will be created within"
    path::String = pwd()

    "[OPTION] Prefix for run folder where data is stored, e.g. 'run_'"
    run_prefix::String = "run"

    "[OPTION] run identification, added between run_prefix and run_number"
    id::String = ""

    "[OPTION] run identification number, automatically determined if overwrite=false"
    run_number::Int = 1

    "[OPTION] run numbers digits"
    run_digits::Int = 4

    "[OPTION] Overwrite an existing run folder?"
    overwrite::Bool = false

    "[OPTION] name of the output jld2 file"
    filename::String = "output.jld2"

    "[OPTION] also write restart file if output=true?"
    write_restart::Bool = true

    "[OPTION] also write parameters txt file if output=true?"
    write_parameters_txt::Bool = true

    "[OPTION] also write progress txt file if output=true?"
    write_progress_txt::Bool = true

    "[OPTION] output frequency, time step"
    output_dt::Second = Second(DEFAULT_OUTPUT_DT)

    "[OPTION] will reopen and resave the file to merge everything in one big vector. Turn off if the file is too large for memory."
    merge_output::Bool = true

    "[OPTION] which variable groups to save, e.g. (:prognostic,) or (:prognostic, :grid). Use (:all,) for all groups."
    groups::Tuple{Vararg{Symbol}} = (:all,)

    "[DERIVED] shared output writer state (run folder, counters, output frequency)"
    core::OutputWriterCore = OutputWriterCore()

    jld2_file::Union{JLDFile, Nothing} = nothing
end

function Base.show(io::IO, output::JLD2Output)
    println(io, "JLD2Output")
    println(io, "├ status: $(output.active ? "active" : "inactive/uninitialized")")
    println(io, "├ write restart file: $(output.write_restart) (if active)")
    println(io, "├ path: $(joinpath(output.run_path, output.filename))")
    println(io, "├ groups: $(output.groups)")
    return println(io, "└ frequency: $(output.output_dt)")
end

"""$(TYPEDSIGNATURES)
Initialize JLD2 `output` by creating a JLD2 file.
To be called just before the first timesteps."""
function initialize!(
        output::JLD2Output,
        vars::Variables,
        model::AbstractModel,
    )
    initialize!(output.core, output, model) || return nothing

    # CREATE JLD2 FILE
    (; run_path, filename) = output
    jld2_file = jldopen(joinpath(run_path, filename), "w")
    output.jld2_file = jld2_file

    # write initial condition
    output_jld2!(output, Simulation(vars, model))
    return nothing
end

Base.close(output::JLD2Output) = close(output.jld2_file)

function output!(output::JLD2Output, simulation::AbstractSimulation)
    output!(output.core, output) || return nothing
    return output_jld2!(output, simulation)
end

function output_jld2!(output::JLD2Output, simulation::AbstractSimulation)
    output.output_counter += 1
    i = output.output_counter
    snapshot = filter_groups(simulation.variables, output)
    output.jld2_file["$i"] = on_architecture(CPU(), snapshot)
    return nothing
end

function finalize!(
        output::JLD2Output,
        simulation::AbstractSimulation,
    )
    if output.merge_output && output.output_counter > 0
        merge_output(output)
    else
        close(output)
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
We can't directly push to arrays in a JLD2 file or have extendable
dimensions. This routine rewrites the file to a single vector.
Might be turned off if the file doesn't fit into the memory or speed
is a concern.
"""
function merge_output(output::JLD2Output)
    (; output_counter, jld2_file, run_path, filename) = output

    output_vector = Vector{typeof(jld2_file["1"])}(undef, output_counter)

    for i in 1:output_counter
        output_vector[i] = jld2_file["$i"]
    end

    # close and overwrite old file
    close(jld2_file)

    return jldsave(joinpath(run_path, filename); output_vector)
end

"""Output writer that stores Variables snapshots in memory instead of writing
to a file. Access the stored snapshots via `output.output` after the simulation.
Otherwise follows the same logic as [`JLD2Output`](@ref). Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct ArrayOutput <: AbstractVariablesOutput

    active::Bool = false

    "[OPTION] also write restart file if output=true?"
    write_restart::Bool = false

    "[OPTION] also write parameters txt file if output=true?"
    write_parameters_txt::Bool = false

    "[OPTION] also write progress txt file if output=true?"
    write_progress_txt::Bool = false

    "[OPTION] output frequency, time step"
    output_dt::Second = Second(DEFAULT_OUTPUT_DT)

    "[OPTION] which variable groups to save, e.g. (:prognostic,) or (:prognostic, :grid). Use (:all,) for all groups."
    groups::Tuple{Vararg{Symbol}} = (:all,)

    "[DERIVED] shared output writer state (run folder, counters, output frequency)"
    core::OutputWriterCore = OutputWriterCore()

    "[DERIVED] vector of stored variable snapshots, populated during the simulation"
    output::Vector = []
end

# ArrayOutput has no run folder or file — stub out path-related option fields
# that AbstractOutput methods like set! may try to assign.
Base.getproperty(output::ArrayOutput, s::Symbol) = begin
    s in (:path, :run_prefix, :id, :run_number, :run_digits, :overwrite) && return _arrayoutput_default(s)
    s === :core && return getfield(output, :core)
    s in OUTPUT_WRITER_CORE_FIELDS && return getproperty(getfield(output, :core), s)
    return getfield(output, s)
end

_arrayoutput_default(s::Symbol) = s === :run_number ? 1 : s === :run_digits ? 4 : s === :overwrite ? false : ""

function Base.setproperty!(output::ArrayOutput, s::Symbol, v)
    # silently ignore path-related option fields that ArrayOutput doesn't use
    s in (:path, :run_prefix, :id, :run_number, :run_digits, :overwrite) && return v
    s === :core && return setfield!(output, :core, v)
    s in OUTPUT_WRITER_CORE_FIELDS && return setproperty!(getfield(output, :core), s, v)
    return setfield!(output, s, v)
end

function Base.show(io::IO, output::ArrayOutput)
    n = length(output.output)
    println(io, "ArrayOutput")
    println(io, "├ status: $(output.active ? "active" : "inactive/uninitialized")")
    println(io, "├ snapshots stored: $n")
    println(io, "├ groups: $(output.groups)")
    return println(io, "└ frequency: $(output.output_dt)")
end

"""$(TYPEDSIGNATURES)
Initialize `ArrayOutput` by pre-allocating the full output array and storing the initial condition."""
function initialize!(
        output::ArrayOutput,
        vars::Variables,
        model::AbstractModel,
    )
    initialize!(getfield(output, :core), output, model; create_folder = false) || return nothing

    # compute total number of output snapshots: IC + one per output_every_n_steps
    n_timesteps = vars.prognostic.clock.n_timesteps
    n_outputs = n_timesteps ÷ output.output_every_n_steps + 1   # +1 for the IC

    # pre-allocate the full vector with deep copies
    ic = deepcopy(filter_groups(vars, output))
    output_vector = Vector{typeof(ic)}(undef, n_outputs)
    output_vector[1] = ic
    setfield!(output, :output, output_vector)

    # counter already reset to 0 by initialize!, set to 1 for IC
    output.output_counter = 1
    return nothing
end

Base.close(::ArrayOutput) = nothing

function output!(output::ArrayOutput, simulation::AbstractSimulation)
    output!(output.core, output) || return nothing
    output.output_counter += 1
    i = output.output_counter
    getfield(output, :output)[i] = deepcopy(filter_groups(simulation.variables, output))
    return nothing
end

function finalize!(::ArrayOutput, ::AbstractSimulation)
    return nothing
end
