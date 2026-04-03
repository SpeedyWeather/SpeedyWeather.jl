export JLD2Output, ArrayOutput

"""Abstract supertype for output writers that save `Variables` (or subsets thereof)
at regular intervals. Subtypes must have fields: `active`, `path`, `run_prefix`, `id`,
`run_number`, `run_digits`, `run_folder`, `run_path`, `overwrite`, `write_restart`,
`write_parameters_txt`, `write_progress_txt`, `output_dt`, `output_every_n_steps`,
`timestep_counter`, `output_counter`, `groups`."""
abstract type AbstractVariablesOutput <: AbstractOutput end

# GROUPS FILTERING

const ALL_VARIABLE_GROUPS = (:prognostic, :grid, :tendencies, :dynamics, :parameterizations, :particles, :scratch)

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

# SHARED INITIALIZE

"""$(TYPEDSIGNATURES)
Shared initialization for `AbstractVariablesOutput` subtypes: determine run folder,
compute output frequency, reset counters, add standard callbacks."""
function initialize_variables_output!(
        output::AbstractVariablesOutput,
        model::AbstractModel,
    )
    output.active || return false

    # GET RUN ID, CREATE FOLDER (only for outputs that write to disk)
    if hasproperty(output, :filename)
        determine_run_folder!(output)
        create_run_folder!(output)
    end

    # OUTPUT FREQUENCY
    output.output_every_n_steps = max(
        1, round(
            Int,
            Millisecond(output.output_dt).value / model.time_stepping.Δt_millisec.value
        )
    )
    output.output_dt = Second(round(Int, output.output_every_n_steps * model.time_stepping.Δt_sec))

    # RESET COUNTERS
    output.timestep_counter = 0
    output.output_counter = 0

    # CALLBACKS
    output.write_parameters_txt && add!(model.callbacks, :parameters_txt => ParametersTxt())
    output.write_progress_txt && add!(model.callbacks, :progress_txt => ProgressTxt())
    output.write_restart && add!(model.callbacks, :restart_file => RestartFile())

    return true
end

# SHARED OUTPUT STEPPING

"""$(TYPEDSIGNATURES)
Shared output stepping logic: increment timestep counter and check whether output
should be written on this step. Returns `true` if output should proceed."""
function should_output!(output::AbstractVariablesOutput)
    output.timestep_counter += 1
    (; active, output_every_n_steps, timestep_counter) = output
    active || return false
    timestep_counter % output_every_n_steps == 0 || return false
    return true
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

    "[DERIVED] folder name where data is stored, determined at initialize!"
    run_folder::String = ""

    "[DERIVED] full path to folder where data is stored, determined at initialize!"
    run_path::String = ""

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

    # TIME STEPS AND COUNTERS (initialize later)
    output_every_n_steps::Int = 0           # output frequency
    timestep_counter::Int = 0               # time step counter
    output_counter::Int = 0                 # output step counter

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
    initialize_variables_output!(output, model) || return nothing

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
    should_output!(output) || return nothing
    return output_jld2!(output, simulation)
end

function output_jld2!(output::JLD2Output, simulation::AbstractSimulation)
    output.output_counter += 1
    i = output.output_counter
    output.jld2_file["$i"] = filter_groups(simulation.variables, output)
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

    # TIME STEPS AND COUNTERS (initialize later)
    output_every_n_steps::Int = 0           # output frequency
    timestep_counter::Int = 0               # time step counter
    output_counter::Int = 0                 # output step counter

    "[DERIVED] vector of stored variable snapshots, populated during the simulation"
    output::Vector = []
end

# ArrayOutput has no run folder or file, but AbstractOutput methods
# like set! need these fields. Provide stubs.
Base.getproperty(output::ArrayOutput, s::Symbol) = begin
    if s === :run_folder
        return ""
    elseif s === :run_path
        return ""
    elseif s === :path
        return ""
    else
        return getfield(output, s)
    end
end

function Base.setproperty!(output::ArrayOutput, s::Symbol, v)
    # silently ignore path-related assignments from set!
    s in (:run_folder, :run_path, :path) && return v
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
Initialize `ArrayOutput` by clearing any previous output and setting up counters."""
function initialize!(
        output::ArrayOutput,
        vars::Variables,
        model::AbstractModel,
    )
    initialize_variables_output!(output, model) || return nothing

    # clear and store initial condition
    empty!(getfield(output, :output))
    output.output_counter += 1
    snapshot = deepcopy(filter_groups(vars, output))
    push!(getfield(output, :output), snapshot)
    return nothing
end

Base.close(::ArrayOutput) = nothing

function output!(output::ArrayOutput, simulation::AbstractSimulation)
    should_output!(output) || return nothing
    output.output_counter += 1
    snapshot = deepcopy(filter_groups(simulation.variables, output))
    push!(getfield(output, :output), snapshot)
    return nothing
end

function finalize!(::ArrayOutput, ::AbstractSimulation)
    return nothing
end
