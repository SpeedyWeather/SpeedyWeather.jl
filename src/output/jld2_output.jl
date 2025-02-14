export JLD2Output

"""Output writer for a JLD2 file that saves the PrognosticVariables
and DiagnosticVariables structs directly to a JLD2 file. All internal
scalings and units are still applied to these outputs. Fields are 
$(TYPEDFIELDS)"""
@kwdef mutable struct JLD2Output <: AbstractOutput

    # FILE OPTIONS
    active::Bool = false
    
    "[OPTION] path to output folder, run_???? will be created within"
    path::String = pwd()
    
    "[OPTION] run identification number/string"
    id::String = "0001"
    run_path::String = ""                   # will be determined in initalize!
    
    "[OPTION] name of the output jld2 file"
    filename::String = "output.jld2"
    
    "[OPTION] also write restart file if output==true?"
    write_restart::Bool = true
    pkg_version::VersionNumber = isnothing(pkgversion(SpeedyWeather)) ? v"0.0.0" : pkgversion(SpeedyWeather)

    "[OPTION] output frequency, time step"
    output_dt::Second = Second(DEFAULT_OUTPUT_DT)

    "[OPTION] will reopen and resave the file to merge everything in one big vector. Turn off if the file is too large for memory."
    merge_output::Bool = true

    "[OPTION] output the DiagnosticVariables as well" 
    output_diagnostic::Bool = true 

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
    println(io, "└ frequency: $(output.output_dt)")
end

"""$(TYPEDSIGNATURES)
Initialize JLD2 `output` by creating a JLD2 file. 
To be called just before the first timesteps."""
function initialize!(   
    output::JLD2Output,
    feedback::AbstractFeedback,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    
    output.active || return nothing     # exit immediately for no output
    
    # GET RUN ID, CREATE FOLDER
    # get new id only if not already specified
    output.id = get_run_id(output.path, output.id)
    output.run_path = create_output_folder(output.path, output.id) 
    
    feedback.id = output.id             # synchronize with feedback struct
    feedback.run_path = output.run_path
    feedback.progress_meter.desc = "Weather is speedy: run $(output.id) "
    feedback.output = true              # if output=true set feedback.output=true too!

    # OUTPUT FREQUENCY
    output.output_every_n_steps = max(1, round(Int,
            Millisecond(output.output_dt).value/model.time_stepping.Δt_millisec.value))
    output.output_dt = Second(round(Int, output.output_every_n_steps*model.time_stepping.Δt_sec))

    # RESET COUNTERS
    output.timestep_counter = 0         # time step counter
    output.output_counter = 0           # output step counter

    # CREATE JLD2 FILE
    (; run_path, filename) = output

    jld2_file = jldopen(joinpath(run_path, filename), "w") 
    output.jld2_file = jld2_file

    # write iniital condition 
    output!(output, Simulation(progn, diagn, model))

    # also export parameters into run????/parameters.txt
    parameters_txt = open(joinpath(output.run_path, "parameters.txt"), "w")
    for property in propertynames(model)
        println(parameters_txt, "model.$property")
        println(parameters_txt, getfield(model, property,), "\n")
    end
    close(parameters_txt)
end

Base.close(output::JLD2Output) = close(output.jld2_file)

function output!(output::JLD2Output, simulation::AbstractSimulation)
    output.output_counter += 1      # output counter increases when writing time
    i = output.output_counter

    (; active, jld2_file, output_diagnostic) = output 
    active || return nothing    # escape immediately for no jld2 output

    if output_diagnostic
        jld2_file["$i"] = (simulation.prognostic_variables, simulation.diagnostic_variables)
    else 
        jld2_file["$i"] = simulation.prognostic_variables
    end 
end 

function finalize!(
    output::JLD2Output,
    simulation::AbstractSimulation,
)   
    if output.merge_output
        merge_output(output)
    else  
        close(output)
    end 
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

    jldsave(joinpath(run_path, filename); output_vector)
end 
