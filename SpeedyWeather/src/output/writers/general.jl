abstract type AbstractOutput <: AbstractModelComponent end
abstract type AbstractOutputVariable <: AbstractModelComponent end

# default number format for output
const DEFAULT_OUTPUT_NF = Float32
const DEFAULT_OUTPUT_INTERVAL = Hour(6)
const OUTPUT_VARIABLES_DICT = Dict{Symbol, AbstractOutputVariable}
OutputVariablesDict() = OUTPUT_VARIABLES_DICT()

const DEFAULT_MISSING_VALUE = NaN
const DEFAULT_COMPRESSION_LEVEL = 1
const DEFAULT_SHUFFLE = false
const DEFAULT_KEEPBITS = 15

"""$(TYPEDSIGNATURES)
Add `outputvariables` to a dictionary defining the variables subject to NetCDF output."""
function add!(D::OUTPUT_VARIABLES_DICT, outputvariables::AbstractOutputVariable...)
    for outputvariable in outputvariables   # loop over all variables in arguments
        key = Symbol(outputvariable.name)   # use name as key::Symbol
        D[key] = outputvariable
    end
    return D
end

"""$(TYPEDSIGNATURES)
Add `outputvariables` to the dictionary in `output::NetCDFOutput`, i.e. at `output.variables`."""
function add!(output::AbstractOutput, outputvariables::AbstractOutputVariable...)
    add!(output.variables, outputvariables...)
    return output
end

"""$(TYPEDSIGNATURES)
Add `outputvariables` to the dictionary in `output::NetCDFOutput` of `model`, i.e. at `model.output.variables`."""
function add!(model::AbstractModel, outputvariables::AbstractOutputVariable...)
    add!(model.output, outputvariables...)
    return model.output
end

# also allow for tuples by splatting, similar to add!(::AbstractModel, ::Tuple)
add!(output::AbstractOutput, tuple::Tuple) = add!(output, tuple...)

"""$(TYPEDSIGNATURES)
Delete output variables from `output` by their (short name) (Symbol or String), corresponding
to the keys in the dictionary."""
function Base.delete!(output::AbstractOutput, keys::Union{String, Symbol}...)
    for key in keys
        delete!(output.variables, Symbol(key))
    end
    return output
end

"""$(TYPEDSIGNATURES)
Add default variables to output for a `Barotropic` model: Vorticity, zonal and meridional velocity."""
function add_default!(
        output_variables::OUTPUT_VARIABLES_DICT,
        Model::Type{<:Barotropic},
    )
    return add!(output_variables, VorticityOutput(), ZonalVelocityOutput(), MeridionalVelocityOutput())
end

"""$(TYPEDSIGNATURES)
Add default variables to output for a `ShallowWater` model, same as for a `Barotropic` model but also
the interface displacement."""
function add_default!(
        variables::Dict{Symbol, AbstractOutputVariable},
        Model::Type{<:ShallowWater},
    )
    add_default!(variables, Barotropic)
    return add!(variables, InterfaceDisplacementOutput())
end

"""$(TYPEDSIGNATURES)
Add default variables to output for a `PrimitiveDry` model, same as for a `Barotropic` model but also
the surface pressure and temperature."""
function add_default!(
        variables::Dict{Symbol, AbstractOutputVariable},
        Model::Type{<:PrimitiveDry},
    )
    add_default!(variables, Barotropic)
    return add!(variables, MeanSeaLevelPressureOutput(), TemperatureOutput())
end

"""$(TYPEDSIGNATURES)
Add default variables to output for a `PrimitiveWet` model, same as for a `PrimitiveDry` model but also
the specific humidity."""
function add_default!(
        variables::Dict{Symbol, AbstractOutputVariable},
        Model::Type{<:PrimitiveWet},
    )
    add_default!(variables, PrimitiveDry)
    return add!(variables, HumidityOutput())
end

function set!(output::AbstractOutput; active, reset_path = true)
    output.active = active
    if reset_path
        output.run_folder = ""
        output.run_path = ""
    end
    return nothing
end

# fallback for nothing output
set!(::Nothing; active, reset_path = true) = nothing

# fallback for nothing output
initialize!(::Nothing, ::Union{AbstractFeedback, Nothing}, ::Variables, ::AbstractModel) = nothing

# in case of no output nothing to close
Base.close(::Nothing) = nothing

"""$(TYPEDSIGNATURES)
Writes the variables from `vars` of time step `i` at time `time` into `output.netcdf_file`.
Simply escapes for no netcdf output or if output shouldn't be written on this time step.
Interpolates onto output grid and resolution as specified in `output`, converts to output
number format, truncates the mantissa for higher compression and applies lossless compression."""
function output!(output::AbstractOutput, simulation::AbstractSimulation)
    output!(output.core, output) || return nothing

    (; clock) = simulation.variables.prognostic
    output!(output, clock.time)                                         # increase counter, write time
    output!(output, output.variables, simulation)                       # write variables
    return output
end

# fallback for nothing output
function output!(::Nothing, ::AbstractSimulation)
    return nothing
end

"""$(TYPEDSIGNATURES)
Loop over every variable in `output.variables` to call the respective `output!` method
to write into the `output.netcdf_file`."""
function output!(
        output::AbstractOutput,
        output_variables::OUTPUT_VARIABLES_DICT,
        simulation::AbstractSimulation,
    )
    for var in values(output_variables)
        output!(output, var, simulation)
    end
    return
end

function finalize!(
        output::AbstractOutput,
        simulation::AbstractSimulation,
    )
    if output.active    # only finalize if active otherwise output.netcdf_file is nothing
        for var in values(output.variables)
            finalize!(output, var, simulation)
        end
    end
    return close(output)
end

finalize!(::Nothing, ::AbstractSimulation) = nothing

# default finalize method for output variables
function finalize!(
        output::AbstractOutput,
        var::AbstractOutputVariable,
        args...
    )
    # do nothing unless finalize! is defined for var <: AbstractOutputVariable
    return nothing
end

"""$(TYPEDSIGNATURES)
Checks existing folders in `path` and determine `run_number`by counting up.
E.g. if folder `run_0001` exists then `run_number` is 2.
Does not create a folder for the returned run id."""
function determine_run_folder!(output::AbstractOutput)
    (; run_prefix, id, run_digits) = output
    fmt = Printf.Format("%0$(run_digits)d")

    if !output.overwrite    # if not overwrite determine the next run number
        # try current run folder name, reset run number to 1
        output.run_number = 1
        run_folder = run_folder_name(run_prefix, id, output.run_number; fmt)
        while run_folder in readdir(output.path)    # if run folder already exists, increase run_number
            output.run_number += 1
            run_folder = run_folder_name(run_prefix, id, output.run_number; fmt)
        end
    end # else use the run_number exactly as specified in output.run_number

    return output.run_folder = run_folder_name(run_prefix, id, output.run_number; fmt)
end

"""$(TYPEDSIGNATURES)
Creates a new folder `prefix_id_number` with the identification `id`. Also returns the full path
`run_path` of that folder."""
function create_run_folder!(output::AbstractOutput)
    run_path = joinpath(output.path, output.run_folder)

    # actually create the folder
    # unless overwrite is true and the folder exists
    # otherwise throws an error if the folder already exists
    (output.overwrite && isdir(run_path)) || mkdir(run_path)
    return output.run_path = run_path
end

"""$(TYPEDSIGNATURES)
Concatenate the run folder name from `prefix`, `id` and `number` to
e.g. "run_0001" or "run_shallow_convection_0001"."""
function run_folder_name(prefix::String, id::String, number::Int; fmt = Printf.Format("%04d"))
    number_fmt = Printf.format(fmt, number)
    length(id) == 0 && return join((prefix, number_fmt), "_")   # otherwise "run__0001" with double underscore
    return join((prefix, id, number_fmt), "_")
end

"""$(TYPEDSIGNATURES)
Returns the full path of the output file after it was created."""
get_full_output_file_path(output::AbstractOutput) = joinpath(output.run_path, output.filename)

"""$(TYPEDSIGNATURES)
Returns the full path of the output file for a `simulation`. Throws an error if output is not active."""
get_output_path(simulation::AbstractSimulation) = get_output_path(simulation.model)
get_output_path(model::AbstractModel) = get_output_path(model.output)
function get_output_path(output::AbstractOutput)
    output.active || error("Output is not active")
    return joinpath(output.run_path, output.filename)
end

"""$(TYPEDSIGNATURES)
Loads a `var_name` trajectory of the model `M` that has been saved in
a netCDF file during the time stepping."""
function load_trajectory(var_name::Union{Symbol, String}, model::AbstractModel)
    @assert model.output.active "Output is turned off"
    return Array(NCDataset(get_full_output_file_path(model.output))[string(var_name)])
end

"""
$(TYPEDSIGNATURES)
Returns the output time step of the model `M`."""
function get_interval(output::AbstractOutput)
    return output.interval
end

# Fallback for when output is nothing
get_interval(::Nothing) = Millisecond(0)
