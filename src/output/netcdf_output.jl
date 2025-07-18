abstract type AbstractOutput <: AbstractModelComponent end
abstract type AbstractOutputVariable end

# default number format for output
const DEFAULT_OUTPUT_NF = Float32
const DEFAULT_OUTPUT_DT = Hour(6)
const OUTPUT_VARIABLES_DICT = Dict{Symbol, AbstractOutputVariable}
OutputVariablesDict() = OUTPUT_VARIABLES_DICT()

const DEFAULT_MISSING_VALUE = NaN
const DEFAULT_COMPRESSION_LEVEL = 1
const DEFAULT_SHUFFLE = false
const DEFAULT_KEEPBITS = 15

export NetCDFOutput

"""Output writer for a netCDF file with (re-)gridded variables.
Interpolates non-rectangular grids. Fields are
$(TYPEDFIELDS)"""
@kwdef mutable struct NetCDFOutput{
    Field2D,
    Field3D,
    Interpolator,
} <: AbstractOutput

    # FILE OPTIONS
    active::Bool = false

    "[OPTION] path to output parent folder, run folders will be created within"
    path::String = pwd()
    
    "[OPTION] Prefix for run folder where data is stored, e.g. 'run_'"
    run_prefix::String = "run"

    "[OPTION] run identification, added between run_prefix and run_number"
    id::String = ""

    "[OPTION] run identification number"
    run_number::Int = 1

    "[OPTION] run numbers digits"
    run_digits::Int = 4

    "[DERIVED] folder name where data is stored, determined at initialize!"
    run_folder::String = ""

    "[DERIVED] full path to folder where data is stored, determined at initialize!"
    run_path::String = ""

    "[OPTION] Overwrite an existing run folder?"
    overwrite::Bool = false
    
    "[OPTION] name of the output netcdf file"
    filename::String = "output.nc"
    
    "[OPTION] also write restart file if output==true?"
    write_restart::Bool = true

    "[DERIVED] package version, used for restart files"
    pkg_version::VersionNumber = isnothing(pkgversion(SpeedyWeather)) ? v"0.0.0" : pkgversion(SpeedyWeather)

    # WHAT/WHEN OPTIONS
    "[DERIVD] start date of the simulation, used for time dimension in netcdf file"
    startdate::DateTime = DateTime(2000, 1, 1)

    "[OPTION] output frequency, time step"
    output_dt::Second = Second(DEFAULT_OUTPUT_DT)

    "[OPTION] dictionary of variables to output, e.g. u, v, vor, div, pres, temp, humid"
    variables::OUTPUT_VARIABLES_DICT = OutputVariablesDict()

    # TIME STEPS AND COUNTERS (initialize later)
    output_every_n_steps::Int = 0           # output frequency
    timestep_counter::Int = 0               # time step counter
    output_counter::Int = 0                 # output step counter
    
    # the netcdf file to be written into, will be created
    netcdf_file::Union{NCDataset, Nothing} = nothing

    const interpolator::Interpolator

    # SCRATCH FIELDS TO INTERPOLATE ONTO
    const field2D::Field2D
    const field3D::Field3D
    const field3Dland::Field3D
end

"""
$(TYPEDSIGNATURES)
Constructor for NetCDFOutput based on `S::SpectralGrid` and optionally
the `Model` type (e.g. `ShallowWater`, `PrimitiveWet`) as second positional argument.
The output grid is optionally determined by keyword arguments `output_Grid` (its type, full grid required),
`nlat_half` (resolution) and `output_NF` (number format). By default, uses the full grid
equivalent of the grid and resolution used in `SpectralGrid` `S`."""
function NetCDFOutput(
    SG::SpectralGrid,
    Model::Type{<:AbstractModel} = Barotropic;
    output_grid::AbstractFullGrid = RingGrids.full_grid_type(SG.grid)(SG.grid.nlat_half),
    output_NF::DataType = DEFAULT_OUTPUT_NF,
    output_dt::Period = Second(DEFAULT_OUTPUT_DT),  # only needed for dispatch
    kwargs...)

    # INPUT GRID
    input_grid = SG.grid

    # CREATE INTERPOLATOR
    interpolator = RingGrids.interpolator(output_grid, input_grid, NF=DEFAULT_OUTPUT_NF)
    
    # CREATE FULL FIELDS TO INTERPOLATE ONTO BEFORE WRITING DATA OUT
    (; nlayers, nlayers_soil) = SG
    field2D = Field(output_NF, output_grid)
    field3D = Field(output_NF, output_grid, nlayers)
    field3Dland = Field(output_NF, output_grid, nlayers_soil)

    output = NetCDFOutput(;
        output_dt=Second(output_dt),    # convert to seconds for dispatch
        interpolator,
        field2D,
        field3D,
        field3Dland,
        kwargs...)

    add_default!(output.variables, Model)
    return output
end

function Base.show(io::IO, output::NetCDFOutput{F}) where F
    println(io, "NetCDFOutput{$F}")
    println(io, "├ status: $(output.active ? "active" : "inactive/uninitialized")")
    println(io, "├ write restart file: $(output.write_restart) (if active)")
    println(io, "├ interpolator: $(typeof(output.interpolator))")
    println(io, "├ path: $(joinpath(output.run_path, output.filename)) (overwrite=$(output.overwrite))")
    println(io, "├ frequency: $(output.output_dt)")
    print(io,   "└┐ variables:")
    nvars = length(output.variables)
    for (i, (key, var)) in enumerate(output.variables)
        print(io, "\n $(i==nvars ? "└" : "├") $key: $(var.long_name) [$(var.unit)]")
    end
end

"""
$(TYPEDSIGNATURES)
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
function add!(output::NetCDFOutput, outputvariables::AbstractOutputVariable...)
    add!(output.variables, outputvariables...)
    return output
end

"""$(TYPEDSIGNATURES)
Add `outputvariables` to the dictionary in `output::NetCDFOutput` of `model`, i.e. at `model.output.variables`."""
function add!(model::AbstractModel, outputvariables::AbstractOutputVariable...)
    add!(model.output, outputvariables...)
    return model.output
end

"""$(TYPEDSIGNATURES)
Delete output variables from `output` by their (short name) (Symbol or String), corresponding
to the keys in the dictionary."""
function Base.delete!(output::NetCDFOutput, keys::Union{String, Symbol}...)
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
    add!(output_variables, VorticityOutput(), ZonalVelocityOutput(), MeridionalVelocityOutput())
end

"""$(TYPEDSIGNATURES)
Add default variables to output for a `ShallowWater` model, same as for a `Barotropic` model but also
the interface displacement."""
function add_default!(
    variables::Dict{Symbol, AbstractOutputVariable},
    Model::Type{<:ShallowWater},
)
    add_default!(variables, Barotropic)
    add!(variables, InterfaceDisplacementOutput())
end

"""$(TYPEDSIGNATURES)
Add default variables to output for a `PrimitiveDry` model, same as for a `Barotropic` model but also
the surface pressure and temperature."""
function add_default!(
    variables::Dict{Symbol, AbstractOutputVariable},
    Model::Type{<:PrimitiveDry},
)
    add_default!(variables, Barotropic)
    add!(variables, SurfacePressureOutput(), TemperatureOutput())
end

"""$(TYPEDSIGNATURES)
Add default variables to output for a `PrimitiveWet` model, same as for a `PrimitiveDry` model but also
the specific humidity."""
function add_default!(
    variables::Dict{Symbol, AbstractOutputVariable},
    Model::Type{<:PrimitiveWet},
)
    add_default!(variables, PrimitiveDry)
    add!(variables, HumidityOutput())
end

"""$(TYPEDSIGNATURES)
Initialize NetCDF `output` by creating a netCDF file and storing the initial conditions
of `diagn` (and `progn`). To be called just before the first timesteps."""
function initialize!(   
    output::NetCDFOutput,
    feedback::AbstractFeedback,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    output.active || return nothing             # exit immediately for no output
    
    # GET RUN ID, CREATE FOLDER
    # get new id only if not already specified
    determine_run_folder!(output)
    create_run_folder!(output)

    feedback.run_folder = output.run_folder     # synchronize with feedback struct
    feedback.run_path = output.run_path
    feedback.progress_meter.desc = "Weather is speedy: $(output.run_folder) "
    feedback.output = true              # if output=true set feedback.output=true too!

    # OUTPUT FREQUENCY
    output.output_every_n_steps = max(1, round(Int,
            Millisecond(output.output_dt).value/model.time_stepping.Δt_millisec.value))
    output.output_dt = Second(round(Int, output.output_every_n_steps*model.time_stepping.Δt_sec))

    # RESET COUNTERS
    output.timestep_counter = 0         # time step counter
    output.output_counter = 0           # output step counter

    # CREATE NETCDF FILE, vector of NcVars for output
    (; run_path, filename) = output
    dataset = NCDataset(joinpath(run_path, filename), "c")
    output.netcdf_file = dataset
    
    # DEFINE NETCDF DIMENSIONS TIME and write current (=initial) time
    (; startdate) = output
    time_string = "hours since $(Dates.format(startdate, "yyyy-mm-dd HH:MM:0.0"))"
    defDim(dataset, "time", Inf)        # unlimited time dimension
    defVar(dataset, "time", Float64, ("time",),
           attrib=Dict("units"=>time_string, "long_name"=>"time",
                       "standard_name"=>"time", "calendar"=>"proleptic_gregorian"))
    output!(output, progn.clock.time)   # write initial time

    # DEFINE NETCDF DIMENSIONS SPACE
    lond = get_lond(output.field2D)
    latd = get_latd(output.field2D)
    σ = model.geometry.σ_levels_full
    soil_indices = collect(1:model.spectral_grid.nlayers_soil)
        
    defVar(dataset, "lon", lond, ("lon",), attrib=Dict("units"=>"degrees_east", "long_name"=>"longitude"))
    defVar(dataset, "lat", latd, ("lat",), attrib=Dict("units"=>"degrees_north", "long_name"=>"latitude"))
    defVar(dataset, "layer", σ, ("layer",), attrib=Dict("units"=>"1", "long_name"=>"sigma layer"))
    defVar(dataset, "soil_layer", soil_indices, ("soil_layer",), attrib=Dict("units"=>"1", "long_name"=>"soil layer index"))

    # VARIABLES, define every output variable in the netCDF file and write initial conditions
    output_NF = eltype(output.field2D)
    for (key, var) in output.variables
        define_variable!(dataset, var, output_NF)
        output!(output, var, Simulation(progn, diagn, model))
    end

    # also export parameters into run????/parameters.txt
    parameters_txt = open(joinpath(output.run_path, "parameters.txt"), "w")
    for property in propertynames(model)
        println(parameters_txt, "model.$property")
        println(parameters_txt, getfield(model, property,), "\n")
    end
    close(parameters_txt)
end

Base.close(output::NetCDFOutput) = NCDatasets.close(output.netcdf_file)
Base.close(::Nothing) = nothing     # in case of no netCDF output nothing to close

function define_variable!(
    dataset::NCDataset,
    var::AbstractOutputVariable,
    output_NF::Type{<:AbstractFloat} = DEFAULT_OUTPUT_NF,
)
    missing_value = hasfield(typeof(var), :missing_value) ? var.missing_value : DEFAULT_MISSING_VALUE
    attributes = Dict("long_name"=>var.long_name, "units"=>var.unit, "_FillValue"=>output_NF(missing_value))

    # land variables have a different vertical dimension
    all_dims = is_land(var) ? ("lon", "lat", "soil_layer", "time") : ("lon", "lat", "layer", "time")
    dims = collect(dim for (dim, this_dim) in zip(all_dims, var.dims_xyzt) if this_dim)

    # pick defaults for compression if not defined
    deflatelevel = hasproperty(var, :compression_level) ? var.compression_level : DEFAULT_COMPRESSION_LEVEL
    shuffle = hasproperty(var, :shuffle) ? var.shuffle : DEFAULT_SHUFFLE

    defVar(dataset, var.name, output_NF, dims, attrib=attributes; deflatelevel, shuffle)
end

"""
$(TYPEDSIGNATURES)
Writes the variables from `progn` or `diagn` of time step `i` at time `time` into `output.netcdf_file`.
Simply escapes for no netcdf output or if output shouldn't be written on this time step.
Interpolates onto output grid and resolution as specified in `output`, converts to output
number format, truncates the mantissa for higher compression and applies lossless compression."""
function output!(output::NetCDFOutput, simulation::AbstractSimulation)
    output.timestep_counter += 1                                        # increase counter
    (; active, output_every_n_steps, timestep_counter ) = output
    active || return nothing                                            # escape immediately for no netcdf output
    timestep_counter % output_every_n_steps == 0 || return nothing      # escape if output not written on this step

    (; clock) = simulation.prognostic_variables
    output!(output, clock.time)                                         # increase counter write time
    output!(output, output.variables, simulation)                       # write variables
end

get_indices(i, variable::AbstractOutputVariable) = get_indices(i, Val.(variable.dims_xyzt)...)
get_indices(i, x::Val{true}, y::Val{true}, z::Val{true}, t::Val{true}) = (:, :, :, i)   # 3D + time
get_indices(i, x::Val{true}, y::Val{true}, z::Val{true}, t::Val{false}) = (:, :, :)     # 3D
get_indices(i, x::Val{true}, y::Val{true}, z::Val{false}, t::Val{true}) = (:, :, i)     # 2D + time
get_indices(i, x::Val{true}, y::Val{true}, z::Val{false}, t::Val{false}) = (:, :)       # 2D

is3D(variable::AbstractOutputVariable) = variable.dims_xyzt[3]
is_land(variable::AbstractOutputVariable) = hasproperty(variable, :is_land) ? variable.is_land : false
hastime(variable::AbstractOutputVariable) = variable.dims_xyzt[4]

"""$(TYPEDSIGNATURES)
Output a `variable` into the netCDF file `output.netcdf_file`.
Interpolates onto the output grid and resolution as specified in `output`.
Method used for all output variables `<: AbstractOutputVariable`
with dispatch over the second argument. Interpolates, scales,
custom transform, bitrounding and writes to file."""
function output!(
    output::NetCDFOutput,
    variable::AbstractOutputVariable,
    simulation::AbstractSimulation,
)
    # escape immediately after first call if variable doesn't have a time dimension
    ~hastime(variable) && output.output_counter > 1 && return nothing

    # interpolate 2D/3D variables
    var = is3D(variable) ? (is_land(variable) ? output.field3Dland : output.field3D) : output.field2D
    raw = path(variable, simulation)
    RingGrids.interpolate!(var, raw, output.interpolator)

    # unscale if variable.unscale == true and exists
    if hasproperty(variable, :unscale)
        if variable.unscale
            unscale!(var, simulation.diagnostic_variables.scale[])
        end
    end

    if hasproperty(variable, :transform)    # transform (e.g. scale, offset, exp, etc) if defined
        @. var = variable.transform(var)
    end

    if hasproperty(variable, :keepbits)     # round mantissabits for compression
        round!(var, variable.keepbits)
    end

    i = output.output_counter               # output time step i to write
    indices = get_indices(i, variable)      # returns (:, :, i) for example, depending on dims
    output.netcdf_file[variable.name][indices...] = var     # actually write to file
    return nothing
end

"""
$(TYPEDSIGNATURES)
Write the current time `time::DateTime` to the netCDF file in `output`."""
function output!(
    output::NetCDFOutput,
    time::DateTime,
)
    output.output_counter += 1      # output counter increases when writing time
    i = output.output_counter

    (; netcdf_file, startdate ) = output
    time_passed = Millisecond(time-startdate)
    time_hrs = time_passed.value/3600_000       # [ms] -> [hrs]
    netcdf_file["time"][i] = time_hrs
    NCDatasets.sync(netcdf_file)

    return nothing
end

"""$(TYPEDSIGNATURES)
Loop over every variable in `output.variables` to call the respective `output!` method
to write into the `output.netcdf_file`."""
function output!(
    output::NetCDFOutput,
    output_variables::OUTPUT_VARIABLES_DICT,
    simulation::AbstractSimulation,
)
    for var in values(output_variables)
        output!(output, var, simulation)
    end
end

function Base.show(io::IO, outputvariable::AbstractOutputVariable)
    print(io, "$(typeof(outputvariable)) <: SpeedyWeather.AbstractOutputVariable")
    for field in propertynames(outputvariable)
        value = getfield(outputvariable, field)
        print(io, "\n├ $field::$(typeof(value)) = $value")
    end
end

function finalize!(
    output::NetCDFOutput,
    simulation::AbstractSimulation,
)
    if output.active    # only finalize if active otherwise output.netcdf_file is nothing
        for var in values(output.variables)
            finalize!(output, var, simulation)
        end
    end
    close(output)
end

# default finalize method for output variables
function finalize!(
    output::NetCDFOutput,
    var::AbstractOutputVariable,
    args...
)
    # do nothing unless finalize! is defined for var <: AbstractOutputVariable
    return nothing
end

"""Dummy output variable that doesn't do anything."""
struct NoOutputVariable <: AbstractOutputVariable end
output!(output::NetCDFOutput, variable::NoOutputVariable, args...) = nothing

include("variables/output_variables.jl")

"""$(TYPEDSIGNATURES)
Checks existing folders in `path` and determine `run_number`by counting up.
E.g. if folder run_0001 exists then run_number is 2.
Does not create a folder for the returned run id."""
function determine_run_folder!(output::NetCDFOutput)
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

    output.run_folder = run_folder_name(run_prefix, id, output.run_number; fmt)
end

"""$(TYPEDSIGNATURES)
Creates a new folder `prefix_id_number` with the identification `id`. Also returns the full path
`run_path` of that folder."""
function create_run_folder!(output::NetCDFOutput)
    run_path = joinpath(output.path, output.run_folder)

    # actually create the folder
    # unless overwrite is true and the folder exists
    # otherwise throws an error if the folder already exists
    (output.overwrite && isdir(run_path)) || mkdir(run_path)
    output.run_path = run_path
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
Loads a `var_name` trajectory of the model `M` that has been saved in
a netCDF file during the time stepping."""
function load_trajectory(var_name::Union{Symbol, String}, model::AbstractModel) 
    @assert model.output.active "Output is turned off"
    return Array(NCDataset(get_full_output_file_path(model.output))[string(var_name)])
end