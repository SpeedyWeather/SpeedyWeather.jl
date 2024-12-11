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
    Grid2D,
    Grid3D,
    Interpolator,
} <: AbstractOutput

    # FILE OPTIONS
    active::Bool = false
    
    "[OPTION] path to output folder, run_???? will be created within"
    path::String = pwd()
    
    "[OPTION] run identification number/string"
    id::String = "0001"
    run_path::String = ""                   # will be determined in initalize!
    
    "[OPTION] name of the output netcdf file"
    filename::String = "output.nc"
    
    "[OPTION] also write restart file if output==true?"
    write_restart::Bool = true
    pkg_version::VersionNumber = isnothing(pkgversion(SpeedyWeather)) ? v"0.0.0" : pkgversion(SpeedyWeather)

    # WHAT/WHEN OPTIONS
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

    # SCRATCH GRIDS TO INTERPOLATE ONTO
    const grid2D::Grid2D
    const grid3D::Grid3D
end

"""
$(TYPEDSIGNATURES)
Constructor for NetCDFOutput based on `S::SpectralGrid` and optionally
the `Model` type (e.g. `ShallowWater`, `PrimitiveWet`) as second positional argument.
The output grid is optionally determined by keyword arguments `output_Grid` (its type, full grid required),
`nlat_half` (resolution) and `output_NF` (number format). By default, uses the full grid
equivalent of the grid and resolution used in `SpectralGrid` `S`."""
function NetCDFOutput(
    S::SpectralGrid,
    Model::Type{<:AbstractModel} = Barotropic;
    output_Grid::Type{<:AbstractFullGrid} = RingGrids.full_grid_type(S.Grid),
    nlat_half::Integer = S.nlat_half, 
    output_NF::DataType = DEFAULT_OUTPUT_NF,
    output_dt::Period = Second(DEFAULT_OUTPUT_DT),  # only needed for dispatch
    kwargs...)

    # INPUT GRID
    input_Grid = S.Grid
    input_nlat_half = S.nlat_half

    # OUTPUT GRID
    nlon = RingGrids.get_nlon(output_Grid, nlat_half)
    nlat = RingGrids.get_nlat(output_Grid, nlat_half)
    npoints = nlon*nlat
    (; nlayers) = S

    # CREATE INTERPOLATOR
    interpolator = DEFAULT_INTERPOLATOR(DEFAULT_OUTPUT_NF, input_Grid, input_nlat_half, npoints)

    # CREATE GRIDS TO 
    output_Grid2D = RingGrids.nonparametric_type(output_Grid){output_NF, 1}
    output_Grid3D = RingGrids.nonparametric_type(output_Grid){output_NF, 2}
    grid2D = output_Grid2D(undef, nlat_half)
    grid3D = output_Grid3D(undef, nlat_half, nlayers)

    output = NetCDFOutput(;
        output_dt=Second(output_dt),    # convert to seconds for dispatch
        interpolator,
        grid2D,
        grid3D,
        kwargs...)

    add_default!(output.variables, Model)
    return output
end

function Base.show(io::IO, output::NetCDFOutput{Grid}) where Grid
    println(io, "NetCDFOutput{$Grid}")
    println(io, "├ status: $(output.active ? "active" : "inactive/uninitialized")")
    println(io, "├ write restart file: $(output.write_restart) (if active)")
    println(io, "├ interpolator: $(typeof(output.interpolator))")
    println(io, "├ path: $(joinpath(output.run_path, output.filename))")
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
    return nothing
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
    output::NetCDFOutput{Grid2D, Grid3D, Interpolator},
    feedback::AbstractFeedback,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
) where {Grid2D, Grid3D, Interpolator}
    
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
    Grid = typeof(output.grid2D)
    nlat_half = output.grid2D.nlat_half
    lond = get_lond(Grid, nlat_half)
    latd = get_latd(Grid, nlat_half)

    # INTERPOLATION: PRECOMPUTE LOCATION INDICES
    latds, londs = RingGrids.get_latdlonds(Grid, nlat_half)
    RingGrids.update_locator!(output.interpolator, latds, londs)
        
    σ = model.geometry.σ_levels_full
    defVar(dataset, "lon", lond, ("lon",), attrib=Dict("units"=>"degrees_east", "long_name"=>"longitude"))
    defVar(dataset, "lat", latd, ("lat",), attrib=Dict("units"=>"degrees_north", "long_name"=>"latitude"))
    defVar(dataset, "layer", σ, ("layer",), attrib=Dict("units"=>"1", "long_name"=>"sigma layer"))

    # VARIABLES, define every output variable in the netCDF file and write initial conditions
    output_NF = eltype(output.grid2D)
    for (key, var) in output.variables
        define_variable!(dataset, var, output_NF)
        output!(output, var, progn, diagn, model)
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

    all_dims = ("lon", "lat", "layer", "time")
    dims = collect(dim for (dim, this_dim) in zip(all_dims, var.dims_xyzt) if this_dim)

    # pick defaults for compression if not defined
    deflatelevel = hasfield(typeof(var), :compression_level) ? var.compression_level : DEFAULT_COMPRESSION_LEVEL
    shuffle = hasfield(typeof(var), :shuffle) ? var.shuffle : DEFAULT_SHUFFLE

    defVar(dataset, var.name, output_NF, dims, attrib=attributes; deflatelevel, shuffle)
end

"""
$(TYPEDSIGNATURES)
Writes the variables from `progn` or `diagn` of time step `i` at time `time` into `output.netcdf_file`.
Simply escapes for no netcdf output or if output shouldn't be written on this time step.
Interpolates onto output grid and resolution as specified in `output`, converts to output
number format, truncates the mantissa for higher compression and applies lossless compression."""
function output!(
    output::NetCDFOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    output.timestep_counter += 1                                    # increase counter
    (; active, output_every_n_steps, timestep_counter ) = output
    active || return nothing                                        # escape immediately for no netcdf output
    timestep_counter % output_every_n_steps == 0 || return nothing  # escape if output not written on this step

    output!(output, progn.clock.time)                               # increase counter write time
    output!(output, output.variables, progn, diagn, model)          # write variables
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
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    for (key, var) in output_variables
        output!(output, var, progn, diagn, model)
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
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    if output.active    # only finalize if active otherwise output.netcdf_file is nothing
        for (key, var) in output.variables
            finalize!(output, var, progn, diagn, model)
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

# define all variables for output
include("variables/dynamics.jl")
include("variables/precipitation.jl")   # collected as PrecipitationOutput()
include("variables/boundaries.jl")      # BoundaryOutput()
include("variables/radiation.jl")       # RadiationOutput()
include("variables/stochastic.jl")      # RandomPatternOutput()
include("variables/surface_fluxes.jl")  # SurfaceFluxesOutput()
include("variables/ocean.jl")           # SeaSurfaceTemperatureOutput()

# collect all together for conveneince
AllOutputVariables() = (
    PrecipitationOutput()...,
    BoundaryOutput()...,
    RadiationOutput()...,
    RandomPatternOutput(),
    SurfaceFluxesOutput()...,
    SeaSurfaceTemperatureOutput(),
)

"""
$(TYPEDSIGNATURES)
Checks existing `run_????` folders in `path` to determine a 4-digit `id` number
by counting up. E.g. if folder run_0001 exists it will return the string "0002".
Does not create a folder for the returned run id.
"""
function get_run_id(path::String, id::String)
    # if run_???? folder doesn't exist yet don't change the id
    run_id = string("run_", run_id_to_string(id))
    !(run_id in readdir(path)) && return id

    # otherwise pull list of existing run_???? folders via readdir
    pattern = r"run_\d\d\d\d"               # run_???? in regex
    runlist = filter(x->startswith(x, pattern), readdir(path))
    runlist = filter(x->endswith(  x, pattern), runlist)
    existing_runs = [parse(Int, id[5:end]) for id in runlist]

    # get the run id from existing folders
    if length(existing_runs) == 0           # if no runfolder exists yet
        run_id = 1                          # start with run_0001
    else
        run_id = maximum(existing_runs)+1   # next run gets id +1
    end
    
    return @sprintf("%04d", run_id)
end

"""
$(TYPEDSIGNATURES)
Creates a new folder `run_*` with the identification `id`. Also returns the full path
`run_path` of that folder.
"""
function create_output_folder(path::String, id::Union{String, Int})
    run_id = string("run_", run_id_to_string(id))
    run_path = joinpath(path, run_id)
    @assert !(run_id in readdir(path)) "Run folder $run_path already exists."
    mkdir(run_path)             # actually create the folder
    return run_path
end

run_id_to_string(run_id::Integer) = @sprintf("%04d", run_id)
run_id_to_string(run_id::String) = run_id

"""
$(TYPEDSIGNATURES)
Returns the full path of the output file after it was created.
"""
get_full_output_file_path(output::AbstractOutput) = joinpath(output.run_path, output.filename)

"""
$(TYPEDSIGNATURES)
Loads a `var_name` trajectory of the model `M` that has been saved in a netCDF file during the time stepping.
"""
function load_trajectory(var_name::Union{Symbol, String}, model::AbstractModel) 
    @assert model.output.active "Output is turned off"
    return Array(NCDataset(get_full_output_file_path(model.output))[string(var_name)])
end