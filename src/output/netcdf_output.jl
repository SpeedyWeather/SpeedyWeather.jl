abstract type AbstractOutput end
abstract type AbstractOutputVariable end

# default number format for output
const DEFAULT_OUTPUT_NF = Float32
const DEFAULT_OUTPUT_DT = Hour(6)
const OUTPUT_VARIABLES_DICT = Dict{Symbol, AbstractOutputVariable}
OutputVariablesDict() = OUTPUT_VARIABLES_DICT()

# default compression
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
    pkg_version::VersionNumber = pkgversion(SpeedyWeather)

    # WHAT/WHEN OPTIONS
    startdate::DateTime = DateTime(2000, 1, 1)

    "[OPTION] output frequency, time step"
    output_dt::Second = Second(DEFAULT_OUTPUT_DT)

    "[OPTION] which variables to output, u, v, vor, div, pres, temp, humid"
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

function NetCDFOutput(
    S::SpectralGrid,
    Model::Type{<:AbstractModel} = Barotropic;
    output_Grid::Type{<:AbstractFullGrid} = RingGrids.full_grid_type(S.Grid),
    nlat_half::Integer = S.nlat_half, 
    output_NF::DataType = DEFAULT_OUTPUT_NF,
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
    println(io, "├ interpolator: $(typeof(output.interpolator))")
    println(io, "├ path: $(joinpath(output.run_path, output.filename))")
    println(io, "├ frequency: $(output.output_dt)")
    println(io, "└┐ variables:")
    for (key, var) in output.variables
        println(io, " ├ $key: $(var.long_name) [$(var.unit)]")
    end
end

function add!(D::OUTPUT_VARIABLES_DICT, outputvariables::AbstractOutputVariable...)
    for outputvariable in outputvariables
        key = Symbol(outputvariable.name)   # use name as key::Symbol
        D[key] = outputvariable
    end
end

add!(output::NetCDFOutput, outputvariables::AbstractOutputVariable...) =
    add!(output.variables, outputvariables...)
add!(model::AbstractModel, outputvariables::AbstractOutputVariable...) =
    add!(model.output, outputvariables...)

function add_default!(
    output_variables::OUTPUT_VARIABLES_DICT,
    Model::Type{<:Barotropic},
)
    add!(output_variables, VorticityOutput(), ZonalVelocityOutput(), MeridionalVelocityOutput())
end

function add_default!(
    variables::Dict{Symbol, AbstractOutputVariable},
    Model::Type{<:ShallowWater},
)
    add_default!(variables, Barotropic)
    add!(variables, InterfaceDisplacementOutput())
end

function add_default!(
    variables::Dict{Symbol, AbstractOutputVariable},
    Model::Type{<:PrimitiveDry},
)
    add_default!(variables, Barotropic)
    add!(variables, SurfacePressureOutput(), TemperatureOutput())
end

function add_default!(
    variables::Dict{Symbol, AbstractOutputVariable},
    Model::Type{<:PrimitiveWet},
)
    add_default!(variables, PrimitiveDry)
    add!(variables, HumidityOutput())
end

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
    all_dims = ["lon", "lat", "layer", "time"]

    # VARIABLES, define every output variable in the netCDF file and write initial conditions
    output_NF = eltype(output.grid2D)
    for (key, var) in output.variables
        attributes = Dict("long_name"=>var.long_name, "units"=>var.unit, "_FillValue"=>output_NF(var.missing_value))
        dims = collect(dim for (dim, this_dim) in zip(all_dims, var.dims_xyzt) if this_dim)

        # pick defaults if compression 
        deflatelevel = hasfield(typeof(var), :compression_level) ? var.compression_level : DEFAULT_COMPRESSION_LEVEL
        shuffle = hasfield(typeof(var), :shuffle) ? var.shuffle : DEFAULT_SHUFFLE

        defVar(dataset, var.name, output_NF, dims, attrib=attributes; deflatelevel, shuffle)
        output!(output, var, progn, diagn, model)
    end

    # # also export parameters into run????/parameters.txt
    # parameters_txt = open(joinpath(output.run_path, "parameters.txt"), "w")
    # println(parameters_txt, model.spectral_grid)
    # println(parameters_txt, model.planet)
    # println(parameters_txt, model.atmosphere)
    # println(parameters_txt, model.time_stepping)
    # println(parameters_txt, model.output)
    # println(parameters_txt, model.initial_conditions)
    # println(parameters_txt, model.horizontal_diffusion)
    # model isa Union{ShallowWater, PrimitiveEquation} && println(parameters_txt, model.implicit)
    # model isa Union{ShallowWater, PrimitiveEquation} && println(parameters_txt, model.orography)
    # close(parameters_txt)
end

Base.close(output::NetCDFOutput) = NCDatasets.close(output.netcdf_file)
Base.close(::Nothing) = nothing

"""
$(TYPEDSIGNATURES)
Writes the variables from `diagn` of time step `i` at time `time` into `outputter.netcdf_file`.
Simply escapes for no netcdf output of if output shouldn't be written on this time step.
Interpolates onto output grid and resolution as specified in `outputter`, converts to output
number format, truncates the mantissa for higher compression and applies lossless compression."""
function output!(
    output::NetCDFOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    output.timestep_counter += 1                                 # increase counter
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

## VORTICITY -------------

@kwdef mutable struct VorticityOutput <: AbstractOutputVariable
    name::String = "vor"
    unit::String = "s^-1"
    long_name::String = "relative vorticity"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 5
end

function output!(
    output::NetCDFOutput,
    variable::VorticityOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    vor = output.grid3D
    (; vor_grid) = diagn.grid
    RingGrids.interpolate!(vor, vor_grid, output.interpolator)

    unscale!(vor, diagn.scale[])    # was vor*radius, back to vor
    round!(vor, variable.keepbits)
    i = output.output_counter
    output.netcdf_file[variable.name][:, :, :, i] = vor
    return nothing
end

## U velocity -------------

@kwdef mutable struct ZonalVelocityOutput <: AbstractOutputVariable
    name::String = "u"
    unit::String = "m/s"
    long_name::String = "zonal wind"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

function output!(
    output::NetCDFOutput,
    variable::ZonalVelocityOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    u = output.grid3D
    (; u_grid) = diagn.grid
    RingGrids.interpolate!(u, u_grid, output.interpolator)

    round!(u, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, :, i] = u
    return nothing
end

## V velocity -------------

@kwdef mutable struct MeridionalVelocityOutput <: AbstractOutputVariable
    name::String = "v"
    unit::String = "m/s"
    long_name::String = "meridional wind"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

function output!(
    output::NetCDFOutput,
    variable::MeridionalVelocityOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    v = output.grid3D
    (; v_grid) = diagn.grid
    RingGrids.interpolate!(v, v_grid, output.interpolator)

    round!(v, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, :, i] = v
    return nothing
end

## DIVERGENCE -------------

@kwdef mutable struct DivergenceOutput <: AbstractOutputVariable
    name::String = "div"
    unit::String = "s^-1"
    long_name::String = "divergence"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 5
end

function output!(
    output::NetCDFOutput,
    variable::DivergenceOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    div = output.grid3D
    (; div_grid) = diagn.grid
    RingGrids.interpolate!(div, div_grid, output.interpolator)

    unscale!(div, diagn.scale[])    # was vor*radius, back to vor
    round!(div, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, :, i] = div
    return nothing
end

## INTERFACE DISPLACEMENT -------------

@kwdef mutable struct InterfaceDisplacementOutput <: AbstractOutputVariable
    name::String = "eta"
    unit::String = "m"
    long_name::String = "interface displacement"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

function output!(
    output::NetCDFOutput,
    variable::InterfaceDisplacementOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    eta = output.grid2D
    (; pres_grid) = diagn.grid
    RingGrids.interpolate!(eta, pres_grid, output.interpolator)

    round!(eta, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, i] = eta
    return nothing
end

## SURFACE PRESSURE -------------

@kwdef mutable struct SurfacePressureOutput <: AbstractOutputVariable
    name::String = "pres"
    unit::String = "hPa"
    long_name::String = "surface pressure"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 12
end

function output!(
    output::NetCDFOutput,
    variable::SurfacePressureOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    pres = output.grid2D
    (; pres_grid) = diagn.grid
    RingGrids.interpolate!(pres, pres_grid, output.interpolator)

    @inbounds for ij in eachindex(pres)
        pres[ij] = exp(pres[ij]) / 100    # from log(Pa) to hPa
    end

    round!(pres, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, i] = pres
    return nothing
end

## TEMPERATURE -------------

@kwdef mutable struct TemperatureOutput <: AbstractOutputVariable
    name::String = "temp"
    unit::String = "degC"
    long_name::String = "temperature"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 10
end

function output!(
    output::NetCDFOutput,
    variable::TemperatureOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    temp = output.grid3D
    (; temp_grid) = diagn.grid
    RingGrids.interpolate!(temp, temp_grid, output.interpolator)
    temp .-= 273.15             # convert from K to ˚C

    round!(temp, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, :, i] = temp
    return nothing
end

## HUMIDITY -------------

@kwdef mutable struct HumidityOutput <: AbstractOutputVariable
    name::String = "humid"
    unit::String = "kg/kg"
    long_name::String = "specific humidity"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

function output!(
    output::NetCDFOutput,
    variable::HumidityOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    humid = output.grid3D
    (; humid_grid) = diagn.grid
    RingGrids.interpolate!(humid, humid_grid, output.interpolator)

    round!(humid, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, :, i] = humid
    return nothing
end

## OROGRAPHY -------------

@kwdef mutable struct OrographyOutput <: AbstractOutputVariable
    name::String = "orography"
    unit::String = "m"
    long_name::String = "orography"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, false)
    missing_value::Float64 = NaN
    compression_level::Int = DEFAULT_COMPRESSION_LEVEL
    shuffle::Bool = DEFAULT_SHUFFLE
    keepbits::Int = Base.significand_bits(NF)
end

function output!(
    output::NetCDFOutput,
    variable::OrographyOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    # escape immediately when initialization counter > 1 to not write orography again
    output.output_counter > 1 || return nothing

    orog = output.grid2D
    (; orography) = model.orography
    RingGrids.interpolate!(orog, orography, output.interpolator)

    round!(orog, variable.keepbits)
    output.netcdf_file[variable.name][:, :] = orog
    return nothing
end

## PRECIPITATION -------------

@kwdef mutable struct ConvectivePrecipitationOutput <: AbstractOutputVariable
    name::String = "precip_conv"
    unit::String = "mm/hr"
    long_name::String = "convective precipitation"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

function output!(
    output::NetCDFOutput,
    variable::ConvectivePrecipitationOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    precip = output.grid2D
    (; precip_convection) = diagn.grid.physics
    RingGrids.interpolate!(precip, precip_convection, output.interpolator)
    
    # after output set precip accumulator back to zero
    precip_convection .= 0

    # convert from [m] to [mm/hr] rain rate over output time step (e.g. 6hours)
    s = (1000*Hour(1)/output.output_dt)
    precip .*= s

    round!(precip, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, i] = precip
    return nothing
end

@kwdef mutable struct LargeScalePrecipitationOutput <: AbstractOutputVariable
    name::String = "precip_cond"
    unit::String = "mm/hr"
    long_name::String = "large-scale precipitation"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

function output!(
    output::NetCDFOutput,
    variable::LargeScalePrecipitationOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    precip = output.grid2D
    (; precip_large_scale) = diagn.grid.physics
    RingGrids.interpolate!(precip, precip_large_scale, output.interpolator)

    # after output set precip accumulator back to zero
    precip_large_scale .= 0

    # convert from [m] to [mm/hr] rain rate over output time step (e.g. 6hours)
    s = (1000*Hour(1)/output.output_dt)
    precip .*= s

    round!(precip, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, i] = precip
    return nothing
end

## CLOUDS -------------

@kwdef mutable struct CloudTopOutput <: AbstractOutputVariable
    name::String = "cloud_top"
    unit::String = "m"
    long_name::String = "cloud top height"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

function output!(
    output::NetCDFOutput,
    variable::CloudTopOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    cloud = output.grid2D
    (; cloud_top) = diagn.grid.physics
    RingGrids.interpolate!(cloud, cloud_top, output.interpolator)

    round!(cloud, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, i] = cloud
    return nothing
end

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