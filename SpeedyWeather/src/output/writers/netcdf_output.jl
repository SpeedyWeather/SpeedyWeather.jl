export NetCDFOutput

"""Output writer for a netCDF file with (re-)gridded variables.
Interpolates non-rectangular grids. Fields are
$(TYPEDFIELDS)"""
@kwdef mutable struct NetCDFOutput{
        Field2D,
        Field3D,
        Interpolator,
        DT,
        S,
    } <: AbstractOutput

    # FILE OPTIONS
    active::Bool = false

    "[OPTION] path to output parent folder, run folders will be created within"
    path::String = pwd()

    "[OPTION] Prefix for run folder where data is stored, e.g. 'run_'"
    run_prefix::String = "run"

    "[OPTION] run identification, added between run_prefix and run_number"
    id::String = ""

    "[OPTION] run identification number, automatically determined if overwrite=false"
    run_number::Int = 1

    "[OPTION] run numbers digits"
    run_digits::Int = 4

    "[DERIVED] shared output writer state (run folder, counters, output frequency)"
    core::OutputWriterCore = OutputWriterCore()

    "[OPTION] Overwrite an existing run folder?"
    overwrite::Bool = false

    "[OPTION] name of the output netcdf file"
    filename::String = "output.nc"

    "[OPTION] also write restart file if output=true?"
    write_restart::Bool = true

    "[OPTION] also write parameters txt file if output=true?"
    write_parameters_txt::Bool = true

    "[OPTION] also write progress txt file if output=true?"
    write_progress_txt::Bool = true

    # WHAT/WHEN OPTIONS
    "[DERIVD] start date of the simulation, used for time dimension in netcdf file"
    startdate::DT = DateTime(2000, 1, 1)

    "[OPTION] output frequency, time step"
    interval::S = Second(DEFAULT_OUTPUT_INTERVAL)

    "[OPTION] dictionary of variables to output, e.g. u, v, vor, div, pres, temp, humid"
    variables::OUTPUT_VARIABLES_DICT = OutputVariablesDict()

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
the `Model` type (e.g. `ShallowWater`, `PrimitiveWet`) as second positional argument. In case a 
non-default number of soil layers is used, it also needs the respective `nlayers_soil` to allocate those outputs.
The output grid is optionally determined by keyword arguments `output_Grid` (its type, full grid required),
`nlat_half` (resolution) and `output_NF` (number format, only used for variables not coordinates).
By default, uses the full grid equivalent of the grid and resolution used in `SpectralGrid` `S`."""
function NetCDFOutput(
        SG::SpectralGrid,
        Model::Type{<:AbstractModel} = Barotropic;
        nlayers_soil = DEFAULT_NLAYERS_SOIL,
        output_grid::AbstractFullGrid = on_architecture(CPU(), RingGrids.full_grid_type(SG.grid)(SG.grid.nlat_half)),
        output_NF::DataType = DEFAULT_OUTPUT_NF,
        interval::Period = Second(DEFAULT_OUTPUT_INTERVAL),  # only needed for dispatch
        kwargs...
    )

    # INPUT GRID (but on CPU)
    input_grid = on_architecture(CPU(), SG.grid)

    # CREATE INTERPOLATOR
    interpolator = RingGrids.interpolator(output_grid, input_grid, NF = DEFAULT_OUTPUT_NF)

    # CREATE FULL FIELDS TO INTERPOLATE ONTO BEFORE WRITING DATA OUT
    (; nlayers) = SG

    field2D = Field(output_NF, output_grid)
    field3D = Field(output_NF, output_grid, nlayers)
    field3Dland = Field(output_NF, output_grid, nlayers_soil)

    output = NetCDFOutput(;
        interval = Second(interval),    # convert to seconds for dispatch
        interpolator,
        field2D,
        field3D,
        field3Dland,
        kwargs...
    )

    add_default!(output.variables, Model)
    return output
end

function Base.show(io::IO, output::NetCDFOutput{F}) where {F}

    F_str = string("{", F, "}")
    type_param_str = length(F_str) > 30 ? string(first(F_str, 30), "...}") : F_str
    active = output.active ? "active" : "inactive/uninitialized"

    println(io, styled"{warning:NetCDFOutput}{note:$type_param_str}")
    println(io, styled"├ {info:status}: $active")
    println(io, styled"├ {info:write restart file} = $(output.write_restart) (if active)")

    interp_type_str = string(typeof(output.interpolator))
    interp_type_str_short = length(interp_type_str) > 70 ? string(first(interp_type_str, 70), "...}") : interp_type_str

    println(io, styled"├ {info:interpolator}::$interp_type_str_short")
    println(io, styled"├ {info:path} = $(joinpath(output.run_path, output.filename)) (overwrite=$(output.overwrite))")
    println(io, styled"├ {info:interval} = $(output.interval)")
    print(io, styled"└ {info:variables}")
    nvars = length(output.variables)
    for (i, (key, var)) in enumerate(output.variables)
        print(io, "\n  $(i == nvars ? "└" : "├") ", styled"{magenta:$key}: $(var.long_name) ", styled"{note:[$(var.unit)]}")
    end
    return nothing
end

"""$(TYPEDSIGNATURES)
Initialize NetCDF `output` by creating a netCDF file and storing the initial conditions
of `vars`. To be called just before the first timesteps."""
function initialize!(
        output::NetCDFOutput,
        vars::Variables,
        model::AbstractModel,
    )
    output.active || return nothing             # exit immediately for no output

    # only checked for models that have a land component and output variables that
    # actually use the soil vertical dimension (and its scratch field `field3Dland`)
    if hasfield(typeof(model), :land) && !isnothing(model.land) &&
            any(var -> is_land(var) && is3D(var), values(output.variables))
        @assert get_nlayers(model.land) == size(output.field3Dland, 2) "$(size(output.field3Dland, 2))" *
            " soil layers initialized for output, but $(get_nlayers(model.land)) soil layers initialized for model." *
            " Please construct NetCDFOutput with the same `nlayers_soil` as the model."
    end

    # SHARED INITIALIZATION (run folder, output frequency, counters, callbacks)
    initialize!(output.core, output, model)

    # CREATE NETCDF FILE, vector of NcVars for output
    (; run_path, filename) = output
    dataset = NCDataset(joinpath(run_path, filename), "c")
    output.netcdf_file = dataset

    # DEFINE NETCDF DIMENSIONS TIME and write current (=initial) time
    (; startdate) = output
    time_string = "hours since $(Dates.format(startdate, "yyyy-mm-dd HH:MM:0.0"))"
    defDim(dataset, "time", Inf)        # unlimited time dimension
    defVar(
        dataset, "time", Float64, ("time",),
        attrib = Dict(
            "units" => time_string, "long_name" => "time",
            "standard_name" => "time", "calendar" => "proleptic_gregorian"
        )
    )
    output!(output, vars.prognostic.clock.time)   # write initial time

    # DEFINE NETCDF DIMENSIONS SPACE
    # explictly move to CPU and convert to common format as determined by RingGrids.get_lond (Float64 default)
    lond = get_lond(output.field2D)
    latd = get_latd(output.field2D)
    σ = convert.(eltype(lond), on_architecture(CPU(), model.geometry.σ_levels_full))
    soil_indices = collect(1:get_soil_layers(model))

    defVar(dataset, "lon", lond, ("lon",), attrib = Dict("units" => "degrees_east", "long_name" => "longitude"))
    defVar(dataset, "lat", latd, ("lat",), attrib = Dict("units" => "degrees_north", "long_name" => "latitude"))
    defVar(dataset, "layer", σ, ("layer",), attrib = Dict("units" => "1", "long_name" => "sigma layer"))
    defVar(dataset, "soil_layer", soil_indices, ("soil_layer",), attrib = Dict("units" => "1", "long_name" => "soil layer index"))

    # VARIABLES, remove output variables not existent in simulation.variables
    simulation = Simulation(vars, model)
    nonexisting_vars = [key for (key, var) in output.variables if isnothing(path_or_nothing(var, simulation))]
    if !isempty(nonexisting_vars)
        @warn "Some output.variables do not exist in simulation. Deleting: $(join(nonexisting_vars, ", "))"
    end
    delete!(output, nonexisting_vars...)

    # then define every output variable in the netCDF file and write initial conditions
    for (key, var) in output.variables
        define_variable!(dataset, var, eltype(output.field2D))
        output!(output, var, simulation)
    end

    return nothing
end

Base.close(output::NetCDFOutput) = NCDatasets.close(output.netcdf_file)

function define_variable!(
        dataset::NCDataset,
        var::AbstractOutputVariable,
        output_NF::Type{<:AbstractFloat} = DEFAULT_OUTPUT_NF,
    )
    # hook for custom output variables to lazily define their own (vertical) dimension
    define_dimension!(dataset, var)

    missing_value = hasfield(typeof(var), :missing_value) ? var.missing_value : DEFAULT_MISSING_VALUE
    attributes = Dict("long_name" => var.long_name, "units" => var.unit, "_FillValue" => output_NF(missing_value))

    # the vertical dimension depends on the variable, e.g. "layer" or "soil_layer"
    all_dims = ("lon", "lat", vertical_dimension(var), "time")
    dims = collect(dim for (dim, this_dim) in zip(all_dims, var.dims_xyzt) if this_dim)

    # pick defaults for compression if not defined
    deflatelevel = hasproperty(var, :compression_level) ? var.compression_level : DEFAULT_COMPRESSION_LEVEL
    shuffle = hasproperty(var, :shuffle) ? var.shuffle : DEFAULT_SHUFFLE

    return defVar(dataset, var.name, output_NF, dims, attrib = attributes; deflatelevel, shuffle)
end

"""$(TYPEDSIGNATURES)
Length of dimension `name` in `dataset` or `nothing` if not defined.
Used by custom output variables to lazily define their own dimension in
[`define_dimension!`](@ref); a Zarr store equivalent is defined in the
Zarr extension."""
get_dimension(dataset::NCDataset, name::String) =
    haskey(dataset.dim, name) ? dataset.dim[name] : nothing

"""$(TYPEDSIGNATURES)
Define a coordinate in `dataset`: a dimension `name` of `length(values)` with
`values` as its coordinate variable and `attribs` as its attributes.
Used by custom output variables to lazily define their own dimension in
[`define_dimension!`](@ref); a Zarr store equivalent is defined in the
Zarr extension."""
define_coordinate!(dataset::NCDataset, name::String, values::AbstractVector; attribs = Dict{String, String}()) =
    defVar(dataset, name, values, (name,), attrib = attribs)

"""$(TYPEDSIGNATURES)
Backend-specific write of an interpolated, post-processed field `field` for `variable`
into `output`. Implementations exist for [`NetCDFOutput`](@ref) and `ZarrOutput`."""
function write_array!(
        output::NetCDFOutput,
        variable::AbstractOutputVariable,
        field,
    )
    i = output.output_counter
    indices = get_indices(i, variable)
    output.netcdf_file[variable.name][indices...] = field
    return nothing
end

"""
$(TYPEDSIGNATURES)
Write the current time `time::DateTime` to the netCDF file in `output`."""
function output!(
        output::NetCDFOutput,
        time::DateTime,
    )
    i = output.output_counter

    (; netcdf_file, startdate) = output
    time_passed = Millisecond(time - startdate)
    time_hrs = time_passed.value / 3600_000       # [ms] -> [hrs]
    netcdf_file["time"][i] = time_hrs
    NCDatasets.sync(netcdf_file)

    return nothing
end
