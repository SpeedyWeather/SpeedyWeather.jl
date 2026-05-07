module SpeedyWeatherZarrExt

using SpeedyWeather
using Zarr
using DocStringExtensions

import SpeedyWeather: ZarrOutput, AbstractOutput, AbstractOutputVariable,
    AbstractSimulation, AbstractModel, Barotropic, OutputWriterCore,
    OUTPUT_VARIABLES_DICT, OutputVariablesDict, DEFAULT_NLAYERS_SOIL,
    DEFAULT_OUTPUT_NF, DEFAULT_OUTPUT_INTERVAL, DEFAULT_MISSING_VALUE,
    DEFAULT_COMPRESSION_LEVEL, DEFAULT_KEEPBITS,
    Variables, Simulation, SpectralGrid, Field,
    initialize!, finalize!, output!, write_array!, set!, add!, add_default!,
    is3D, is_land, hastime, get_indices, scale!, get_soil_layers,
    get_lond, get_latd, on_architecture, CPU,
    AbstractFullGrid

import SpeedyWeather.RingGrids
import SpeedyWeather: round!
import SpeedyWeather.Dates: Dates, DateTime, Period, Second, Millisecond

# default Zarr compressor: BloscCompressor with the same default level we use
# for NetCDF output. Users can override by passing compressor=... to ZarrOutput.
default_zarr_compressor() = Zarr.BloscCompressor(clevel = DEFAULT_COMPRESSION_LEVEL)

resolve_compressor(c) = c
resolve_compressor(::Nothing) = default_zarr_compressor()

"""$(TYPEDSIGNATURES)
Constructor for [`ZarrOutput`](@ref) (extension version, available once `Zarr.jl`
is loaded). Behaves like the [`SpeedyWeather.NetCDFOutput`](@ref) constructor:
builds the output grid, the interpolator, scratch fields and pre-fills the
`output.variables` dictionary with the defaults for `Model`."""
function ZarrOutput(
        SG::SpectralGrid,
        Model::Type{<:AbstractModel} = Barotropic;
        nlayers_soil = DEFAULT_NLAYERS_SOIL,
        output_grid::AbstractFullGrid = on_architecture(
            CPU(), RingGrids.full_grid_type(SG.grid)(SG.grid.nlat_half)
        ),
        output_NF::DataType = DEFAULT_OUTPUT_NF,
        interval::Period = Second(DEFAULT_OUTPUT_INTERVAL),
        compressor = nothing,
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

    # Concrete type parameters: pick the compressor's type (defaulting to
    # the type of `default_zarr_compressor()` for the `nothing` case so that
    # users who later swap the compressor stay within `Union{C, Nothing}`)
    # and the Zarr group type for path-based stores.
    C = typeof(resolve_compressor(compressor))
    Z = Zarr.ZGroup{Zarr.DirectoryStore}

    interval_sec = Second(interval)
    DT = DateTime
    S = typeof(interval_sec)
    F2 = typeof(field2D)
    F3 = typeof(field3D)
    Itp = typeof(interpolator)

    output = ZarrOutput{F2, F3, Itp, DT, S, C, Z}(;
        interval = interval_sec,
        interpolator,
        field2D,
        field3D,
        field3Dland,
        compressor,
        kwargs...
    )

    add_default!(output.variables, Model)
    return output
end

"""$(TYPEDSIGNATURES)
Initialize `ZarrOutput` by creating a Zarr group on disk and storing the initial
conditions of `vars`. Mirrors the layout of [`SpeedyWeather.NetCDFOutput`](@ref):
the store has dimensions `lon`, `lat`, `layer`, `soil_layer`, `time` plus one
chunked array per output variable."""
function initialize!(
        output::ZarrOutput,
        vars::Variables,
        model::AbstractModel,
    )
    output.active || return nothing

    # only checked for models that have a land component
    if hasfield(typeof(model), :land) && !isnothing(model.land)
        @assert SpeedyWeather.get_nlayers(model.land) == size(output.field3Dland, 2) "$(size(output.field3Dland, 2))" *
            " soil layers initialized for output, but $(SpeedyWeather.get_nlayers(model.land)) soil layers initialized for model." *
            " Please construct ZarrOutput with the same `nlayers_soil` as the model."
    end

    # SHARED INITIALIZATION (run folder, output frequency, counters, callbacks)
    initialize!(output.core, output, model)

    # Total number of output snapshots: IC + one per `output_every_n_steps`.
    n_outputs = vars.prognostic.clock.n_timesteps ÷ output.output_every_n_steps + 1

    # CREATE ZARR GROUP (the Zarr store is a *directory*, not a single file)
    (; run_path, filename) = output
    store_path = joinpath(run_path, filename)
    # remove a stale store at the same location when overwrite is on
    output.overwrite && isdir(store_path) && rm(store_path; recursive = true)

    g = Zarr.zgroup(store_path)
    output.zarr_group = g

    # Coordinate arrays. Zarr has no notion of "dimensions" the way NetCDF
    # does — by convention we store the coordinate values as 1D arrays in
    # the group and tag every variable with its `_ARRAY_DIMENSIONS` attribute,
    # which is what Xarray/NetCDF-Zarr uses for dataset round-tripping.
    lond = get_lond(output.field2D)
    latd = get_latd(output.field2D)
    σ = on_architecture(CPU(), model.geometry.σ_levels_full)
    soil_indices = collect(1:get_soil_layers(model))

    write_coordinate!(
        g, "lon", collect(lond);
        attrs = Dict("units" => "degrees_east", "long_name" => "longitude", "_ARRAY_DIMENSIONS" => ["lon"])
    )
    write_coordinate!(
        g, "lat", collect(latd);
        attrs = Dict("units" => "degrees_north", "long_name" => "latitude", "_ARRAY_DIMENSIONS" => ["lat"])
    )
    write_coordinate!(
        g, "layer", collect(σ);
        attrs = Dict("units" => "1", "long_name" => "sigma layer", "_ARRAY_DIMENSIONS" => ["layer"])
    )
    write_coordinate!(
        g, "soil_layer", collect(soil_indices);
        attrs = Dict("units" => "1", "long_name" => "soil layer index", "_ARRAY_DIMENSIONS" => ["soil_layer"])
    )

    # TIME: full-length, chunked by `output.time_chunk`.
    (; startdate) = output
    time_string = "hours since $(Dates.format(startdate, "yyyy-mm-dd HH:MM:0.0"))"
    Zarr.zcreate(
        Float64, g, "time", n_outputs;
        chunks = (max(output.time_chunk, 1),),
        attrs = Dict(
            "units" => time_string, "long_name" => "time",
            "standard_name" => "time", "calendar" => "proleptic_gregorian",
            "_ARRAY_DIMENSIONS" => ["time"]
        ),
    )
    output!(output, vars.prognostic.clock.time)   # write initial time

    # VARIABLES (pre-allocated to full length along the time axis)
    for (key, var) in output.variables
        define_variable!(g, output, var, n_outputs, eltype(output.field2D))
        output!(output, var, Simulation(vars, model))
    end

    return nothing
end

"""$(TYPEDSIGNATURES)
Helper: write a 1D coordinate array `data` to the Zarr group `g` under `name`.
A single chunk is used since the coordinates are small."""
function write_coordinate!(g::Zarr.ZGroup, name::AbstractString, data::AbstractVector; attrs = Dict())
    n = length(data)
    if n == 0
        # Zarr requires chunk size > 0; create a length-0 array with a
        # nominal chunk size of 1.
        z = Zarr.zcreate(eltype(data), g, name, 0; chunks = (1,), attrs = attrs)
    else
        z = Zarr.zcreate(eltype(data), g, name, n; chunks = (n,), attrs = attrs)
        z[:] = data
    end
    return z
end

"""$(TYPEDSIGNATURES)
Define a Zarr array for output `var` in the Zarr group `g`. Shape and chunk
shape are derived from `var.dims_xyzt` and the output grid; the time axis is
pre-allocated to its final length `n_outputs`. Unwritten chunks read back as
the array's `fill_value`."""
function define_variable!(
        g::Zarr.ZGroup,
        output::ZarrOutput,
        var::AbstractOutputVariable,
        n_outputs::Int,
        output_NF::Type{<:AbstractFloat} = DEFAULT_OUTPUT_NF,
    )
    missing_value = hasfield(typeof(var), :missing_value) ? var.missing_value : DEFAULT_MISSING_VALUE

    # Shape per dimension; `false` means the dimension is collapsed away.
    nlon = length(get_lond(output.field2D))
    nlat = length(get_latd(output.field2D))
    nz = is_land(var) ? size(output.field3Dland, 2) : size(output.field3D, 2)
    full_shape = (nlon, nlat, nz, n_outputs)

    # Spatial chunking: 0 (default) or any non-positive value ⇒ full extent.
    # Otherwise clamp to the dimension size so users can't request chunks
    # larger than the array (Zarr requires chunk ≤ shape).
    cx = output.lon_chunk > 0 ? min(output.lon_chunk, nlon) : nlon
    cy = output.lat_chunk > 0 ? min(output.lat_chunk, nlat) : nlat
    cz = output.z_chunk   > 0 ? min(output.z_chunk,   nz)   : nz
    full_chunks = (cx, cy, cz, max(output.time_chunk, 1))
    all_dims = is_land(var) ? ("lon", "lat", "soil_layer", "time") : ("lon", "lat", "layer", "time")

    # Pick out the active dims as flagged by var.dims_xyzt.
    active = var.dims_xyzt
    shape = Tuple(d for (d, on) in zip(full_shape, active) if on)
    chunks = Tuple(c for (c, on) in zip(full_chunks, active) if on)
    # Zarr stores shape/chunks in row-major (C order) but Julia arrays are
    # column-major; Zarr.jl reverses at the metadata boundary. The
    # `_ARRAY_DIMENSIONS` attribute (read by xarray/Xarray-Zarr) must therefore
    # also be in row-major order — i.e. the reverse of our Julia-side dims.
    dims = String[string(d) for (d, on) in zip(all_dims, active) if on]
    reverse!(dims)

    compressor = resolve_compressor(output.compressor)

    attrs = Dict{String, Any}(
        "long_name" => var.long_name,
        "units" => var.unit,
        "_ARRAY_DIMENSIONS" => dims,
    )

    # Zarr stores fill_value in metadata (.zarray); only add a JSON-safe
    # `_FillValue` attribute if the value can be serialized (NaN can't).
    fill = output_NF(missing_value)
    if isfinite(fill)
        attrs["_FillValue"] = fill
    end

    return Zarr.zcreate(
        output_NF, g, var.name, shape...;
        chunks = chunks,
        fill_value = fill,
        compressor = compressor,
        attrs = attrs,
    )
end

"""$(TYPEDSIGNATURES)
Write a single output time step for `variable` to the Zarr store in `output`.
The generic [`output!`](@ref) handles interpolation, transforms and bitrounding;
this method just performs the Zarr-specific store-side write."""
function write_array!(
        output::ZarrOutput,
        variable::AbstractOutputVariable,
        field,
    )
    z = output.zarr_group[variable.name]
    if hastime(variable)
        i = output.output_counter                # current write index
        indices = get_indices(i, variable)
        z[indices...] = parent_array(field)
    else
        # static fields are only written once — just dump the array.
        z[:] = parent_array(field)
    end
    return nothing
end

"""$(TYPEDSIGNATURES)
Write the current time `time::DateTime` to the Zarr store in `output`."""
function output!(output::ZarrOutput, time::DateTime)
    output.output_counter += 1
    i = output.output_counter

    (; startdate) = output
    time_passed = Millisecond(time - startdate)
    time_hrs = time_passed.value / 3600_000
    output.zarr_group["time"][i] = time_hrs
    return nothing
end

"""Pull out the parent (Array) of a Field for direct copy into a Zarr array."""
parent_array(var) = Array(parent(var))

Base.close(output::ZarrOutput) = nothing

end # module
