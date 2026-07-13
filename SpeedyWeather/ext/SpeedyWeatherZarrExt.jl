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
    AbstractFullGrid, run_folder_name

import SpeedyWeather.RingGrids
import SpeedyWeather: round!
import SpeedyWeather.Printf
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
chunked array per output variable.

When ensemble output is on (`output.ensemble_index > 0`) an additional `ensemble`
dimension of length `ensemble_size` is added, chunked with size 1 so that each member
writes disjoint chunk files. All members share a deterministic run folder and hence
one store; member 1 (the *creator*) builds the store schema, coordinates and the shared
`time` axis, then drops a readiness marker; members `2..ensemble_size` (the *writers*) wait
for that marker, open the existing store and write only their own ensemble slice."""
function initialize!(
        output::ZarrOutput,
        vars::Variables,
        model::AbstractModel,
    )
    output.active || return nothing

    # only checked for models that have a land component and output variables that
    # actually use the soil vertical dimension (and its scratch field `field3Dland`)
    if hasfield(typeof(model), :land) && !isnothing(model.land) &&
            any(var -> is_land(var) && is3D(var), values(output.variables))
        @assert SpeedyWeather.get_nlayers(model.land) == size(output.field3Dland, 2) "$(size(output.field3Dland, 2))" *
            " soil layers initialized for output, but $(SpeedyWeather.get_nlayers(model.land)) soil layers initialized for model." *
            " Please construct ZarrOutput with the same `nlayers_soil` as the model."
    end

    ensemble = output.ensemble_index > 0
    if ensemble
        @assert 1 <= output.ensemble_index <= output.ensemble_size "ensemble_index=$(output.ensemble_index)" *
            " must be in 1..ensemble_size=$(output.ensemble_size). Set ensemble_size to the total number of members."
    end

    # SHARED INITIALIZATION (output frequency, counters, callbacks). For ensemble output
    # we manage the (shared, deterministic) run folder ourselves so that all members
    # resolve to the same store path instead of auto-incrementing into separate folders.
    initialize!(output.core, output, model; create_folder = !ensemble)
    ensemble && setup_ensemble_run_folder!(output)

    # Total number of output snapshots: IC + one per `output_every_n_steps`.
    n_outputs = vars.prognostic.clock.n_time_steps ÷ output.output_every_n_steps + 1

    # The Zarr store is a *directory*, not a single file.
    (; run_path, filename) = output
    store_path = joinpath(run_path, filename)

    # WRITER member (ensemble_index > 1): wait for the creator to build the store, then
    # open it and write only this member's ensemble slice. The schema, coordinates and
    # shared time axis are owned by the creator.
    if ensemble && output.ensemble_index != 1
        wait_for_ensemble_store(output, store_path)
        output.zarr_group = Zarr.zopen(store_path, "w")
        for (key, var) in output.variables
            output!(output, var, Simulation(vars, model))
        end
        return nothing
    end

    # CREATOR member (ensemble_index == 1) or non-ensemble output: build the full store.
    # Remove a stale readiness marker + store first when overwrite is on.
    ensemble && rm(ensemble_marker_path(output); force = true)
    output.overwrite && isdir(store_path) && rm(store_path; recursive = true)

    g = Zarr.zgroup(store_path)
    output.zarr_group = g

    # Coordinate arrays. Zarr has no notion of "dimensions" the way NetCDF does — by
    # convention we store the coordinate values as 1D arrays in the group and tag every
    # variable with its `_ARRAY_DIMENSIONS` attribute, which is what Xarray/NetCDF-Zarr
    # uses for dataset round-tripping.
    write_zarr_coordinates!(g, output, model)
    if ensemble
        write_coordinate!(
            g, "ensemble", collect(1:output.ensemble_size);
            attrs = Dict("units" => "1", "long_name" => "ensemble member", "_ARRAY_DIMENSIONS" => ["ensemble"])
        )
    end

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

    # readiness marker so writer members (ensemble_index > 1) may proceed
    ensemble && touch(ensemble_marker_path(output))

    return nothing
end

"""$(TYPEDSIGNATURES)
Write the spatial coordinate arrays (`lon`, `lat`, `layer`, `soil_layer`) into the Zarr
group `g` for `output`, shared by the ensemble and non-ensemble store layouts."""
function write_zarr_coordinates!(g::Zarr.ZGroup, output::ZarrOutput, model::AbstractModel)
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
    return nothing
end

"""$(TYPEDSIGNATURES)
Path to the readiness marker file dropped by the ensemble creator (member 1) once the
shared store schema, coordinates and time axis are in place."""
ensemble_marker_path(output::ZarrOutput) = joinpath(output.run_path, ".zarr_ensemble_ready")

"""$(TYPEDSIGNATURES)
Set a deterministic (non-auto-incrementing) run folder for an ensemble `output` so all
members resolve to the same store path, and create it (idempotent, tolerant of concurrent
creation by other members)."""
function setup_ensemble_run_folder!(output::ZarrOutput)
    fmt = Printf.Format("%0$(output.run_digits)d")
    output.run_folder = run_folder_name(output.run_prefix, output.id, output.run_number; fmt)
    run_path = joinpath(output.path, output.run_folder)
    mkpath(run_path)    # idempotent; safe if another member created it first
    output.run_path = run_path
    return run_path
end

"""$(TYPEDSIGNATURES)
Block until the ensemble creator (member 1) has written the readiness marker for the
shared store at `store_path`, polling at 1s intervals and erroring after
`output.ensemble_timeout` seconds."""
function wait_for_ensemble_store(output::ZarrOutput, store_path::AbstractString)
    marker = ensemble_marker_path(output)
    t0 = time()
    while !isfile(marker)
        (time() - t0) > output.ensemble_timeout && error(
            "ZarrOutput ensemble member $(output.ensemble_index) timed out after " *
                "$(output.ensemble_timeout)s waiting for member 1 to create the shared store at " *
                "$store_path. Ensure the member with ensemble_index=1 is running."
        )
        sleep(1)
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
    cz = output.vertical_chunk > 0 ? min(output.vertical_chunk, nz) : nz
    full_chunks = (cx, cy, cz, max(output.time_chunk, 1))
    all_dims = is_land(var) ? ("lon", "lat", "soil_layer", "time") : ("lon", "lat", "layer", "time")

    # Pick out the active dims as flagged by var.dims_xyzt.
    active = var.dims_xyzt
    shape = Tuple(d for (d, on) in zip(full_shape, active) if on)
    chunks = Tuple(c for (c, on) in zip(full_chunks, active) if on)
    dims = String[string(d) for (d, on) in zip(all_dims, active) if on]

    # Ensemble output: append the ensemble axis as the outermost (last, Julia
    # column-major) dimension, chunked with size 1 so members write disjoint chunk files.
    if output.ensemble_index > 0
        shape = (shape..., output.ensemble_size)
        chunks = (chunks..., 1)
        push!(dims, "ensemble")
    end

    # Zarr stores shape/chunks in row-major (C order) but Julia arrays are
    # column-major; Zarr.jl reverses at the metadata boundary. The
    # `_ARRAY_DIMENSIONS` attribute (read by xarray/Xarray-Zarr) must therefore
    # also be in row-major order — i.e. the reverse of our Julia-side dims. This puts
    # `ensemble` first, matching the CF realization/ensemble convention.
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
    # time-varying variables index into the current time slot; static fields are written
    # once (index 1 is ignored by get_indices). The ensemble slot, if any, is appended.
    i = hastime(variable) ? output.output_counter : 1
    indices = get_indices(i, variable, output.ensemble_index)
    z[indices...] = parent_array(field)
    return nothing
end

"""$(TYPEDSIGNATURES)
Write the current time `time::DateTime` to the Zarr store in `output`."""
function output!(output::ZarrOutput, time::DateTime)
    output.ensemble_index > 1 && return nothing
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
