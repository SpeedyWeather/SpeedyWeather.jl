# Utilities shared by every `AbstractZarrOutput` writer: compressor defaults, coordinate
# and time-axis creation, the time-chunk write buffering and the store-independent parts
# of `initialize!`. The two layouts (`ZarrOutput`'s rectangular lon/lat store and
# `HEALPixOutput`'s flat cell store) differ only in `define_variable!`, in their coordinate
# arrays and in the write indices returned by `zarr_write_indices`.

# default Zarr compressor: BloscCompressor with the same default level we use
# for NetCDF output. Users can override by passing compressor=... to the writer.
default_zarr_compressor() = Zarr.BloscCompressor(clevel = DEFAULT_COMPRESSION_LEVEL)

resolve_compressor(c) = c
resolve_compressor(::Nothing) = default_zarr_compressor()

"""$(TYPEDSIGNATURES)
Ensemble member index of a Zarr-backed output writer, or 0 for writers that don't support
ensemble output (currently everything except [`ZarrOutput`](@ref)). Used by the shared write
path to decide whether an ensemble slot has to be appended to the write indices and whether
this member owns the shared time axis."""
zarr_ensemble_index(::AbstractZarrOutput) = 0

"""$(TYPEDSIGNATURES)
Indices to write the output slice of `variable` at output step `i` into its Zarr array.
Layout hook of the shared write path: the default is the rectangular `lon`/`lat` layout
(with the optional ensemble slot appended), [`HEALPixOutput`](@ref) overrides it with the
flat-cell equivalent."""
zarr_write_indices(output::AbstractZarrOutput, i, variable::AbstractOutputVariable) =
    get_indices(i, variable, zarr_ensemble_index(output))

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
Length of the coordinate array `name` in the Zarr group `g` or `nothing` if not
defined. Zarr-store equivalent of `get_dimension_length(::NCDataset, name)` so that
custom output variables can define their own dimension in `define_dimension!`
with one method for all output backends."""
get_dimension_length(g::Zarr.ZGroup, name::String) = haskey(g, name) ? length(g[name]) : nothing

"""$(TYPEDSIGNATURES)
Define a coordinate in the Zarr group `g`: a 1D array `name` with `values` and
`attribs` as attributes, tagged with its own `_ARRAY_DIMENSIONS`. Zarr-store
equivalent of `define_coordinate!(::NCDataset, ...)` so that custom output
variables can define their own dimension in `define_dimension!` with one
method for all output backends."""
define_coordinate!(g::Zarr.ZGroup, name::String, values::AbstractVector; attribs = Dict{String, String}()) =
    write_coordinate!(g, name, values; attrs = merge(Dict{String, Any}("_ARRAY_DIMENSIONS" => [name]), attribs))

"""$(TYPEDSIGNATURES)
Create the `time` coordinate of the Zarr group `g` for `output`: a full-length array of
`n_outputs` (initial conditions plus one entry per output step), chunked by
`output.time_chunk` and CF-tagged as hours since `output.startdate`."""
function create_time_axis!(g::Zarr.ZGroup, output::AbstractZarrOutput, n_outputs::Integer)
    (; startdate) = output
    time_string = "hours since $(Dates.format(startdate, "yyyy-mm-dd HH:MM:0.0"))"
    return Zarr.zcreate(
        Float64, g, "time", n_outputs;
        chunks = (max(output.time_chunk, 1),),
        attrs = Dict(
            "units" => time_string, "long_name" => "time",
            "standard_name" => "time", "calendar" => "proleptic_gregorian",
            "_ARRAY_DIMENSIONS" => ["time"]
        ),
    )
end

"""$(TYPEDSIGNATURES)
Assert that the number of soil layers allocated for `output` matches the number of soil
layers of `model.land`. Only checked for models that have a land component and for output
variables that actually use the soil vertical dimension (and hence the `field3Dland`
scratch field)."""
function assert_soil_layers(output::AbstractZarrOutput, model::AbstractModel)
    if hasfield(typeof(model), :land) && !isnothing(model.land) &&
            any(var -> is_land(var) && is3D(var), values(output.variables))
        @assert get_nlayers(model.land) == size(output.field3Dland, 2) "$(size(output.field3Dland, 2))" *
            " soil layers initialized for output, but $(get_nlayers(model.land)) soil layers initialized for model." *
            " Please construct $(nameof(typeof(output))) with the same `nlayers_soil` as the model."
    end
    return nothing
end

"""$(TYPEDSIGNATURES)
Create the Zarr array for output `var` in group `g` with the given Julia-order `shape` and
`chunks` and the row-major `dims` for its `_ARRAY_DIMENSIONS` attribute (Zarr stores
shape/chunks in row-major C order but Julia arrays are column-major, and Zarr.jl reverses
at the metadata boundary — so `dims` must be the reverse of the Julia-side dimension
order). Shared by all `AbstractZarrOutput` layouts, which differ only in how they derive
shape, chunks and dims in their `define_variable!` methods."""
function zcreate_output_variable!(
        g::Zarr.ZGroup,
        output::AbstractZarrOutput,
        var::AbstractOutputVariable,
        shape::Tuple,
        chunks::Tuple,
        dims::Vector{String},
        output_NF::Type{<:AbstractFloat},
        missing_value,
    )
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
        compressor = resolve_compressor(output.compressor),
        attrs = attrs,
    )
end

"""$(TYPEDSIGNATURES)
Return `true` if writes for `output` should be buffered along the time dimension.
Buffering only pays off when a Zarr chunk spans more than one time step
(`time_chunk > 1`); with `time_chunk == 1` every time slice is already its own chunk
and is written directly (see [`write_array!`](@ref))."""
time_buffered(output::AbstractZarrOutput) = output.time_chunk > 1

"""$(TYPEDSIGNATURES)
Write a single output time step for `variable` to the Zarr store in `output`.

A Zarr chunk is the atomic unit of (de)compression and I/O, so writing one time
slice into a chunk that spans `time_chunk > 1` steps would force a read-modify-write
of the whole chunk on every step. To avoid this, time-varying variables are buffered
`time_chunk` slices deep (see `output.time_buffers`) and flushed one chunk at a time via
[`flush_time_chunk!`](@ref); a chunk-aligned write is fulfilled by Zarr's
single-chunk write operation with no read. Any remaning chunks in the buffer are flushed
to disk on `close` (see [`flush_partial_time_chunks!`](@ref)).

Buffering is skipped entirely (see [`time_buffered`](@ref)) for static, non-temporal
variables and when `time_chunk == 1`: in these cases, [`write_slice!`](@ref) writes the
data directly to disk."""
function write_array!(
        output::AbstractZarrOutput,
        variable::AbstractOutputVariable,
        field,
    )
    data = parent_array(field)

    # Skip buffering for static variables (written once at index 1) and when
    # time_chunk == 1 (every time slice is its own chunk): write the slice directly.
    if !hastime(variable) || !time_buffered(output)
        i = hastime(variable) ? output.output_counter : 1
        write_slice!(output, variable, data, i)
        return nothing
    end

    # Buffered path: copy this slice into the variable's time buffer (lazily allocated
    # to the chunk shape (spatial..., time_chunk)) at its offset within the current chunk.
    time_chunk = output.time_chunk
    i = output.output_counter
    offset = mod1(i, time_chunk)    # 1..time_chunk position within the current chunk
    buffer = get!(output.time_buffers, variable.name) do
        Array{eltype(data)}(undef, size(data)..., time_chunk)
    end
    selectdim(buffer, ndims(buffer), offset) .= data

    # A full chunk has been collected (this slice closes it, so it is chunk-aligned):
    # flush the whole buffer in one chunk-aligned write.
    offset == time_chunk && flush_time_chunk!(output, variable, buffer, i, time_chunk)
    return nothing
end

"""$(TYPEDSIGNATURES)
Write a single slice `data` for `variable` into its Zarr array at time index `i`
(ignored for static variables). The write indices, including the ensemble slot if any,
come from [`zarr_write_indices`](@ref). This is the unbuffered write path used when
[`time_buffered`](@ref) is `false`."""
function write_slice!(output::AbstractZarrOutput, variable::AbstractOutputVariable, data, i::Integer)
    z = output.zarr_group[variable.name]
    z[zarr_write_indices(output, i, variable)...] = data
    return nothing
end

"""$(TYPEDSIGNATURES)
Write the `count` buffered time slices ending at time index `i_last` for `variable`
from `buffer` to the Zarr store in `output` in a single write covering the time range
`(i_last - count + 1):i_last`. A full-chunk flush (`count == time_chunk`) passes the
buffer `Array` through so Zarr writes it directly."""
function flush_time_chunk!(
        output::AbstractZarrOutput,
        variable::AbstractOutputVariable,
        buffer::AbstractArray,
        i_last::Integer,
        count::Integer,
    )
    t0 = i_last - count + 1
    z = output.zarr_group[variable.name]
    indices = zarr_write_indices(output, t0:i_last, variable)
    if count == size(buffer, ndims(buffer))
        z[indices...] = buffer
    else
        z[indices...] = selectdim(buffer, ndims(buffer), 1:count)
    end
    return nothing
end

"""$(TYPEDSIGNATURES)
Flush the trailing partial time chunk of every time-varying variable to the Zarr store
in `output`, called on `close` after the run. All time-varying variables are written in
lockstep every output step, so they share the same `output_counter`; when it is not a
multiple of `time_chunk` the last `output_counter % time_chunk` slices are still buffered
and are written here. A no-op when `time_chunk == 1` or the final chunk was already full."""
function flush_partial_time_chunks!(output::AbstractZarrOutput)
    (!time_buffered(output) || isempty(output.time_buffers)) && return nothing
    time_chunk = output.time_chunk
    i = output.output_counter
    remainder = i % time_chunk
    remainder == 0 && return nothing    # last chunk already flushed by write_array!
    for variable in values(output.variables)
        haskey(output.time_buffers, variable.name) || continue
        flush_time_chunk!(output, variable, output.time_buffers[variable.name], i, remainder)
    end
    return nothing
end

"""$(TYPEDSIGNATURES)
Write the current time `time::DateTime` to the Zarr store in `output`. Ensemble writer
members (index > 1) share the creator's time axis and skip this."""
function output!(output::AbstractZarrOutput, time::DateTime)
    zarr_ensemble_index(output) > 1 && return nothing
    i = output.output_counter

    (; startdate) = output
    time_passed = Millisecond(time - startdate)
    time_hrs = time_passed.value / 3600_000
    output.zarr_group["time"][i] = time_hrs
    return nothing
end

"""Pull out the parent (Array) of a Field for direct copy into a Zarr array."""
parent_array(var) = Array(parent(var))

function Base.close(output::AbstractZarrOutput)
    flush_partial_time_chunks!(output)
    return nothing
end

# Implemented in the GeoMakie extension, just define here how to load a zarr dataset
SpeedyWeather.animate(::Type{<:Zarr.ZGroup}, path::String; kwargs...) =
    animate(zopen(path); kwargs...)
