module SpeedyWeatherTerrariumZarrExt

# Zarr output support for TerrariumOutputVariable (defined in
# SpeedyWeatherTerrariumExt). Only loads once both Terrarium.jl and Zarr.jl
# are loaded (in either order), since it needs types/functions from both.

using SpeedyWeather
using Terrarium
using Zarr
using DocStringExtensions

import SpeedyWeather: ZarrOutput, is3D, hastime, write_array!, define_variable!
import SpeedyWeather.RingGrids

const TerrariumExt = Base.get_extension(SpeedyWeather, :SpeedyWeatherTerrariumExt)
const TerrariumOutputVariable = TerrariumExt.TerrariumOutputVariable
const terrarium_output_field! = TerrariumExt.terrarium_output_field!
const SOIL_DEPTH_DIM = TerrariumExt.SOIL_DEPTH_DIM

# SpeedyWeatherZarrExt's compressor resolution helper, reused so Terrarium
# variables get the same default (and user-overridden) compressor as every
# other Zarr output variable.
const ZarrExt = Base.get_extension(SpeedyWeather, :SpeedyWeatherZarrExt)
const resolve_compressor = ZarrExt.resolve_compressor

"""$(TYPEDSIGNATURES)
Output `variable` into the Zarr store of `output`: gather, interpolate,
transform (via `terrarium_output_field!`, shared with `NetCDFOutput`) and
write the Zarr-specific array slice."""
function SpeedyWeather.output!(
        output::ZarrOutput,
        variable::TerrariumOutputVariable,
        simulation::SpeedyWeather.AbstractSimulation,
    )
    # escape immediately after first call if variable doesn't have a time dimension
    ~hastime(variable) && output.output_counter > 1 && return nothing
    out = terrarium_output_field!(output, variable, simulation)
    isnothing(out) && return nothing
    write_array!(output, variable, out)
    return nothing
end

"""$(TYPEDSIGNATURES)
Define `variable` as a Zarr array in the group `g` of `output`. 3D variables
use their own vertical dimension `soil_depth` (shared between all Terrarium
output variables, and with the netCDF backend) with the depths of the
Terrarium soil layer centres [m, positive down] as coordinates."""
function SpeedyWeather.define_variable!(
        g::Zarr.ZGroup,
        output::ZarrOutput,
        variable::TerrariumOutputVariable,
        n_outputs::Int,
        output_NF::Type{<:AbstractFloat} = SpeedyWeather.DEFAULT_OUTPUT_NF,
    )
    if is3D(variable)   # lazily define the shared soil depth coordinate array
        if !haskey(g, SOIL_DEPTH_DIM)
            Zarr.zcreate(
                Float64, g, SOIL_DEPTH_DIM, variable.nlayers;
                chunks = (variable.nlayers,),
                attrs = Dict(
                    "units" => "m", "long_name" => "depth of soil layer centre",
                    "positive" => "down", "_ARRAY_DIMENSIONS" => [SOIL_DEPTH_DIM],
                ),
            )
            g[SOIL_DEPTH_DIM][:] = variable.depths
        elseif length(g[SOIL_DEPTH_DIM]) != variable.nlayers
            error("Terrarium output variable $(variable.name) has $(variable.nlayers) soil layers, " *
                "but the $SOIL_DEPTH_DIM dimension already has $(length(g[SOIL_DEPTH_DIM])).")
        end
    end

    nlon = length(RingGrids.get_lond(output.field2D))
    nlat = length(RingGrids.get_latd(output.field2D))
    nz = is3D(variable) ? variable.nlayers : 1
    full_shape = (nlon, nlat, nz, n_outputs)

    cx = output.lon_chunk > 0 ? min(output.lon_chunk, nlon) : nlon
    cy = output.lat_chunk > 0 ? min(output.lat_chunk, nlat) : nlat
    cz = output.vertical_chunk > 0 ? min(output.vertical_chunk, nz) : nz
    full_chunks = (cx, cy, cz, max(output.time_chunk, 1))
    all_dims = ("lon", "lat", SOIL_DEPTH_DIM, "time")

    active = variable.dims_xyzt
    shape = Tuple(d for (d, on) in zip(full_shape, active) if on)
    chunks = Tuple(c for (c, on) in zip(full_chunks, active) if on)
    dims = String[string(d) for (d, on) in zip(all_dims, active) if on]

    if output.ensemble_index > 0
        shape = (shape..., output.ensemble_size)
        chunks = (chunks..., 1)
        push!(dims, "ensemble")
    end
    reverse!(dims)   # Zarr stores dims in row-major (C) order, see SpeedyWeatherZarrExt

    compressor = resolve_compressor(output.compressor)
    fill = output_NF(variable.missing_value)
    attrs = Dict{String, Any}(
        "long_name" => variable.long_name, "units" => variable.unit,
        "_ARRAY_DIMENSIONS" => dims,
    )
    isfinite(fill) && (attrs["_FillValue"] = fill)

    return Zarr.zcreate(
        output_NF, g, variable.name, shape...;
        chunks = chunks, fill_value = fill, compressor = compressor, attrs = attrs,
    )
end

end # module
