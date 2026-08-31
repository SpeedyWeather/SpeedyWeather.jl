export HEALPixOutput

"""Output writer that writes a SpeedyWeather simulation to a Zarr store on a
[`HEALPixGrid`](@ref), keeping the horizontal dimension *flat* (unravelled into a
single `cell` dimension) exactly as a `Field`'s data is.
Data on any other grid is interpolated onto HEALPix; a simulation that already runs
on the target HEALPix grid is written without interpolation.

Alongside the variables the store holds one flat vector per cell with its latitude
(`lat`), its longitude (`lon`) and its ring index (`ring`), so consumers need no
knowledge of SpeedyWeather's grid machinery to place the data on the sphere. Cells
are in HEALPix RING order, i.e. the flat index `ij` is the standard RING pixel index
`ij - 1` (0-based) and the store can be read by healpy/cuHPX directly.

The actual implementation lives in the `SpeedyWeatherZarrExt` extension and is only
available once `Zarr.jl` is loaded:

```julia
using Zarr
using SpeedyWeather
spectral_grid = SpectralGrid(truncation = 31, nlayers = 8)
output = HEALPixOutput(spectral_grid, PrimitiveWet, nside = 16)
```

Type parameters: `Field2D`, `Field3D` are the scratch field types on the output
HEALPix grid, `Interpolator` is the interpolator type (or `Nothing` when the model
grid already is the output grid and interpolation is skipped), `DT` and `S` are the
start-date and output-step types, `C` is the Zarr compressor type (or `Nothing` for
the Zarr default), and `Z` is the Zarr group type once `initialize!` has been called
(`Nothing` before). Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct HEALPixOutput{
        Field2D,
        Field3D,
        Interpolator,
        DT,
        S,
        C,
        Z,
    } <: AbstractZarrOutput

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

    "[OPTION] name of the output zarr store (a directory)"
    filename::String = "output_healpix.zarr"

    "[OPTION] also write restart file if output=true?"
    write_restart::Bool = true

    "[OPTION] also write parameters txt file if output=true?"
    write_parameters_txt::Bool = true

    "[OPTION] also write progress txt file if output=true?"
    write_progress_txt::Bool = true

    # WHAT/WHEN OPTIONS
    "[DERIVED] start date of the simulation, used for time dimension in zarr store"
    startdate::DT = DateTime(2000, 1, 1)

    "[OPTION] output frequency, time step"
    interval::S = Second(DEFAULT_OUTPUT_INTERVAL)

    "[OPTION] dictionary of variables to output, e.g. u, v, vor, div, pres, temp, humid"
    variables::OUTPUT_VARIABLES_DICT = OutputVariablesDict()

    "[OPTION] number of time steps per chunk along the time dimension"
    time_chunk::Int = 1

    "[OPTION] chunk size along the flat cell dimension. 0 (default) means one chunk = all cells."
    cell_chunk::Int = 0

    "[OPTION] chunk size along the vertical (layer / soil_layer). 0 (default) means full extent."
    vertical_chunk::Int = 0

    "[OPTION] Zarr compressor (extension-typed). `nothing` keeps the Zarr default."
    compressor::Union{C, Nothing} = nothing

    "[DERIVED] the Zarr group to be written into, created on initialize!"
    zarr_group::Union{Z, Nothing} = nothing

    "[DERIVED] per-variable time-chunk write buffers keyed by variable name, see `ZarrOutput`"
    time_buffers::Dict{String, Any} = Dict{String, Any}()

    "[DERIVED] interpolator onto the output HEALPix grid, `nothing` if the model already runs on it"
    const interpolator::Interpolator
    const land_fraction::Field2D

    # SCRATCH FIELDS ON THE OUTPUT HEALPIX GRID TO INTERPOLATE (OR COPY) ONTO
    const field2D::Field2D
    const field3D::Field3D
    const field3Dland::Field3D
end

"""$(TYPEDSIGNATURES)
Stub constructor for [`HEALPixOutput`](@ref). Errors with a helpful message until the
`Zarr.jl` extension is loaded, at which point the extension installs the real
constructor."""
function HEALPixOutput(SG::SpectralGrid, args...; kwargs...)
    Base.get_extension(@__MODULE__, :SpeedyWeatherZarrExt) === nothing && error(
        "HEALPixOutput requires Zarr.jl to be loaded. Add `using Zarr` (or " *
            "`import Zarr`) before constructing a HEALPixOutput."
    )
    # When the extension is loaded its `HEALPixOutput` method takes precedence and
    # this fallback is unreachable; the throw guards against being called via a
    # generic dispatch path before the extension has registered its method.
    throw(MethodError(HEALPixOutput, (SG, args...)))
end

"""$(TYPEDSIGNATURES)
`true` if `output` writes its data without interpolating, i.e. the simulation already
runs on the output's HEALPix grid."""
skips_interpolation(output::HEALPixOutput) = isnothing(output.interpolator)

function Base.show(io::IO, output::HEALPixOutput{F}) where {F}

    F_str = string("{", F, "}")
    type_param_str = length(F_str) > 30 ? string(first(F_str, 30), "...}") : F_str
    active = output.active ? "active" : "inactive/uninitialized"

    grid = output.field2D.grid
    nlat = get_nlat(grid)
    nside = RingGrids.nside_healpix(grid.nlat_half)
    npix = get_npoints(grid)

    println(io, styled"{warning:HEALPixOutput}{note:$type_param_str}")
    println(io, styled"├ {info:status}: $active")
    println(io, styled"├ {info:write restart file} = $(output.write_restart) (if active)")
    println(io, styled"├ {info:output grid} = $nlat-ring HEALPixGrid {note:(nside=$nside, $npix cells)}")

    if skips_interpolation(output)
        println(io, styled"├ {info:interpolation} = skipped (model already on this grid)")
    else
        interp_type_str = string(typeof(output.interpolator))
        interp_type_str_short = length(interp_type_str) > 70 ? string(first(interp_type_str, 70), "...}") : interp_type_str
        println(io, styled"├ {info:interpolator}::$interp_type_str_short")
    end

    println(io, styled"├ {info:path} = $(joinpath(output.run_path, output.filename)) (overwrite=$(output.overwrite))")
    println(io, styled"├ {info:interval} = $(output.interval)")
    println(io, styled"├ {info:time chunk} = $(output.time_chunk)")
    println(io, styled"├ {info:chunks (cell,vertical)} = ($(output.cell_chunk),$(output.vertical_chunk)) (0 ⇒ full extent)")
    print(io, styled"└ {info:variables}")
    nvars = length(output.variables)
    for (i, (key, var)) in enumerate(output.variables)
        print(io, "\n  $(i == nvars ? "└" : "├") ", styled"{magenta:$key}: $(var.long_name) ", styled"{note:[$(var.unit)]}")
    end
    return nothing
end
