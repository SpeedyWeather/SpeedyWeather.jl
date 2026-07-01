export ZarrOutput

"""Output writer that writes a SpeedyWeather simulation to a Zarr store. Mirrors the
behaviour of [`NetCDFOutput`](@ref) but writes a chunked, optionally compressed Zarr
hierarchy on disk. The actual implementation lives in the `SpeedyWeatherZarrExt`
extension and is only available once `Zarr.jl` is loaded:

```julia
using Zarr
using SpeedyWeather
output = ZarrOutput(spectral_grid)
```

Type parameters: `Field2D`, `Field3D` are the scratch field types, `Interpolator`
is the interpolator type, `DT` and `S` are the start-date and output-step types,
`C` is the Zarr compressor type (or `Nothing` for the Zarr default), and `Z` is
the Zarr group type once `initialize!` has been called (`Nothing` before)."""
@kwdef mutable struct ZarrOutput{
        Field2D,
        Field3D,
        Interpolator,
        DT,
        S,
        C,
        Z,
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

    "[OPTION] name of the output zarr store (a directory)"
    filename::String = "output.zarr"

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

    "[OPTION] chunk size along longitude. 0 (default) means one chunk = full lon extent."
    lon_chunk::Int = 0

    "[OPTION] chunk size along latitude. 0 (default) means one chunk = full lat extent."
    lat_chunk::Int = 0

    "[OPTION] chunk size along the vertical (layer / soil_layer). 0 (default) means full extent."
    vertical_chunk::Int = 0

    # ENSEMBLE OPTIONS
    "[OPTION] this writer's ensemble member index. 0 (default) disables ensemble output; >0 adds an ensemble dimension and this member writes into slot `ensemble_index`. Members are indexed 1..ensemble_size; member 1 creates the shared store."
    ensemble_index::Int = 0

    "[OPTION] total number of ensemble members, sizes the ensemble dimension. Must be ≥ ensemble_index when ensemble output is on."
    ensemble_size::Int = 0

    "[OPTION] seconds a non-creator ensemble member waits for the creator (member 1) to finish building the shared store before erroring."
    ensemble_timeout::Int = 600

    "[OPTION] Zarr compressor (extension-typed). `nothing` keeps the Zarr default."
    compressor::Union{C, Nothing} = nothing

    "[DERIVED] the Zarr group to be written into, created on initialize!"
    zarr_group::Union{Z, Nothing} = nothing

    const interpolator::Interpolator

    # SCRATCH FIELDS TO INTERPOLATE ONTO
    const field2D::Field2D
    const field3D::Field3D
    const field3Dland::Field3D
end

"""$(TYPEDSIGNATURES)
Stub constructor for [`ZarrOutput`](@ref). Errors with a helpful message until the
`Zarr.jl` extension is loaded, at which point the extension installs the real
constructor."""
function ZarrOutput(SG::SpectralGrid, args...; kwargs...)
    Base.get_extension(@__MODULE__, :SpeedyWeatherZarrExt) === nothing && error(
        "ZarrOutput requires Zarr.jl to be loaded. Add `using Zarr` (or " *
            "`import Zarr`) before constructing a ZarrOutput."
    )
    # When the extension is loaded its `ZarrOutput` method takes precedence and
    # this fallback is unreachable; the throw guards against being called via a
    # generic dispatch path before the extension has registered its method.
    throw(MethodError(ZarrOutput, (SG, args...)))
end

function Base.show(io::IO, output::ZarrOutput{F}) where {F}

    F_str = string("{", F, "}")
    type_param_str = length(F_str) > 30 ? string(first(F_str, 30), "...}") : F_str
    active = output.active ? "active" : "inactive/uninitialized"

    println(io, styled"{warning:ZarrOutput}{note:$type_param_str}")
    println(io, styled"├ {info:status}: $active")
    println(io, styled"├ {info:write restart file} = $(output.write_restart) (if active)")

    interp_type_str = string(typeof(output.interpolator))
    interp_type_str_short = length(interp_type_str) > 70 ? string(first(interp_type_str, 70), "...}") : interp_type_str

    println(io, styled"├ {info:interpolator}::$interp_type_str_short")
    println(io, styled"├ {info:path} = $(joinpath(output.run_path, output.filename)) (overwrite=$(output.overwrite))")
    println(io, styled"├ {info:interval} = $(output.interval)")
    println(io, styled"├ {info:time chunk} = $(output.time_chunk)")
    println(io, styled"├ {info:spatial chunks (lon,lat,vertical)} = ($(output.lon_chunk),$(output.lat_chunk),$(output.vertical_chunk)) (0 ⇒ full extent)")
    ensemble_str = output.ensemble_index > 0 ? "member $(output.ensemble_index)/$(output.ensemble_size)" : "off"
    println(io, styled"├ {info:ensemble} = $ensemble_str")
    print(io, styled"└ {info:variables}")
    nvars = length(output.variables)
    for (i, (key, var)) in enumerate(output.variables)
        print(io, "\n  $(i == nvars ? "└" : "├") ", styled"{magenta:$key}: $(var.long_name) ", styled"{note:[$(var.unit)]}")
    end
    return nothing
end
