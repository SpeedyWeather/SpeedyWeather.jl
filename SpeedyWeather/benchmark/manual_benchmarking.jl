#=
Run the SpeedyWeather benchmark suite on a chosen architecture and write the
result into a shared `README.md` without overwriting other architectures.

Usage:

    julia --project=SpeedyWeather/benchmark manual_benchmarking.jl                # CPU, auto-labelled
    julia --project=SpeedyWeather/benchmark manual_benchmarking.jl gpu            # CUDA GPU
    julia --project=SpeedyWeather/benchmark manual_benchmarking.jl amdgpu         # AMDGPU (HIP graphs forced on)
    julia --project=SpeedyWeather/benchmark manual_benchmarking.jl reactant-cpu   # Reactant on CPU
    julia --project=SpeedyWeather/benchmark manual_benchmarking.jl reactant-gpu   # Reactant on CUDA GPU

An optional second argument multiplies the number of timesteps per timed run
(default 1), for longer, more robust publication-ready timings, e.g.

    julia --project=SpeedyWeather/benchmark manual_benchmarking.jl gpu 10         # 10× longer timed runs

Debug mode (`--debug`, anywhere in the arguments) is a quick check whether the
performance has regressed: it only runs the PrimitiveWet resolution sweep
(`:benchmark201`) with 8 layers up to truncation 128, prints the results and leaves the
`README.md` and `assets/benchmark_results.json` untouched. With `--output=FILE`
the results are additionally written to FILE as JSON, e.g. to compare revisions

    julia --project=SpeedyWeather/benchmark manual_benchmarking.jl --debug
    julia --project=SpeedyWeather/benchmark manual_benchmarking.jl gpu --debug --output=main.json

The CPU label is auto-derived from `Sys.ARCH` (`cpu-arm` for aarch64/arm64,
`cpu-x86` otherwise). GPU runs require a CUDA-capable device — `using CUDA`
is loaded automatically. `amdgpu` runs require a ROCm-capable device —
`using AMDGPU` is loaded automatically — and explicitly pass
`gpu_graphs = true` to every `SpectralTransform` built for the `BenchmarkSuite`
resolution/model suites, since `default_gpu_graphs` for AMDGPU's
`ROCBackend` is `false` (HIP graphs are opt-in, verified only on datacenter/CDNA
hardware). The Reactant variants additionally `using Reactant` and 
`Reactant.set_default_backend("cpu" | "gpu")`; for those archs the suites fall 
back to `MatrixSpectralTransform` since the FFT-based default is not yet covered
by Reactant's MLIR pipeline.

Per-architecture results are persisted to `assets/benchmark_results.json`.
The `README.md` (kept at the benchmark root) is then regenerated from all
stored architectures: each gets its own section, plus a cross-architecture
overview table of the PrimitiveWet resolution sweep. The docs page
`docs/src/benchmarks.md` is generated at doc-build time from the same JSON
and additionally renders comparison figures — see `docs/make.jl`.
=#

# Flags (`--debug`, `--output=FILE`) may appear anywhere, the remaining positional
# arguments are the architecture and the timestep multiplier
const FLAGS = filter(startswith("--"), ARGS)
const POSITIONAL_ARGS = filter(!startswith("--"), ARGS)
const DEBUG_MODE = "--debug" in FLAGS
const DEBUG_MAX_TRUNCATION = 128
const DEBUG_NLAYERS = 8

# resolve before the `cd` below so that a relative path is relative to the caller's directory
const OUTPUT_PATH = let i = findfirst(startswith("--output="), FLAGS)
    isnothing(i) ? nothing : abspath(chopprefix(FLAGS[i], "--output="))
end

let unknown_flags = filter(flag -> flag != "--debug" && !startswith(flag, "--output="), FLAGS)
    isempty(unknown_flags) || error("Unknown flag(s) $(join(unknown_flags, ", ")). Use --debug and --output=FILE.")
end
isnothing(OUTPUT_PATH) || DEBUG_MODE || error("--output=FILE is only supported together with --debug.")

import Pkg
cd(@__DIR__)
Pkg.activate(".")

const ARCH_ARG = length(POSITIONAL_ARGS) >= 1 ? lowercase(POSITIONAL_ARGS[1]) : ""

# Conditionally load backend packages BEFORE `using SpeedyWeather` so that
# SpeedyWeather's extensions register the right array types and overloads.
# CUDA must be loaded for any Reactant run (even the CPU backend) — Reactant's
# KernelAbstractions kernel-raising path always needs it.
if ARCH_ARG == "gpu" || ARCH_ARG == "reactant-cpu" || ARCH_ARG == "reactant-gpu"
    using CUDA
end
if ARCH_ARG == "amdgpu"
    using AMDGPU
end
if ARCH_ARG == "reactant-cpu" || ARCH_ARG == "reactant-gpu"
    using Reactant
    Reactant.set_default_backend(ARCH_ARG == "reactant-gpu" ? "gpu" : "cpu")
end

using SpeedyWeather, Dates, Printf, InteractiveUtils, JSON3
import SpeedyWeather.SpeedyTransforms: prettymemory

function pick_architecture(arg::AbstractString)
    if arg == "gpu"
        return (SpeedyWeather.GPU(), "gpu-nvidia")
    elseif arg == "amdgpu"
        return (SpeedyWeather.GPU(), "gpu-amd")
    elseif arg == "reactant-cpu"
        return error("Reactant benchmarks are not yet supported. Sorry")
        #return (SpeedyWeather.ReactantDevice(), "reactant-cpu")
    elseif arg == "reactant-gpu"
        return error("Reactant benchmarks are not yet supported. Sorry")
        #return (SpeedyWeather.ReactantDevice(), "reactant-gpu")
    elseif arg == "cpu" || isempty(arg)
        arch_str = String(Sys.ARCH)
        label = (startswith(arch_str, "aarch") || arch_str == "arm64") ? "cpu-arm" : "cpu-x86"
        return (SpeedyWeather.CPU(), label)
    else
        error("Unknown architecture argument `$arg`. Use one of: cpu, gpu, amdgpu, reactant-cpu, reactant-gpu.")
    end
end

const ARCH, ARCH_LABEL = pick_architecture(ARCH_ARG)
@info "Running benchmarks for architecture: $ARCH_LABEL ($(typeof(ARCH)))"

# Optional second argument: a multiplier on the number of timesteps per timed run.
# Lengthens every run for more robust (e.g. publication-ready) timings; defaults
# to 1. Example: `julia manual_benchmarking.jl gpu 10`.
const TIMESTEP_MULTIPLIER = length(POSITIONAL_ARGS) >= 2 ? parse(Float64, POSITIONAL_ARGS[2]) : 1.0
TIMESTEP_MULTIPLIER > 0 || error("timestep multiplier must be > 0, got $TIMESTEP_MULTIPLIER")
@info "Timestep multiplier: $TIMESTEP_MULTIPLIER"

# Free (device) memory between benchmark runs, called outside the timed region
free_memory!() = SpeedyWeather.SpeedyTransforms.clear_fourier_graph_cache!(ARCH)

include("benchmark_suite.jl")
include("define_benchmarks.jl")

if DEBUG_MODE
    resolution_suite = benchmarks[:benchmark201]
    debug_runs = findall(
        i -> resolution_suite.truncation[i] <= DEBUG_MAX_TRUNCATION && resolution_suite.nlayers[i] == DEBUG_NLAYERS,
        1:resolution_suite.nruns,
    )
    benchmarks = Dict{Symbol, AbstractBenchmarkSuite}(:benchmark201 => select_runs(resolution_suite, debug_runs))
    @info "Debug mode: running only :benchmark201 with L = $DEBUG_NLAYERS, T ≤ $DEBUG_MAX_TRUNCATION ($(length(debug_runs)) runs)"
end

for suite in values(benchmarks)
    suite.architecture = ARCH
    # BenchmarkSuiteDynamics has no timed-run step count, so skip it there.
    hasproperty(suite, :timestep_multiplier) && (suite.timestep_multiplier = TIMESTEP_MULTIPLIER)
    # AMDGPU's `default_gpu_graphs` is `false` (HIP graphs are opt-in,
    # so force it on explicitly for this arch;
    # leave every other arch on its own `SpectralTransform` default (`nothing`).
    hasproperty(suite, :gpu_graphs) && (suite.gpu_graphs = ARCH_LABEL == "gpu-amd" ? true : nothing)
end

# Deterministic suite ordering (sort by key for reproducible README sections).
const SUITE_KEYS = sort!(collect(keys(benchmarks)))

for key in SUITE_KEYS
    @info "→ Running suite $key: $(benchmarks[key].title)"
    run_benchmark_suite!(benchmarks[key])
end


function arch_markdown()
    io = IOBuffer()
    for key in SUITE_KEYS
        write_results(io, benchmarks[key])
    end
    return String(take!(io))
end

function machine_info()
    io = IOBuffer()
    write(io, "```julia\njulia> versioninfo()\n")
    versioninfo(io)
    write(io, "```\n")
    if ARCH_LABEL == "gpu-amd"
        write(io, "\n```julia\njulia> AMDGPU.versioninfo()\n")
        try
            AMDGPU.versioninfo(io)
        catch err
            write(io, "(AMDGPU.versioninfo() failed: $err)\n")
        end
        write(io, "```\n")
    elseif ARCH isa SpeedyWeather.GPU || ARCH_LABEL == "reactant-gpu"
        write(io, "\n```julia\njulia> CUDA.versioninfo()\n")
        try
            CUDA.versioninfo(io)
        catch err
            write(io, "(CUDA.versioninfo() failed: $err)\n")
        end
        write(io, "```\n")
    end
    if ARCH isa SpeedyWeather.ReactantDevice
        backend = ARCH_LABEL == "reactant-gpu" ? "gpu" : "cpu"
        write(io, "\nReactant backend: `$backend` (via `Reactant.set_default_backend(\"$backend\")`).\n")
    end
    return String(take!(io))
end

# Pull the structured numbers needed for the cross-architecture comparison
# from the PrimitiveWet resolution suite (:benchmark201).
function overview_data()
    suite = benchmarks[:benchmark201]
    # JSON has no Inf/NaN literal — emit `null` for any non-finite metric (e.g. the SYPD of an
    # unstable config) so JSON3 can write the file instead of erroring. The README generator
    # already renders missing/non-finite entries as "—".
    json_safe(x) = (x isa Number && isfinite(x)) ? x : nothing
    return Dict(
        "truncation" => collect(suite.truncation),
        "nlayers" => collect(suite.nlayers),
        "nlat" => collect(suite.nlat),
        "sypd" => map(json_safe, suite.SYPD),
        "memory" => map(json_safe, suite.memory),
        "dt" => map(json_safe, suite.Δt),
        "spectral_transform" => string.(suite.spectral_transform),
    )
end

arch_record = Dict(
    "meta" => Dict(
        "speedyweather_version" => string(SpeedyWeather.pkgversion(SpeedyWeather)),
        "timestamp" => Dates.format(Dates.now(), Dates.RFC1123Format),
        "arch_type" => string(typeof(ARCH)),
        "machine_info" => machine_info(),
    ),
    "markdown" => arch_markdown(),
    "overview" => overview_data(),
)

# Debug mode: print the results and write them to `--output` (if given) to compare
# revisions, but leave benchmark_results.json and README.md untouched
if DEBUG_MODE
    print(arch_record["markdown"])
    if !isnothing(OUTPUT_PATH)
        meta = arch_record["meta"] = Dict{String, Any}(arch_record["meta"])
        meta["arch_label"] = ARCH_LABEL
        meta["debug_max_truncation"] = DEBUG_MAX_TRUNCATION
        meta["debug_nlayers"] = DEBUG_NLAYERS
        meta["timestep_multiplier"] = TIMESTEP_MULTIPLIER
        # where the monorepo packages were loaded from, to verify the benchmarked revision
        meta["package_dirs"] = Dict(
            string(name) => pkgdir(getfield(SpeedyWeather, name))
                for name in (:SpeedyWeather, :SpeedyWeatherInternals, :LowerTriangularArrays, :RingGrids, :SpeedyTransforms)
                if isdefined(SpeedyWeather, name)
        )
        mkpath(dirname(OUTPUT_PATH))
        open(io -> JSON3.pretty(io, arch_record), OUTPUT_PATH, "w")
        @info "Wrote debug results → $OUTPUT_PATH"
    end
    exit()
end

# Merge into the JSON store

const ASSETS_DIR = joinpath(@__DIR__, "assets")
mkpath(ASSETS_DIR)
const RESULTS_JSON = joinpath(ASSETS_DIR, "benchmark_results.json")

function load_results()
    isfile(RESULTS_JSON) || return Dict{String, Any}()
    try
        return Dict{String, Any}(JSON3.read(read(RESULTS_JSON, String), Dict{String, Any}))
    catch err
        @warn "Could not parse $RESULTS_JSON; starting fresh" exception = err
        return Dict{String, Any}()
    end
end

all_results = load_results()
all_results[ARCH_LABEL] = arch_record

open(RESULTS_JSON, "w") do io
    JSON3.pretty(io, all_results)
end
@info "Wrote results for $ARCH_LABEL → $RESULTS_JSON"

# Regenerate README.md from the JSON store

# Stable ordering for arch sections in the README.
const ARCH_ORDER = ["cpu-arm", "cpu-x86", "gpu-nvidia", "gpu-amd", "reactant-cpu", "reactant-gpu"]

function sorted_arch_labels(results)
    known = filter(in(keys(results)), ARCH_ORDER)
    extra = sort!([k for k in keys(results) if !(k in ARCH_ORDER)])
    return vcat(known, extra)
end

function write_preamble(md)
    write(md, "# Benchmarks\n\n")
    write(md, "Performance benchmarks for SpeedyWeather.jl, collected across multiple architectures. ")
    write(md, "Each architecture's results live in its own section below; the overview table at the top compares the headline ")
    write(md, "PrimitiveWet resolution sweep across all archs that have been benchmarked so far.\n\n")

    write(md, "All simulations are benchmarked over several seconds (wallclock time) without output. ")
    write(md, "Benchmarking excludes initialization and is started just before the main time loop and finishes right after. ")
    write(md, "The benchmarking results here are not very robust; timings that change by ±50% are not uncommon. ")
    write(md, "Proper benchmarking for performance optimization uses the minimum or median of many executions, while we run ")
    write(md, "a simulation for several time steps which effectively represents the mean, susceptible to outliers. ")
    write(md, "However, this is what a user will experience in most situations anyway and the following therefore presents a ")
    write(md, "rough idea of how fast a SpeedyWeather simulation will run, and how much memory it requires.\n\n")

    write(md, "### Explanation\n\n")
    write(md, "Abbreviations in the tables below are as follows; omitted columns use defaults.\n")
    write(md, "- NF: Number format, default: $(SpeedyWeather.DEFAULT_NF)\n")
    write(md, "- T: Spectral resolution, maximum degree of spherical harmonics (1-based), default: T$(SpeedyWeather.DEFAULT_TRUNCATION)\n")
    write(md, "- L: Number of vertical layers, default: $(SpeedyWeather.DEFAULT_NLAYERS) (for 3D models)\n")
    write(md, "- Grid: Horizontal grid, default: $(SpeedyWeather.DEFAULT_GRID)\n")
    write(md, "- Rings: Grid-point resolution, number of latitude rings pole to pole\n")
    write(md, "- Dynamics: With dynamics?, default: true\n")
    write(md, "- Physics: With physical parameterizations?, default: true (for primitive equation models)\n")
    write(md, "- Δt: time step [s].\n")
    write(md, "- SYPD: Speed of simulation, simulated years per wallclock day.\n")
    write(md, "- Memory: Memory footprint of simulation, variables and constants.\n\n")

    write(md, "### Running the benchmarks\n\n")
    write(md, "Reproduce the benchmark suite by running, from `SpeedyWeather/benchmark`:\n\n")
    write(md, "```\n")
    write(md, "julia --project=. manual_benchmarking.jl                # CPU (auto-labelled cpu-arm or cpu-x86)\n")
    write(md, "julia --project=. manual_benchmarking.jl gpu            # CUDA GPU\n")
    write(md, "julia --project=. manual_benchmarking.jl amdgpu         # AMDGPU (HIP graphs forced on)\n")
    write(md, "julia --project=. manual_benchmarking.jl reactant-cpu   # Reactant on CPU\n")
    write(md, "julia --project=. manual_benchmarking.jl reactant-gpu   # Reactant on CUDA GPU\n")
    write(md, "julia --project=. manual_benchmarking.jl --debug        # quick regression check: PrimitiveWet, L$DEBUG_NLAYERS, T ≤ $DEBUG_MAX_TRUNCATION, not stored\n")
    write(md, "julia regression/regression.jl check                   # debug mode on main, latest release and latest benchmarked revision\n")
    write(md, "```\n\n")
    write(md, "Each run updates only its own architecture's section in this `README.md`; results for other architectures are preserved via `benchmark_results.json`.\n\n")
    return
end

# Map the stored transform symbol ("default"/"matrix") to the short label used
# in user-facing tables and figures.
transform_short_label(s) = String(s) == "matrix" ? "MT" : "LT+FFT"

function write_overview(md, all_results, labels)
    write(md, "## Overview: PrimitiveWet resolution across architectures\n\n")
    write(md, "Simulated years per wallclock day (SYPD) for the `PrimitiveWetModel` resolution sweep, ")
    write(md, "one column per architecture. Each (T, L) configuration is reported for both the standard ")
    write(md, "Legendre transform and fast Fourier transform (LT+FFT) and the single matrix transform (MT). ")
    write(md, "Empty cells mean the architecture has not yet been benchmarked or that suite was skipped. ")
    write(md, "Comparison figures across architectures are available on the documentation's `Benchmarks` page.\n\n")

    # Pick the union of (truncation, nlayers, transform) rows from all archs.
    rows = Tuple{Int, Int, String}[]
    for label in labels
        ov = get(all_results[label], "overview", nothing)
        ov === nothing && continue
        truncations = ov["truncation"]
        nlayers = ov["nlayers"]
        transforms = get(ov, "spectral_transform", fill("default", length(truncations)))
        for i in eachindex(truncations)
            r = (Int(truncations[i]), Int(nlayers[i]), String(transforms[i]))
            r in rows || push!(rows, r)
        end
    end
    # LT+FFT (default) before MT (matrix) within the same (T, L) group.
    sort!(rows; by = r -> (r[1], r[2], r[3] == "default" ? 0 : 1))

    header = "| T | L | Transform | " * join(labels, " | ") * " |"
    sep = "| --- | --- | --- | " * join(fill("---", length(labels)), " | ") * " |"
    write(md, header * "\n")
    write(md, sep * "\n")

    for (t, l, tr) in rows
        cells = String[]
        for label in labels
            ov = get(all_results[label], "overview", nothing)
            cell = "—"
            if ov !== nothing
                truncations = ov["truncation"]
                nlayers = ov["nlayers"]
                sypd = ov["sypd"]
                transforms = get(ov, "spectral_transform", fill("default", length(truncations)))
                for i in eachindex(truncations)
                    if Int(truncations[i]) == t && Int(nlayers[i]) == l && String(transforms[i]) == tr
                        s = sypd[i]
                        cell = (s isa Number && isfinite(s)) ? format_sypd(s) : "—"
                        break
                    end
                end
            end
            push!(cells, cell)
        end
        write(md, "| $t | $l | $(transform_short_label(tr)) | " * join(cells, " | ") * " |\n")
    end
    write(md, "\n")
    return
end

function write_arch_section(md, label, record)
    meta = record["meta"]
    write(md, "## Architecture: `$label`\n\n")
    write(md, "Created for SpeedyWeather.jl v$(meta["speedyweather_version"]) on $(meta["timestamp"]).\n\n")
    write(md, "### Machine details\n\n")
    write(md, meta["machine_info"])
    write(md, "\n")
    write(md, record["markdown"])
    write(md, "\n")
    return
end

const README_PATH = joinpath(@__DIR__, "README.md")
labels = sorted_arch_labels(all_results)

open(README_PATH, "w") do md
    write_preamble(md)
    write_overview(md, all_results, labels)
    for label in labels
        write_arch_section(md, label, all_results[label])
    end
end
@info "Regenerated $README_PATH for architectures: $(join(labels, ", "))"
