#=
Run the SpeedyWeather benchmark suite on a chosen architecture and write the
result into a shared `README.md` without overwriting other architectures.

Usage:

    julia --project=SpeedyWeather/benchmark manual_benchmarking.jl                       # CPU, auto-labelled
    julia --project=SpeedyWeather/benchmark manual_benchmarking.jl gpu                   # CUDA GPU
    julia --project=SpeedyWeather/benchmark manual_benchmarking.jl amdgpu                # AMDGPU (HIP graphs forced on)
    julia --project=SpeedyWeather/benchmark manual_benchmarking.jl reactant-cpu          # Reactant on CPU
    julia --project=SpeedyWeather/benchmark manual_benchmarking.jl reactant-nvidia-gpu   # Reactant on NVIDIA GPU (CUDA)
    julia --project=SpeedyWeather/benchmark manual_benchmarking.jl reactant-amd-gpu      # Reactant on AMD GPU (ROCm)
    julia --project=SpeedyWeather/benchmark manual_benchmarking.jl reactant-tpu          # Reactant on TPU

Every configuration of the `BenchmarkSuite`s is timed up to 5 times from the same
initial state, as long as all its timed runs together take at most 30 s
(`repeats` and `max_repeat_time` in `benchmark_suite.jl`), and the median SYPD is
reported with the spread of the runs. The option `--repeats N` times every
configuration exactly N times instead, without the 30 s budget (not the
reverse-mode suite, which keeps its up to 3 samples).

The option `--length X` (or a second positional argument, as before) multiplies
the number of timesteps per timed run of the main benchmark, the PrimitiveWet
resolution sweep (`:benchmark201`), for longer, more robust publication-ready
timings (default 1, the other suites are not affected), e.g.

    julia --project=SpeedyWeather/benchmark manual_benchmarking.jl gpu --repeats 5 --length 10   # 5 runs each, 10× longer benchmark201

The option `--suites` restricts the run to specific suites, given as comma-separated
keys, e.g. `--suites benchmark201,benchmark900` (default: all suites). A partial run
replaces the architecture's stored results with just the suites it ran.

Options take a value, as `--repeats 5` or `--repeats=5`. The settings are stated at
the top of each architecture's section and stored in its `meta`. The arguments are
parsed with ArgParse.jl, `manual_benchmarking.jl --help` lists them all.

The flag `--reverse` (Reactant only, anywhere among the arguments) additionally
runs the reverse-mode (AD) benchmarks of `reverse_benchmark_suite.jl`, which
take hours, e.g.

    julia --project=SpeedyWeather/benchmark manual_benchmarking.jl reactant-nvidia-gpu --reverse

The flag `--hi-res` adds T342 and T512 at 8, 16 and 24 layers to the main
benchmark (`:benchmark201`), on GPU and TPU architectures only (ignored with a
warning on CPU, too slow).

The flag `--SPEEDY` additionally runs a suite (`:benchmark900`) of the PrimitiveWet
model over the resolutions of the main benchmark (`:benchmark201`, default
transform) in a configuration closer to the original SPEEDY.f90: seasonal
prescribed climatologies of sea surface temperature, land surface temperature and
soil moisture instead of the slab ocean and the land bucket model, constant
instead of learned surface roughness lengths, and a constant wind slowdown (0.95)
and gusts (5 m/s) as in SPEEDY.f90 instead of the learned wind slowdown and
Beljaars gusts (`speedy_components` in `benchmark_suite.jl`). Its section and its
numbers (the `speedy` field of the record) are part of the architecture's results,
like those of `--reverse`. On every architecture, including Reactant.

A run replaces all results stored for its label, including those of suites that a
flag added in an earlier run. Flags and options can be combined in any order, e.g.
`reactant-nvidia-gpu --length 10 --reverse --hi-res`; `--SPEEDY` can also be written `--speedy`. Truncations are 1-based throughout
(T32, …, T256) as in the model code.

The CPU label is auto-derived from `Sys.ARCH` (`cpu-arm` for aarch64/arm64,
`cpu-x86` otherwise). GPU runs require a CUDA-capable device — `using CUDA`
is loaded automatically. `amdgpu` runs require a ROCm-capable device —
`using AMDGPU` is loaded automatically — and explicitly pass
`gpu_graphs = true` to every `SpectralTransform` built for the `BenchmarkSuite`
resolution/model suites, since `default_gpu_graphs` for AMDGPU's
`ROCBackend` is `false` (HIP graphs are opt-in, verified only on datacenter/CDNA
hardware). The Reactant variants additionally `using Reactant` and
`Reactant.set_default_backend("cpu" | "gpu" | "tpu")`, "gpu" for both
`reactant-nvidia-gpu` and `reactant-amd-gpu`: the GPU vendor is that of the
installed Reactant_jll build, which is checked against the label (CUDA or ROCm)
so that results aren't stored under the wrong one. Reactant_jll picks its ROCm
build on a machine with an AMD GPU driver and no CUDA device, or set the
preference `gpu = "rocm"` of Reactant_jll (in `LocalPreferences.toml`) before
instantiating. The time loop is compiled once per configuration (`reactant_runner.jl`),
the compile time is reported separately and not included in the SYPD.
`reactant-tpu` runs on a TPU host (e.g. a Cloud TPU VM, Reactant downloads libtpu
itself); CUDA.jl is loaded there too (and with `reactant-amd-gpu`), without an NVIDIA GPU,
as Reactant's kernel raising needs it. The matrix transform only goes up to T128 on TPU, as T171 needs more 
memory than some TPU chips have.

Per-architecture results are persisted to `assets/benchmark_results.json`.
The `README.md` (kept at the benchmark root) is then regenerated from all
stored architectures: each gets its own section, plus a cross-architecture
overview table of the PrimitiveWet resolution sweep. The docs page
`docs/src/benchmarks.md` is generated at doc-build time from the same JSON
and additionally renders comparison figures — see `docs/make.jl`.
=#

import Pkg
cd(@__DIR__)
Pkg.activate(".")

# Command-line arguments, `manual_benchmarking.jl --help` lists them. Parsed before any backend is
# loaded, so that wrong arguments fail within seconds
using ArgParse

const ARCHITECTURES = ("cpu", "gpu", "amdgpu", "reactant-cpu", "reactant-nvidia-gpu", "reactant-amd-gpu", "reactant-tpu")
const ARGUMENT_SETTINGS = ArgParseSettings(
    description = "Run the SpeedyWeather benchmark suite on one architecture and update README.md " *
        "and assets/benchmark_results.json, see the comment at the top of manual_benchmarking.jl.",
)
@add_arg_table! ARGUMENT_SETTINGS begin
    "architecture"
    help = "cpu (labelled cpu-arm or cpu-x86), gpu, amdgpu, reactant-cpu, reactant-nvidia-gpu, reactant-amd-gpu or reactant-tpu"
    default = "cpu"
    range_tester = architecture -> lowercase(architecture) in ARCHITECTURES
    "multiplier"
    help = "same as --length, as before"
    arg_type = Float64
    required = false
    range_tester = multiplier -> multiplier > 0
    "--length"
    help = "multiply the number of time steps of the timed runs of the PrimitiveWet resolution sweep (benchmark201) by X"
    arg_type = Float64
    metavar = "X"
    range_tester = multiplier -> multiplier > 0
    "--repeats"
    help = "exactly N timed runs per configuration, instead of up to 5 within 30 s (not the reverse-mode suite)"
    arg_type = Int
    metavar = "N"
    range_tester = repeats -> repeats >= 1
    "--suites"
    help = "comma-separated suite keys to run, e.g. `benchmark201,benchmark900` (benchmark800 needs --reverse, " *
        "benchmark900 --SPEEDY). Default: all suites. A partial run replaces this architecture's stored results " *
        "with just these suites."
    arg_type = String
    default = ""
    "--reverse"
    help = "also run the reverse-mode (AD) benchmarks, Reactant only, takes hours"
    action = :store_true
    "--hi-res"
    help = "add T342 and T512 to the PrimitiveWet resolution sweep, GPU and TPU only"
    action = :store_true
    dest_name = "hi_res"
    "--SPEEDY", "--speedy"
    help = "also run the PrimitiveWet resolution sweep in a SPEEDY.f90-like configuration (prescribed seasonal climatologies, constant surface roughness, wind slowdown and gusts)"
    action = :store_true
    dest_name = "speedy"
end
const ARGUMENTS = parse_args(ARGUMENT_SETTINGS)

const RUN_REVERSE = ARGUMENTS["reverse"]
const RUN_HI_RES = ARGUMENTS["hi_res"]
const RUN_SPEEDY = ARGUMENTS["speedy"]

const ARCH_ARG = lowercase(ARGUMENTS["architecture"])
# Reactant's backend per architecture argument, the GPU vendor is that of the installed Reactant_jll build
const REACTANT_BACKENDS = Dict("reactant-cpu" => "cpu", "reactant-nvidia-gpu" => "gpu", "reactant-amd-gpu" => "gpu", "reactant-tpu" => "tpu")
# the platform of Reactant's GPU client that each GPU label requires
const REACTANT_GPU_PLATFORMS = Dict("reactant-nvidia-gpu" => "cuda", "reactant-amd-gpu" => "rocm")
const ON_REACTANT = haskey(REACTANT_BACKENDS, ARCH_ARG)
RUN_REVERSE && !ON_REACTANT &&
    error("The reverse-mode benchmarks (--reverse) need Reactant, use reactant-cpu, reactant-nvidia-gpu, reactant-amd-gpu or reactant-tpu.")

# `--length X` (or a second positional argument): a multiplier on the number of time steps per timed
# run of the main benchmark (:benchmark201) only, for longer, more robust (e.g. publication-ready)
# timings; defaults to 1. Example: `julia manual_benchmarking.jl gpu --length 10`.
isnothing(ARGUMENTS["multiplier"]) || isnothing(ARGUMENTS["length"]) ||
    error("Give the multiplier either as `--length X` or as the second argument, not both.")
const TIMESTEP_MULTIPLIER = something(ARGUMENTS["length"], ARGUMENTS["multiplier"], 1.0)

# `--repeats N`: exactly N timed runs per configuration of the `BenchmarkSuite`s, without the time
# budget. `nothing` for the default: up to `repeats` runs within `max_repeat_time` (benchmark_suite.jl)
const REPEATS = ARGUMENTS["repeats"]
const SUITES_ARG = ARGUMENTS["suites"]

# Conditionally load backend packages BEFORE `using SpeedyWeather` so that
# SpeedyWeather's extensions register the right array types and overloads.
# CUDA must be loaded for any Reactant run (even the CPU, TPU and AMD GPU backends, where it
# loads without an NVIDIA GPU) — Reactant's KernelAbstractions kernel-raising path always needs it.
if ARCH_ARG == "gpu" || ON_REACTANT
    using CUDA
end
if ARCH_ARG == "amdgpu"
    using AMDGPU
end
if ON_REACTANT
    using Reactant
    Reactant.set_default_backend(REACTANT_BACKENDS[ARCH_ARG])       # "cpu", "gpu" or "tpu"
    # "gpu" is CUDA or ROCm, whichever Reactant_jll build is installed: check it matches the label
    if haskey(REACTANT_GPU_PLATFORMS, ARCH_ARG)
        platform = Reactant.XLA.platform_name(Reactant.XLA.default_backend())
        platform == REACTANT_GPU_PLATFORMS[ARCH_ARG] || error(
            "$ARCH_ARG needs Reactant's $(REACTANT_GPU_PLATFORMS[ARCH_ARG]) GPU client but got `$platform`. " *
                "Reactant_jll's GPU build is set by its `gpu` preference (cuda, rocm) or detected from the GPU driver."
        )
    end
end
if RUN_REVERSE
    using Enzyme   # Enzyme.autodiff and make_zero(::Variables) of SpeedyWeatherEnzymeExt
end

using SpeedyWeather, Dates, Printf, InteractiveUtils, JSON3
import SpeedyWeather.SpeedyTransforms: prettymemory

function pick_architecture(arg::AbstractString)
    if arg == "gpu"
        return (SpeedyWeather.GPU(), "gpu-nvidia")
    elseif arg == "amdgpu"
        return (SpeedyWeather.GPU(), "gpu-amd")
    elseif haskey(REACTANT_BACKENDS, arg)       # labelled as the argument
        return (SpeedyWeather.ReactantDevice(), arg)
    elseif arg == "cpu" || isempty(arg)
        arch_str = String(Sys.ARCH)
        label = (startswith(arch_str, "aarch") || arch_str == "arm64") ? "cpu-arm" : "cpu-x86"
        return (SpeedyWeather.CPU(), label)
    else
        error("Unknown architecture argument `$arg`. Use one of: $(join(ARCHITECTURES, ", ")).")
    end
end

const ARCH, ARCH_LABEL = pick_architecture(ARCH_ARG)
# GPU or TPU: higher resolutions (`--hi-res`, reverse mode), larger matrix transforms than on CPU (T128 on TPU)
const ON_ACCELERATOR = ARCH_LABEL in ("gpu-nvidia", "gpu-amd", "reactant-nvidia-gpu", "reactant-amd-gpu", "reactant-tpu")
@info "Running benchmarks for architecture: $ARCH_LABEL ($(typeof(ARCH)))"
RUN_HI_RES && !ON_ACCELERATOR && @warn "--hi-res only adds GPU and TPU configurations (T342, T512 are too slow on CPU), ignored on $ARCH_LABEL."

@info "Timestep multiplier (benchmark201 only, --length): $TIMESTEP_MULTIPLIER"
@info "Timed runs per configuration: $(isnothing(REPEATS) ? "default, up to 5 within 30 s, set exactly with --repeats N" : "$REPEATS (--repeats)")"
@info "Reverse-mode benchmarks: $(RUN_REVERSE ? "on" : "off, enable with --reverse")"
@info "High resolutions T342, T512 in benchmark201 (GPU, TPU only): $(RUN_HI_RES && ON_ACCELERATOR ? "on" : "off")"
@info "SPEEDY.f90-like configuration suite: $(RUN_SPEEDY ? "on" : "off, enable with --SPEEDY")"

# Free (device) memory between benchmark runs, called outside the timed region
free_memory!() = SpeedyWeather.SpeedyTransforms.clear_fourier_graph_cache!(ARCH)

# Per-architecture results store. Read it now, an unreadable store is an error before any suite
# runs (rather than after hours of benchmarking); it is read again for the merge at the end.
const ASSETS_DIR = joinpath(@__DIR__, "assets")
const RESULTS_JSON = joinpath(ASSETS_DIR, "benchmark_results.json")

# A store that exists but can't be parsed is an error, not a fresh start:
# writing a fresh store over it would lose the results of all architectures.
function load_results()
    isfile(RESULTS_JSON) || return Dict{String, Any}()
    try
        return Dict{String, Any}(JSON3.read(read(RESULTS_JSON, String), Dict{String, Any}))
    catch
        error("Could not parse $RESULTS_JSON. Fix it or remove it (which loses the stored results of all architectures).")
    end
end
load_results()

# The reverse-mode suite (`--reverse`), its numbers are also stored in the `reverse` field of the arch record
const REVERSE_SUITE_KEY = :benchmark800
# The SPEEDY.f90-like configuration suite (`--SPEEDY`), its numbers also in the `speedy` field
const SPEEDY_SUITE_KEY = :benchmark900

include("benchmark_suite.jl")
# The Reactant methods of the benchmark hooks need the `@compile` macro, only parse them once Reactant is loaded
ARCH isa SpeedyWeather.ReactantDevice && include("reactant_runner.jl")
RUN_REVERSE && include("reverse_benchmark_suite.jl")
include("define_benchmarks.jl")

for (key, suite) in benchmarks
    suite.architecture = ARCH
    adapt_to_architecture!(suite, ARCH)
    # The multiplier only lengthens the main benchmark, the PrimitiveWet resolution sweep
    key == :benchmark201 && (suite.timestep_multiplier = TIMESTEP_MULTIPLIER)
    # `--repeats N`: exactly N timed runs, no time budget (the reverse suite keeps its own samples)
    if !isnothing(REPEATS) && suite isa BenchmarkSuite
        suite.repeats = REPEATS
        suite.max_repeat_time = Inf
    end
    # AMDGPU's `default_gpu_graphs` is `false` (HIP graphs are opt-in,
    # so force it on explicitly for this arch;
    # leave every other arch on its own `SpectralTransform` default (`nothing`).
    hasproperty(suite, :gpu_graphs) && (suite.gpu_graphs = ARCH_LABEL == "gpu-amd" ? true : nothing)
end

# Deterministic suite order, so README sections are stable across runs. `--suites` narrows the run to
# a subset (in the same order); because each run replaces this architecture's whole stored record, a
# partial run drops the suites it did not run — warn when that is the case.
const ALL_SUITE_KEYS = sort!(collect(keys(benchmarks)))
const SUITE_KEYS = if isempty(strip(SUITES_ARG))
    ALL_SUITE_KEYS
else
    requested = [Symbol(strip(s)) for s in split(SUITES_ARG, ',') if !isempty(strip(s))]
    unknown = setdiff(requested, ALL_SUITE_KEYS)
    isempty(unknown) || error("Unknown suite(s): $(join(unknown, ", ")). Known suites: $(join(ALL_SUITE_KEYS, ", ")).")
    filter(in(requested), ALL_SUITE_KEYS)
end
SUITE_KEYS == ALL_SUITE_KEYS ||
    @warn "Running only $(join(SUITE_KEYS, ", ")): this replaces $ARCH_LABEL's stored results with just these suites."

for key in SUITE_KEYS
    @info "→ Running suite $key: $(benchmarks[key].title)"
    run_benchmark_suite!(benchmarks[key])
end


# The run's settings for repeats (`--repeats`) and the multiplier (`--length`), at the top of the section
function settings_markdown()
    suite = benchmarks[:benchmark201]
    runs = isnothing(REPEATS) ?
        "up to $(suite.repeats) timed runs per configuration, as long as they take at most $(round(Int, suite.max_repeat_time)) s together" :
        "$REPEATS timed run$(REPEATS == 1 ? "" : "s") per configuration (`--repeats $REPEATS`)"
    longer = TIMESTEP_MULTIPLIER == 1 ? "" :
        @sprintf(", the PrimitiveWet resolution sweep %g× longer (`--length %g`)", TIMESTEP_MULTIPLIER, TIMESTEP_MULTIPLIER)
    return "\nSettings: $runs$longer.\n"
end

function arch_markdown()
    io = IOBuffer()
    write(io, settings_markdown())
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
    elseif ARCH isa SpeedyWeather.GPU || ARCH_LABEL == "reactant-nvidia-gpu"
        write(io, "\n```julia\njulia> CUDA.versioninfo()\n")
        try
            CUDA.versioninfo(io)
        catch err
            write(io, "(CUDA.versioninfo() failed: $err)\n")
        end
        write(io, "```\n")
    end
    if ARCH isa SpeedyWeather.ReactantDevice
        backend = REACTANT_BACKENDS[ARCH_ARG]
        devices = Reactant.devices()
        write(io, "\nReactant backend: `$backend` (via `Reactant.set_default_backend(\"$backend\")`), ")
        write(io, "$(length(devices)) device(s): $(join(unique(Reactant.XLA.device_kind.(devices)), ", ")) (one used), ")
        write(io, "Reactant v$(pkgversion(Reactant)), Reactant_jll v$(pkgversion(Reactant.Reactant_jll)).\n")
    end
    return String(take!(io))
end

# JSON has no Inf/NaN literal — emit `null` for any non-finite metric (e.g. the SYPD of an
# unstable config) so JSON3 can write the file instead of erroring. The README generator
# already renders missing/non-finite entries as "—".
json_safe(x) = (x isa Number && isfinite(x)) ? x : nothing

# The structured numbers of a `BenchmarkSuite`: for the cross-architecture comparison from the
# PrimitiveWet resolution suite (:benchmark201, the `overview` field), and the SPEEDY suite (`speedy`)
function suite_data(suite::BenchmarkSuite)
    return Dict(
        "model" => string.(nameof.(suite.model)),
        "truncation" => collect(suite.truncation),
        "nlayers" => collect(suite.nlayers),
        "nlat" => collect(suite.nlat),
        "sypd" => map(json_safe, suite.SYPD),        # median of the timed runs
        "sypd_min" => map(json_safe, suite.SYPD_min),
        "sypd_max" => map(json_safe, suite.SYPD_max),
        "repeats" => collect(suite.nrepeats),
        "memory" => map(json_safe, suite.memory),
        "dt" => map(json_safe, suite.Δt),
        "compile_time" => map(json_safe, suite.compile_time),
        "spectral_transform" => string.(suite.spectral_transform),
    )
end

arch_record = Dict(
    "meta" => Dict(
        "speedyweather_version" => string(SpeedyWeather.pkgversion(SpeedyWeather)),
        "timestamp" => Dates.format(Dates.now(), Dates.RFC1123Format),
        "arch_type" => string(typeof(ARCH)),
        "repeats" => REPEATS,               # `nothing` (null) for the default, up to 5 within 30 s
        "length" => TIMESTEP_MULTIPLIER,
        "machine_info" => machine_info(),
    ),
    "markdown" => arch_markdown(),
)
# the structured numbers of the suites that ran (`--suites`)
:benchmark201 in SUITE_KEYS && (arch_record["overview"] = suite_data(benchmarks[:benchmark201]))
REVERSE_SUITE_KEY in SUITE_KEYS && (arch_record["reverse"] = reverse_data(benchmarks[REVERSE_SUITE_KEY]))
SPEEDY_SUITE_KEY in SUITE_KEYS && (arch_record["speedy"] = suite_data(benchmarks[SPEEDY_SUITE_KEY]))

# Merge into the JSON store, replacing everything stored for this label before
mkpath(ASSETS_DIR)
all_results = load_results()
all_results[ARCH_LABEL] = arch_record

open(RESULTS_JSON, "w") do io
    JSON3.pretty(io, all_results)
end
@info "Wrote results for $ARCH_LABEL → $RESULTS_JSON"

# Regenerate README.md from the JSON store

# Stable ordering for arch sections in the README.
const ARCH_ORDER = ["cpu-arm", "cpu-x86", "gpu-nvidia", "gpu-amd", "reactant-cpu", "reactant-nvidia-gpu", "reactant-amd-gpu", "reactant-tpu"]

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
    write(md, "By default each configuration is timed up to 5 times from the same initial state, as long as all its timed runs together ")
    write(md, "take at most 30 s, so short runs are timed 5 times and runs longer than 15 s once; `--repeats N` times every configuration ")
    write(md, "exactly N times. Each architecture's section states its settings. ")
    write(md, "SYPD is the median of these runs, and the Spread column shows how much they differ. ")
    write(md, "Results still vary between processes and machines, e.g. by ±30% on a laptop, ")
    write(md, "so the following presents a rough idea of how fast a SpeedyWeather simulation will run, and how much memory it requires.\n\n")

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
    write(md, "- SYPD: Speed of simulation, simulated years per wallclock day, median of the timed runs.\n")
    write(md, "- Runs: Number of timed runs, see the settings of each architecture.\n")
    write(md, "- Spread: (slowest − fastest) / median wallclock time of the timed runs, — for a single run.\n")
    write(md, "- Memory: Memory footprint of simulation, variables and constants.\n")
    write(md, "- Compile: One-off compilation of the time step (Reactant only), not included in the SYPD. ")
    write(md, "Compiled code is cached within a session, so configurations compiled before (e.g. in an earlier suite) compile much faster.\n")
    write(md, "- Steps, Forward, Gradient (reverse-mode benchmarks only): wallclock time of that many time steps forward, ")
    write(md, "and of the reverse-mode gradient of a scalar loss after them with respect to the initial state (incl. its forward pass), ")
    write(md, "with the Checkpointing scheme of the time loop. Gradient/Forward is their ratio.\n\n")

    write(md, "### Running the benchmarks\n\n")
    write(md, "Reproduce the benchmark suite by running, from `SpeedyWeather/benchmark`:\n\n")
    write(md, "```\n")
    write(md, "julia --project=. manual_benchmarking.jl                       # CPU (auto-labelled cpu-arm or cpu-x86)\n")
    write(md, "julia --project=. manual_benchmarking.jl gpu                   # CUDA GPU\n")
    write(md, "julia --project=. manual_benchmarking.jl amdgpu                # AMDGPU (HIP graphs forced on)\n")
    write(md, "julia --project=. manual_benchmarking.jl reactant-cpu          # Reactant on CPU\n")
    write(md, "julia --project=. manual_benchmarking.jl reactant-nvidia-gpu   # Reactant on NVIDIA GPU (CUDA)\n")
    write(md, "julia --project=. manual_benchmarking.jl reactant-amd-gpu      # Reactant on AMD GPU (ROCm)\n")
    write(md, "julia --project=. manual_benchmarking.jl reactant-tpu          # Reactant on TPU\n")
    write(md, "```\n\n")
    write(md, "Options and flags, in any order (`manual_benchmarking.jl --help` lists them): `--length X`, e.g. `manual_benchmarking.jl gpu --length 10`, multiplies the number of time steps of the ")
    write(md, "timed runs of the PrimitiveWet resolution sweep (the overview below) for more robust timings; the other suites are not affected. ")
    write(md, "`--repeats N` times every configuration exactly N times instead of up to 5 times within 30 s. ")
    write(md, "`--suites benchmark201,benchmark900` restricts the run to specific suites (default: all); ")
    write(md, "a partial run replaces that architecture's stored results with just the suites it ran. ")
    write(md, "`--reverse` (Reactant only, e.g. `manual_benchmarking.jl reactant-nvidia-gpu --reverse`) additionally runs the ")
    write(md, "reverse-mode benchmarks, which take hours. `--hi-res` (GPU and TPU only) adds T342 and T512 to the PrimitiveWet resolution sweep. ")
    write(md, "`--SPEEDY` additionally runs the PrimitiveWet model over the resolutions of the resolution sweep in a configuration closer to SPEEDY.f90 (seasonal prescribed climatologies of ")
    write(md, "sea surface temperature, land surface temperature and soil moisture, constant surface roughness, wind slowdown and gusts).\n\n")
    write(md, "Each run replaces its own architecture's section in this `README.md` (including suites added by flags in an earlier run); ")
    write(md, "results for other architectures are preserved via `benchmark_results.json`.\n\n")
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
