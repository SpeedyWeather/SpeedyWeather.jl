# One run of the quadrature stability matrix, see
# docs/dev/2026-08/healpix-quadrature-exactness.md, "Verifying `PerOrderQuadrature` on a GPU node".
#
# Three things changed at once in this PR — the grid (dealiasing 3 -> 3.5), the quadrature weights,
# and the Nyquist bin — and the grid change is much the largest in model terms because it also
# changes the timestep through the CFL condition. The cases below separate them:
#
#   A  dealiasing 3.0   EqualAreaQuadrature    old grid, old weights (the status quo)
#   B  dealiasing 3.5   EqualAreaQuadrature    new grid, old weights
#   C  dealiasing 3.5   RingQuadrature         new grid, classical HEALPix weights (healpy/ducc/cuHPX)
#   D  dealiasing 3.5   PerOrderQuadrature     new grid, new weights (what this PR ships)
#   E  dealiasing 3.5   ContractiveQuadrature  new grid, non-expansive weights
#
# A->B isolates the resolution change, B->{C,D,E} isolates the weights. Dropping each ring's
# Nyquist bin follows the *grid* in the implementation rather than the scheme, so it is on in every
# case here and does not confound the comparison.
#
# Usage:
#   julia --project=SpeedyWeather/benchmark healpix_quadrature/run_case.jl \
#       case=D grid=HEALPixGrid trunc=128 diffusion_hours=4 years=10 output=runs/
#
# Writes <output>/<tag>.jld2 with the diagnostics time series and the run's outcome.

using CUDA
using Random
using SpeedyWeather
using SpeedyWeather.SpeedyTransforms
using JLD2
using Dates
using Printf

include(joinpath(@__DIR__, "diagnostics.jl"))

# ---------------------------------------------------------------------------------------------
# arguments

const ARGUMENTS = Dict{String, String}()
for arg in ARGS
    key, value = split(arg, '=', limit = 2)
    ARGUMENTS[String(key)] = String(value)
end

argument(key, default::AbstractString) = get(ARGUMENTS, key, default)
argument(key, default::Int) = parse(Int, get(ARGUMENTS, key, string(default)))
argument(key, default::Float64) = parse(Float64, get(ARGUMENTS, key, string(default)))

const CASES = Dict(
    "A" => (dealiasing = 3.0, Quadrature = SpeedyTransforms.EqualAreaQuadrature),
    "B" => (dealiasing = 3.5, Quadrature = SpeedyTransforms.EqualAreaQuadrature),
    "C" => (dealiasing = 3.5, Quadrature = SpeedyTransforms.RingQuadrature),
    "D" => (dealiasing = 3.5, Quadrature = SpeedyTransforms.PerOrderQuadrature),
    "E" => (dealiasing = 3.5, Quadrature = SpeedyTransforms.ContractiveQuadrature),
)

const GRIDS = Dict(
    "HEALPixGrid" => HEALPixGrid,
    "OctaHEALPixGrid" => OctaHEALPixGrid,
    "OctahedralGaussianGrid" => OctahedralGaussianGrid,
)

case = argument("case", "D")
Grid = GRIDS[argument("grid", "HEALPixGrid")]
truncation = argument("trunc", 128)
nlayers = argument("nlayers", 8)
diffusion_hours = argument("diffusion_hours", 4)
years = argument("years", 10)
days = argument("days", 0)          # non-zero overrides `years`, for calibration runs
output_directory = argument("output", joinpath(@__DIR__, "runs"))
record_every_days = argument("record_days", 1)
netcdf = argument("netcdf", "false") == "true"
# every Nth record, keep the surface vorticity coefficients so two runs on the same spectrum
# (B vs D, say) can be differenced directly; 0 disables
snapshot_every = argument("snapshot_every", 0)
# > 0 seeds a tiny perturbation of the initial vorticity, so the same case can be run as an
# ensemble. Single failure times are not comparable across cases — the configurations differ by far
# more than roundoff, so their trajectories decorrelate within weeks and a blow-up time is one draw
# from a distribution. Replicates are what make an ordering meaningful.
seed = argument("seed", 0)
perturbation = argument("perturbation", 1.0e-6)

haskey(CASES, case) || error("unknown case $case, expected one of $(sort(collect(keys(CASES))))")
(; dealiasing, Quadrature) = CASES[case]

period = days > 0 ? Day(days) : Year(years)
tag = days > 0 ?
    @sprintf("%s_%s_T%d_diff%dh_%dd", case, string(nameof(Grid)), truncation, diffusion_hours, days) :
    @sprintf("%s_%s_T%d_diff%dh_%dy", case, string(nameof(Grid)), truncation, diffusion_hours, years)
seed > 0 && (tag *= @sprintf("_seed%d", seed))
mkpath(output_directory)

# ---------------------------------------------------------------------------------------------
# model

architecture = SpeedyWeather.GPU()
spectral_grid = SpectralGrid(; truncation, nlayers, Grid, dealiasing, architecture)
spectral_transform = SpectralTransform(spectral_grid; Quadrature)

# Deliberately weakened diffusion is the sensitising knob: the useful signal is not "does it blow
# up" with the default settings but at what diffusion strength each scheme first fails.
horizontal_diffusion = HyperDiffusion(spectral_grid, time_scale = Hour(diffusion_hours))

model = if netcdf
    PrimitiveWetModel(
        spectral_grid; spectral_transform, horizontal_diffusion,
        output = NetCDFOutput(spectral_grid, path = output_directory, interval = Day(5)),
    )
else
    PrimitiveWetModel(spectral_grid; spectral_transform, horizontal_diffusion)
end

diagnostics = QuadratureDiagnostics(
    spectral_grid; schedule = Schedule(every = Day(record_every_days)), snapshot_every
)
add!(model, :quadrature_diagnostics => diagnostics)

println("=" ^ 90)
println("case $case: $(nameof(Grid)) T$truncation, dealiasing $dealiasing, $(nameof(Quadrature))")
println("  nlat_half   = $(spectral_grid.grid.nlat_half)")
println("  lmax, mmax  = $(spectral_grid.spectrum.lmax), $(spectral_grid.spectrum.mmax)")
println("  NF          = $(spectral_grid.NF)")
println("  diffusion   = $(diffusion_hours) h")
println("  period      = $period, recording every $record_every_days day(s)")
println("=" ^ 90)
flush(stdout)

simulation = initialize!(model)

if seed > 0
    # relative perturbation of the initial vorticity, one ensemble member per seed
    vorticity = simulation.variables.prognostic.vorticity
    host = on_architecture(SpeedyWeather.CPU(), vorticity)
    Random.seed!(seed)
    scale = maximum(abs, host.data)
    host.data .+= (perturbation * scale) .* randn(eltype(host.data), size(host.data))
    copyto!(vorticity.data, on_architecture(architecture, host).data)
    println("perturbed initial vorticity, seed $seed, relative amplitude $perturbation")
end

# ---------------------------------------------------------------------------------------------
# run

elapsed = @elapsed run!(simulation, period = period, output = netcdf)

# A blow-up does not stop the loop: `time_step!` returns immediately once NaNs are detected and
# the clock still runs to the end of the period. So the outcome comes from the feedback's flag,
# and the moment of failure from the first record that stopped being finite.
clock = simulation.variables.prognostic.clock
records = as_namedtuple(diagnostics)
finished = !model.feedback.nans_detected
first_nonfinite = findfirst(!isfinite, records.max_vorticity)
diverged_at = finished ? nothing :
    (isnothing(first_nonfinite) ? nothing : records.time[first_nonfinite])

jldsave(
    joinpath(output_directory, tag * ".jld2");
    case, grid = string(nameof(Grid)), truncation, nlayers, dealiasing,
    quadrature = string(nameof(Quadrature)), diffusion_hours, period = string(period),
    seed, perturbation,
    nlat_half = spectral_grid.grid.nlat_half, NF = string(spectral_grid.NF),
    timestep_seconds = Dates.value(Second(clock.Δt)),
    finished, diverged_at, elapsed, records...,
)

println()
println("case $case finished=$finished in $(round(elapsed, digits = 1)) s")
finished || println("  diverged, first non-finite record at $diverged_at")
@printf("  mean lnpₛ:  %.12f -> %.12f\n", first(records.mean_lnpressure), last(records.mean_lnpressure))
@printf("  energy:     %.6e -> %.6e\n", first(records.energy), last(records.energy))
@printf("  enstrophy:  %.6e -> %.6e\n", first(records.enstrophy), last(records.enstrophy))
@printf("  max |ζ|:    %.6e -> %.6e\n", first(records.max_vorticity), last(records.max_vorticity))
println("  wrote $(joinpath(output_directory, tag * ".jld2"))")
