# dictionary of all benchmark suites, define with whatever key ::Symbol
benchmarks = Dict{Symbol, AbstractBenchmarkSuite}()

# Models
benchmarks[:benchmark100] = BenchmarkSuite(
    title = "Models, default setups",
    nruns = 4,
    model = [BarotropicModel, ShallowWaterModel, PrimitiveDryModel, PrimitiveWetModel],
)

# BarotropicModel, resolution
benchmarks[:benchmark200] = BenchmarkSuite(
    title = "Shallow water model, resolution",
    nruns = 7,
    model = fill(ShallowWaterModel, 7),
    truncation = [32, 43, 64, 86, 128, 171, 256],
)

## Primitive WET MODELS RESOLUTION
# Resolution sweep: L=8 from T32 → T256, plus L=16 and L=24 at the four highest
# truncations. With `--hi-res` on GPU and TPU additionally T342 and T512 for all three layer
# counts (not on CPU, too slow). Each configuration is run with the default (FFT + Legendre)
# SpectralTransform and additionally with MatrixSpectralTransform up to an
# architecture-dependent maximum truncation, as the dense matrix transform gets
# too large (memory + speed) at high resolution: up to T171 on GPU and TPU (incl. Reactant GPU),
# up to T86 on CPU (incl. Reactant CPU). `ARCH_LABEL`, `ON_ACCELERATOR` and `RUN_HI_RES` are set in manual_benchmarking.jl.
let hi_res = RUN_HI_RES && ON_ACCELERATOR ? [342, 512] : Int[],
        truncations_L8 = vcat([32, 43, 64, 86, 128, 171, 256], hi_res),
        truncations_L16_L24 = vcat([86, 128, 171, 256], hi_res),
        truncations = vcat(truncations_L8, truncations_L16_L24, truncations_L16_L24),
        nlayers = vcat(fill(8, length(truncations_L8)), fill(16, length(truncations_L16_L24)), fill(24, length(truncations_L16_L24)))

    max_truncation_matrix = ON_ACCELERATOR ? 180 : 86
    matrix_idx = findall(t -> t <= max_truncation_matrix, truncations)
    truncations_matrix = truncations[matrix_idx]
    nlayers_matrix = nlayers[matrix_idx]

    n_default = length(truncations)
    n_matrix = length(truncations_matrix)
    benchmarks[:benchmark201] = BenchmarkSuite(
        title = "Primitive wet model, resolution",
        nruns = n_default + n_matrix,
        model = fill(PrimitiveWetModel, n_default + n_matrix),
        truncation = vcat(truncations, truncations_matrix),
        nlayers = vcat(nlayers, nlayers_matrix),
        spectral_transform = vcat(fill(:default, n_default), fill(:matrix, n_matrix)),
    )
end

## NUMBER FORMATS
benchmarks[:benchmark300] = BenchmarkSuite(
    title = "Primitive Equation, Float32 vs Float64",
    nruns = 2,
    NF = [Float32, Float64],
)

## GRIDS
benchmarks[:benchmark400] = BenchmarkSuite(
    title = "Grids",
    nruns = 6,
    truncation = fill(64, 6),
    Grid = [
        FullGaussianGrid, FullClenshawGrid, OctahedralGaussianGrid, OctahedralClenshawGrid,
        HEALPixGrid, OctaHEALPixGrid,
    ],
)

## nlayers
benchmarks[:benchmark500] = BenchmarkSuite(
    title = "Number of vertical layers",
    nruns = 4,
    nlayers = [4, 8, 12, 16],
)

## PHYSICS/DYNAMICS
benchmarks[:benchmark600] = BenchmarkSuite(
    title = "PrimitiveDryModel: Physics or dynamics only",
    nruns = 3,
    model = [PrimitiveDryModel, PrimitiveDryModel, PrimitiveDryModel],
    physics = [true, false, true],
    dynamics = [true, true, false],
)

## PHYSICS/DYNAMICS
benchmarks[:benchmark601] = BenchmarkSuite(
    title = "PrimitiveWetModel: Physics or dynamics only",
    nruns = 3,
    physics = [true, false, true],
    dynamics = [true, true, false],
)

## DYNAMICS, benchmark individual functions
benchmarks[:benchmark700] = BenchmarkSuiteDynamics(
    title = "Individual dynamics functions",
    nruns = 1,
)

## REVERSE MODE (AD) WITH REACTANT, opt-in with `--reverse` (reverse_benchmark_suite.jl)
# Forward and reverse-mode gradient over 50, 100, 200, 2000 time steps of the PrimitiveWet model, L8,
# SpectralTransform. On Reactant CPU only up to T128, above the 2000-step gradient takes too long.
if RUN_REVERSE
    let truncations = ON_ACCELERATOR ? [32, 64, 128, 256] : [32, 64, 128],
            nsteps = [50, 100, 200, 2000]

        benchmarks[REVERSE_SUITE_KEY] = BenchmarkSuiteReverse(
            title = "Primitive wet model, reverse mode (AD) with Reactant",
            nruns = length(truncations) * length(nsteps),
            truncation = repeat(truncations, inner = length(nsteps)),
            nsteps = repeat(nsteps, outer = length(truncations)),
        )
    end
end
