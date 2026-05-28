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
    trunc = [31, 42, 63, 85, 127, 170, 255],
)

## Primitive WET MODELS RESOLUTION
# Resolution sweep: L=8 from T31 → T255, plus L=16 and L=24 at the four highest
# truncations. Each configuration is run with the default (FFT + Legendre)
# SpectralTransform and additionally with MatrixSpectralTransform — except for
# T255 which is too large for the dense matrix transform (memory + speed).
let truncs   = [31, 42, 63, 85, 127, 170, 255, 85, 127, 170, 255, 85, 127, 170, 255],
    nlayers  = [ 8,  8,  8,  8,   8,   8,   8, 16, 16,  16,  16, 24, 24, 24, 24]

    matrix_idx = findall(t -> t < 150, truncs) # no matrix transform for T > 150 (it's too large)
    truncs_matrix  = truncs[matrix_idx]
    nlayers_matrix = nlayers[matrix_idx]

    n_default = length(truncs)
    n_matrix  = length(truncs_matrix)
    benchmarks[:benchmark201] = BenchmarkSuite(
        title = "Primitive wet model, resolution",
        nruns = n_default + n_matrix,
        model = fill(PrimitiveWetModel, n_default + n_matrix),
        trunc = vcat(truncs, truncs_matrix),
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
    trunc = fill(63, 6),
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
