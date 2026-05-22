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
# Resolution sweep: L=8 from T31 → T255, plus L=16 at the three highest truncations
# so that the L=8 vs L=16 comparison plots have enough overlap.
benchmarks[:benchmark201] = BenchmarkSuite(
    title = "Primitive wet model, resolution",
    nruns = 11,
    model = fill(PrimitiveWetModel, 11),
    trunc =   [31, 42, 63, 85, 127, 170, 255, 85, 127, 170, 255],
    nlayers = [ 8,  8,  8,  8,   8,   8,   8, 16, 16,  16,  16],
)

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
