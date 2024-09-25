# dictionary of all benchmark suites, define with whatever key ::Symbol
benchmarks = Dict{Symbol,AbstractBenchmarkSuite}()

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
benchmarks[:benchmark201] = BenchmarkSuite(
    title = "Primitive wet model, resolution",
    nruns = 6,
    model = fill(PrimitiveWetModel, 6),
    trunc = [31, 42, 63, 85, 127, 170],
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
    Grid = [FullGaussianGrid, FullClenshawGrid, OctahedralGaussianGrid, OctahedralClenshawGrid,
            HEALPixGrid, OctaHEALPixGrid],
    )

## NLEV
benchmarks[:benchmark500] = BenchmarkSuite(
    title = "Number of vertical layers",
    nruns = 4,
    nlev = [4, 8, 12, 16],
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
    physics = [true],
    dynamics = [true],
)