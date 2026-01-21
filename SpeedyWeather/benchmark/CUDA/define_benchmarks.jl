benchmarks = Dict{Symbol, AbstractBenchmarkSuiteTimed}()

# Individual transform benchmarks for CPU
benchmarks[:benchmark100] = BenchmarkSuiteTransform(
    title = "Transform benchmarks, CPU",
    nruns = 3,
    trunc = [31, 127, 511],
    nlayers = fill(64, 3),
    NF = fill(Float32, 3),
    model = fill(CPU(), 3),
)

benchmarks[:benchmark101] = BenchmarkSuiteTransform(
    title = "Transform benchmarks, GPU",
    nruns = 3,
    trunc = [31, 127, 511],
    nlayers = fill(64, 3),
    NF = fill(Float32, 3),
    model = fill(GPU(), 3),
)

benchmarks[:benchmark200] = BenchmarkSuiteModel(
    title = "PrimitiveWet benchmarks, CPU",
    nruns = 3,
    trunc = [31, 63, 127, 255, 511],
    nlayers = fill(16, 5),
    Grid = fill(SpeedyWeather.DEFAULT_GRID, 5),
    model = fill(CPU(), 5),
)

benchmarks[:benchmark201] = BenchmarkSuiteModel(
    title = "PrimitiveWet benchmarks, GPU",
    nruns = 3,
    trunc = [31, 63, 127, 255, 511],
    nlayers = fill(16, 5),
    Grid = fill(SpeedyWeather.DEFAULT_GRID, 5),
    model = fill(GPU(), 5),
)

benchmarks[:benchmark300] = BenchmarkSuiteDynamics(
    title = "PrimitiveWet dynamical core benchmarks, CPU",
    nruns = 3,
    trunc = [31, 63, 127, 255, 511, 1023, 1401], # up to 10km resolution
    nlayers = fill(16, 5),
    Grid = fill(SpeedyWeather.DEFAULT_GRID, 5),
    model = fill(CPU(), 5),
)

benchmarks[:benchmark301] = BenchmarkSuiteDynamics(
    title = "PrimitiveWet dynamical core benchmarks, GPU",
    nruns = 3,
    trunc = [31, 63, 127, 255, 511, 1023, 1401], # up to 10km resolution
    nlayers = fill(16, 5),
    Grid = fill(SpeedyWeather.DEFAULT_GRID, 5),
    model = fill(GPU(), 5),
)
