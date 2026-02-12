benchmarks = Dict{Symbol, AbstractBenchmarkSuiteTimed}()

# Individual transform benchmarks for CPU
#=
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
)=#

benchmarks[:benchmark200] = BenchmarkSuiteModel(
    title = "PrimitiveWet benchmarks, matrix transform CPU",
    nruns = 2,
    trunc = [63, 127],
    nlayers = fill(16, 5),
    Grid = fill(SpeedyWeather.DEFAULT_GRID, 5),
    transform_type = fill(:matrix, 5),
    model = fill(CPU(), 5),
)

benchmarks[:benchmark201] = BenchmarkSuiteModel(
    title = "PrimitiveWet benchmarks, matrix transform GPU",
    nruns = 2,
    trunc = [63, 127],
    nlayers = fill(16, 5),
    Grid = fill(SpeedyWeather.DEFAULT_GRID, 5),
    transform_type = fill(:matrix, 5),
    model = fill(GPU(), 5),
)


benchmarks[:benchmark202] = BenchmarkSuiteModel(
    title = "PrimitiveWet benchmarks, FFT+LT transform CPU",
    nruns = 2,
    trunc = [63, 127],
    nlayers = fill(16, 5),
    Grid = fill(SpeedyWeather.DEFAULT_GRID, 5),
    transform_type = fill(:fft, 5),
    model = fill(CPU(), 5),
)

benchmarks[:benchmark203] = BenchmarkSuiteModel(
    title = "PrimitiveWet benchmarks, FFT+LT transform GPU",
    nruns = 2,
    trunc = [63, 127],
    nlayers = fill(16, 5),
    Grid = fill(SpeedyWeather.DEFAULT_GRID, 5),
    transform_type = fill(:fft, 5),
    model = fill(GPU(), 5),
)


benchmarks[:benchmark300] = BenchmarkSuiteDynamicsGPU(
    title = "PrimitiveWet dynamical core benchmarks, matrix transform CPU",
    nruns = 2,
    trunc = [63, 127],
    nlayers = fill(16, 5),
    Grid = fill(SpeedyWeather.DEFAULT_GRID, 5),
    transform_type = fill(:matrix, 5),
    model = fill(CPU(), 5),
)

benchmarks[:benchmark301] = BenchmarkSuiteDynamicsGPU(
    title = "PrimitiveWet dynamical core benchmarks, matrix transform GPU",
    nruns = 2,
    trunc = [63, 127],
    nlayers = fill(16, 5),
    Grid = fill(SpeedyWeather.DEFAULT_GRID, 5),
    transform_type = fill(:matrix, 5),
    model = fill(GPU(), 5),
)

benchmarks[:benchmark302] = BenchmarkSuiteDynamicsGPU(
    title = "PrimitiveWet dynamical core benchmarks, FFT+LT transform CPU",
    nruns = 2,
    trunc = [63, 127],
    nlayers = fill(16, 5),
    Grid = fill(SpeedyWeather.DEFAULT_GRID, 5),
    transform_type = fill(:fft, 5),
    model = fill(CPU(), 5),
)

benchmarks[:benchmark303] = BenchmarkSuiteDynamicsGPU(
    title = "PrimitiveWet dynamical core benchmarks, FFT+LT transform GPU",
    nruns = 2,
    trunc = [63, 127],
    nlayers = fill(16, 5),
    Grid = fill(SpeedyWeather.DEFAULT_GRID, 5),
    transform_type = fill(:fft, 5),
    model = fill(GPU(), 5),
)
