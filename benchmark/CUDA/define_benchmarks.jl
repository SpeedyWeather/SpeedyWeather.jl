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