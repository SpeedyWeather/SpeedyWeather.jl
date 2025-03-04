using CUDA, SpeedyWeather, Dates
using Adapt
using BenchmarkTools

include("benchmark_suite.jl")
include("define_benchmarks.jl")

## RUN BENCHMARKS
for key in keys(benchmarks)
    suite = benchmarks[key]
    run_benchmark_suite!(suite)
end

## WRITE into benchmarks.md
md = open("README.md", "w")
write(md, "# Benchmarks\n")

version = SpeedyWeather.pkgversion(SpeedyWeather)
write(md, "\ncreated for SpeedyWeather.jl v$version on $(Dates.format(Dates.now(), Dates.RFC1123Format)). \n\n")

write(md, "### Machine details\n\n")
write(md, "All benchmark simulatons were run on a server with access to an NVIDIA A100 GPU.\n")
write(md, "```julia\njulia> CUDA.versioninfo()\n")
CUDA.versioninfo(md)
write(md,"```\n\n")

write(md, "### Explanation\n\n")
write(md, "Abbreviations in the tables below are as follows, omitted columns use defaults.\n")
write(md, "- NF: Number format, default: $(SpeedyWeather.DEFAULT_NF)\n")
write(md, "- T: Spectral resolution, maximum degree of spherical harmonics, default: T$(SpeedyWeather.DEFAULT_TRUNC)\n")
write(md, "- L: Number of vertical layers, default: $(SpeedyWeather.DEFAULT_NLAYERS) (for 3D models)\n")
write(md, "- Grid: Horizontal grid, default: $(SpeedyWeather.DEFAULT_GRID)\n")
write(md, "- Rings: Grid-point resolution, number of latitude rings pole to pole\n")
write(md, "- Model: whether run on a CPU or GPU, default: true\n")



write(md, "### Running the benchmarks\n\n")
write(md, "The benchmark suite here can be reproduced by executing:\n\n")
write(md, "```> julia manual_benchmarking.jl```\n\n")
write(md, "inside `the SpeedyWeather.jl/benchmark` folder. It will create this `README.md` which can be pushed ")
write(md, "to the repository for updates or comparison.")

# Write benchmark suites into markdown
for key in keys(benchmarks)
    write_results(md, benchmarks[key])
end
write(md, "\n\n")

write(md, "## Benchmark graphs\n\n")
write(md, "Times and speedup graphs for the benchmarks are shown below.\n\n")
write(md, "![benchmark times](benchmark_times.png \"Times\") \n")
write(md, "![benchmark times](benchmark_speedup.png \"Speedup\") \n")

close(md)