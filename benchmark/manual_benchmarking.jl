using SpeedyWeather, Dates, Printf
import SpeedyWeather.SpeedyTransforms: prettymemory

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
write(md, "All simulations have been benchmarked over several seconds (wallclock time) without output. ")
write(md, "Benchmarking excludes initialization and is started just before the main time loop and finishes right after. ")
write(md, "The benchmarking results here are not very robust, timings that change with +-50% are not uncommon. ")
write(md, "Proper benchmarking for performance optimization uses the minimum or median of many executions, ")
write(md, "while we run a simulation for several time steps which effectively represents the mean, ")
write(md, "susceptible to outliers that slow down the simulation. However, this is what a user will experience ")
write(md, "in most situations anyway and the following therefore presents a rough idea of how fast a ")
write(md, "SpeedyWeather simulation will run, and how much memory it requires.\n\n")

write(md, "### Machine details\n\n")
write(md, "All benchmark simulation were single-threaded on a CPU:\n")

write(md, "```julia\njulia> versioninfo()\n")
versioninfo(md)
write(md,"```\n\n")

write(md, "### Explanation\n\n")
write(md, "Abbreviations in the tables below are as follows, omitted columns use defaults.\n")
write(md, "- NF: Number format, default: $(SpeedyWeather.DEFAULT_NF)\n")
write(md, "- T: Spectral resolution, maximum degree of spherical harmonics, default: T$(SpeedyWeather.DEFAULT_TRUNC)\n")
write(md, "- L: Number of vertical layers, default: $(SpeedyWeather.DEFAULT_NLAYERS) (for 3D models)\n")
write(md, "- Grid: Horizontal grid, default: $(SpeedyWeather.DEFAULT_GRID)\n")
write(md, "- Rings: Grid-point resolution, number of latitude rings pole to pole\n")
write(md, "- Dynamics: With dynamics?, default: true\n")
write(md, "- Physics: With physical parameterizations?, default: true (for primitive equation models)\n")
write(md, "- Δt: time step [s].\n")
write(md, "- SYPD: Speed of simulation, simulated years per wallclock day.\n")
write(md, "- Memory: Memory footprint of simulation, variables and constants.\n\n")

write(md, "### Running the benchmarks\n\n")
write(md, "The benchmark suite here can be reproduced by executing:\n\n")
write(md, "```> julia manual_benchmarking.jl```\n\n")
write(md, "inside `the SpeedyWeather.jl/benchmark` folder. It will create this `README.md` which can be pushed ")
write(md, "to the repository for updates or comparison.")

# Write benchmark suites into markdown
for key in keys(benchmarks)
    write_results(md, benchmarks[key])
end

close(md)