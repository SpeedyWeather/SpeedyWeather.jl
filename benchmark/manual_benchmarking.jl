using SpeedyWeather, Dates, Printf
import SpeedyWeather.SpeedyTransforms: prettymemory

# include("benchmark_suite.jl")
# include("define_benchmarks.jl")

# ## RUN BENCHMARKS
# for key in keys(benchmarks)
#     suite = benchmarks[key]
#     run_benchmark_suite!(suite)
# end

## WRITE into benchmarks.md
md = open("README.md", "w")
write(md, "# Benchmarks\n")

version = SpeedyWeather.pkgversion(SpeedyWeather)
write(md, "\ncreated for SpeedyWeather.jl v$version on $(Dates.format(Dates.now(), Dates.RFC1123Format)). \n\n")
write(md, "All simulations have been benchmarked over several seconds (wallclock time) without output. ")
write(md, "Benchmarking excludes initialization and is started just before the main time loop and finishes right after. ")
write(md, "All simulations single-threaded on a CPU, more details:\n")

write(md, "<details><summary>Machine details</summary>\n")
write(md, "```julia\njulia> versioninfo()\n")
versioninfo(md)
write(md,"```\n</details>")

write(md, "\n")
write(md, "The benchmarking results here are not very robust, timings that change with +-50% are not uncommon. ")
write(md, "Proper benchmarking for performance optimization uses the minimum or median of many executions, ")
write(md, "while we run a simulation for several time steps which effectively represents the mean, ")
write(md, "susceptible to outliers that slow down the simulation. However, this is what a user will experience ")
write(md, "in most situations anyway and the following therefore presents a rough idea of how fast a ")
write(md, "SpeedyWeather simulation will run, and how much memory it requires.\n\n")

write(md, "Explanation\n")
write(md, "- NF: Number format, default: $(SpeedyWeather.DEFAULT_NF)\n")
write(md, "- T: Spectral resolution, maximum degree of spherical harmonics, default: T$(SpeedyWeather.DEFAULT_TRUNC)\n")
write(md, "- L: Number of vertical layers, default: $(SpeedyWeather.DEFAULT_NLEV) (for 3D models)\n")
write(md, "- Grid: Horizontal grid, default: $(SpeedyWeather.DEFAULT_GRID)\n")
write(md, "- Rings: Grid-point resolution, number of latitude rings pole to pole\n")
write(md, "- Dynamics: With dynamics?, default: true\n")
write(md, "- Physics: With physical parameterizations?, default: true (for primitive equation models)\n")
write(md, "- Δt: time step [s].\n")
write(md, "- SYPD: Speed of simulation, simulated years per wallclock day.\n")
write(md, "- Memory: Memory footprint of simulation, variables and constants.\n")

# Write benchmark suites into markdown
for key in keys(benchmarks)

    suite = benchmarks[key]

    write(md, "\n## $(suite.title)\n\n")

    print_NF = any(suite.NF .!= suite.NF[1])
    print_Grid = any(suite.Grid .!= suite.Grid[1])
    print_nlat = any(suite.nlat .!= suite.nlat[1]) 
    print_dynamics = any(suite.dynamics .!= suite.dynamics[1])
    print_physics = any(suite.physics .!= suite.physics[1])

    column_header = "| Model "
    column_header *= print_NF ? "| NF " : ""
    column_header *= "| T "
    column_header *= "| L "
    column_header *= print_Grid ? "| Grid " : ""
    column_header *= print_nlat ? "| Rings " : ""
    column_header *= print_dynamics ? "| Dynamics " : ""
    column_header *= print_physics ? "| Physics " : ""
    column_header *= "| Δt | SYPD | Memory|"

    ncolumns = length(findall('|',column_header)) - 1
    second_row = repeat("| - ", ncolumns) * "|"

    write(md,"$column_header\n")
    write(md,"$second_row\n")

    for i in 1:suite.nruns

        row = "| $(suite.model[i]) "
        row *= print_NF ? "| $(suite.NF[i]) " : ""
        row *= "| $(suite.trunc[i]) "
        row *= "| $(suite.nlev[i]) "
        row *= print_Grid ? "| $(suite.Grid[i]) " : ""
        row *= print_nlat ? "| $(suite.nlat[i]) " : ""
        row *= print_dynamics ? "| $(suite.dynamics[i]) " : ""
        row *= print_physics  ? "| $(suite.physics[i]) " : ""

        Δt = round(Int,suite.Δt[i])
        sypd = suite.SYPD[i]
        SYPD = isfinite(sypd) ? round(Int, sypd) : 0
        memory = prettymemory(suite.memory[i])
        row *= "| $Δt | $SYPD | $memory |"

        write(md,"$row\n")
    end
end

close(md)