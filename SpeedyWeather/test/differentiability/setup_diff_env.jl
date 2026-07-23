# One-time setup of the differentiability test environment.
# RUN THIS ONCE (on the x86 node, Julia 1.12) BEFORE launching the bisection sweep
# or the SLURM array — parallel array tasks must NOT all Pkg.instantiate at once
# (they would race on the shared Manifest / depot).
#
#   julia --project=SpeedyWeather/test/differentiability \
#         SpeedyWeather/test/differentiability/setup_diff_env.jl
#
# Repo root is assumed to be the current working directory.
import Pkg
Pkg.activate(@__DIR__)
Pkg.develop([
    Pkg.PackageSpec(path = joinpath(@__DIR__, "..", "..", "..", "SpeedyWeatherInternals")),
    Pkg.PackageSpec(path = joinpath(@__DIR__, "..", "..", "..", "LowerTriangularArrays")),
    Pkg.PackageSpec(path = joinpath(@__DIR__, "..", "..", "..", "RingGrids")),
    Pkg.PackageSpec(path = joinpath(@__DIR__, "..", "..", "..", "SpeedyTransforms")),
    Pkg.PackageSpec(path = joinpath(@__DIR__, "..", "..", "..", "SpeedyWeather")),
])
Pkg.instantiate()
Pkg.precompile()
using Enzyme
println("Setup done. Julia ", VERSION, "   Enzyme ", pkgversion(Enzyme))
