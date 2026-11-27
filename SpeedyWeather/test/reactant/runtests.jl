"""
Test scripts running SpeedyWeather with Reactant, in two groups that CI runs as separate jobs:

- `unit`: single features on Reactant (`@maybe_jit`, transforms, interpolation, output, restart files,
  compiled time loops, the implicit solver, differentiation)
- `correctness`: whole models on Reactant against the CPU (`test_correctness.jl`), their tendencies
  after one time step and their state after `NSTEPS` time steps

Run from the repository root with the groups to run as arguments, all groups without:

    julia --project=SpeedyWeather/test/reactant SpeedyWeather/test/reactant/runtests.jl unit
    julia --project=SpeedyWeather/test/reactant SpeedyWeather/test/reactant/runtests.jl correctness
"""

cd("SpeedyWeather/test/reactant")

const TEST_GROUPS = ("unit", "correctness")
const GROUPS = isempty(ARGS) ? TEST_GROUPS : Tuple(ARGS)
for group in GROUPS
    group in TEST_GROUPS || error("Unknown Reactant test group \"$group\", choose from $TEST_GROUPS")
end

import Pkg
Pkg.activate(@__DIR__)

using SpeedyWeather
using Reactant
using Test
using LinearAlgebra
using CUDA
using Statistics: mean
using Enzyme
import SpeedyWeather: ReactantDevice, time_step!, time_stepping!, @maybe_jit

const SpeedyWeatherReactantExt = Base.get_extension(SpeedyWeather, :SpeedyWeatherReactantExt)

# Configuration

const TRUNCATION = 32       # spectral truncation
const NSTEPS = 100          # number of time steps to compare
const RTOL = 1.0e-3         # relative tolerance for the comparisons in test_maybe_jit.jl
const ATOL = 1.0e-8         # absolute tolerance for the comparisons in test_maybe_jit.jl

# Tolerances of the CPU-vs-Reactant comparisons in test_correctness.jl, norm-relative for every variable,
# `norm(cpu - reactant) <= rtol * max(norm(cpu), norm(reactant))`, `default` for variables not listed.
# About 10x the largest differences measured (macOS arm64, Reactant 0.2.309, 2026-11-24, see
# docs/dev/2026-11/2026-11-24-reactant-default-primitive-wet-correctness.md). `tendencies`: one time
# step from identical states, rounding differences only (the spectral divergence, temperature and pressure
# tendencies are sums of terms that largely cancel). `time_stepping`: after NSTEPS time steps, rounding
# differences amplified by the dynamics, in the PrimitiveWetModel mostly by convection switching on or
# off at a few grid points
const ONE_LAYER_TOLERANCES = (
    tendencies = (default = 5.0e-5, divergence = 5.0e-4, η = 5.0e-4),
    time_stepping = (default = 5.0e-5, divergence = 1.0e-3),
)

const PRIMITIVE_TOLERANCES = (
    tendencies = (default = 1.0e-4, divergence = 2.0e-2, temperature = 3.0e-3, pressure = 2.0e-3),
    time_stepping = (
        default = 2.0e-2,       # vorticity, v, humidity
        divergence = 1.0e-1, temperature = 5.0e-4, pressure = 2.0e-5, u = 1.0e-2,
    ),
)

correctness_tolerances(::Type{<:Union{BarotropicModel, ShallowWaterModel}}) = ONE_LAYER_TOLERANCES
correctness_tolerances(::Type{<:Union{PrimitiveDryModel, PrimitiveWetModel}}) = PRIMITIVE_TOLERANCES

include("setup.jl")

if "unit" in GROUPS
    include("test_maybe_jit.jl")
    include("utilities.jl")
    include("test_transform.jl")
    include("test_interpolation.jl")
    include("test_output.jl")
    include("test_restart.jl")
    include("test_time_steps.jl")
    include("test_implicit.jl")
    include("test_drag.jl")
    include("test_thermodynamics.jl")
    include("differentation.jl")
end

if "correctness" in GROUPS
    include("test_correctness.jl")
end
