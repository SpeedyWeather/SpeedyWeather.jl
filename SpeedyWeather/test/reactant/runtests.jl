"""
Test scripts running SpeedyWeather with Reactant.
"""

cd("SpeedyWeather/test/reactant")
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
# off at a few grid points (a 1-ulp perturbation of the CPU model grows to the same differences)
const ONE_LAYER_TOLERANCES = (
    tendencies = (default = 5.0e-5, divergence = 5.0e-4, η = 5.0e-4),
    time_stepping = (default = 5.0e-5, divergence = 1.0e-3),
)

const PRIMITIVE_TENDENCY_TOLERANCES = (default = 1.0e-4, divergence = 1.0e-3, temperature = 2.0e-4, pressure = 5.0e-4)

const PRIMITIVE_DRY_TOLERANCES = (
    tendencies = PRIMITIVE_TENDENCY_TOLERANCES,
    time_stepping = (default = 2.0e-4, divergence = 2.0e-3, v = 5.0e-4, temperature = 1.0e-5, pressure = 3.0e-6),
)

const PRIMITIVE_WET_TOLERANCES = (
    tendencies = PRIMITIVE_TENDENCY_TOLERANCES,
    time_stepping = (
        default = 2.0e-2,       # vorticity, v, humidity
        divergence = 1.0e-1, temperature = 5.0e-4, pressure = 2.0e-5, u = 1.0e-2,
    ),
)

correctness_tolerances(::Type{<:Union{BarotropicModel, ShallowWaterModel}}) = ONE_LAYER_TOLERANCES
correctness_tolerances(::Type{<:PrimitiveDryModel}) = PRIMITIVE_DRY_TOLERANCES
correctness_tolerances(::Type{<:PrimitiveWetModel}) = PRIMITIVE_WET_TOLERANCES

include("test_maybe_jit.jl")
include("utilities.jl")
include("setup.jl")
include("test_transform.jl")
include("test_interpolation.jl")
include("test_correctness.jl")
include("test_output.jl")
include("test_restart.jl")
include("test_time_steps.jl")
include("test_implicit.jl")
include("test_drag.jl")

include("differentation.jl")
