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
const RTOL = 1.0e-3         # default relative tolerance for comparison
const ATOL = 1.0e-8         # default absolute tolerance for comparison

const TOLERANCES = (
    default = (rtol = RTOL, atol = ATOL),
    temperature = (rtol = 2.0e-2, atol = 5.0),
    pressure = (rtol = 1.0e-2, atol = 5.0),
    vorticity = (rtol = 5.0e-2, atol = 1.0e-8),
    divergence = (rtol = 1.0e-1, atol = 1.0e-8),
    humidity = (rtol = 5.0e-2, atol = 2.0e-3),
    u = (rtol = 5.0e-2, atol = 1.0e-2),
    v = (rtol = 5.0e-2, atol = 1.0e-3),
    geopotential = (rtol = 1.0e-2, atol = 5.0e2),
    temp_average = (rtol = 1.0e-2, atol = 5.0e-1),
)

include("test_maybe_jit.jl")
include("utilities.jl")
include("setup.jl")
include("test_transform.jl")
include("test_interpolation.jl")
include("test_correctness.jl")
include("test_output.jl")
include("test_time_steps.jl")
include("test_implicit.jl")
include("test_drag.jl")

include("differentation.jl")
