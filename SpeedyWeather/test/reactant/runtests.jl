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
import SpeedyWeather: ReactantDevice, first_timesteps!, later_timestep!

# Configuration

const TRUNC = 31            # spectral truncation
const NSTEPS = 10           # number of time steps to compare
const RTOL = 1.0e-3           # relative tolerance for comparison
const ATOL = 1.0e-8          # absolute tolerance for comparison

include("setup.jl")

#include("test_correctness.jl")
include("differentation.jl")
