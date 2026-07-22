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
using KernelAbstractions
import SpeedyWeather: ReactantDevice, first_timesteps!, later_timestep!, @maybe_jit

const SpeedyWeatherReactantExt = Base.get_extension(SpeedyWeather, :SpeedyWeatherReactantExt)
Reactant.set_default_backend("gpu")

# Configuration

const TRUNC = 31            # spectral truncation
const NSTEPS = 10           # number of time steps to compare
const RTOL = 1.0e-3         # relative tolerance for comparison
const ATOL = 1.0e-8         # absolute tolerance for comparison


nlayers = 8
arch = SpeedyWeather.ReactantDevice()

spectral_grid = SpectralGrid(; architecture = arch, nlayers, trunc = TRUNC)
M = MatrixSpectralTransform(spectral_grid)
initial_conditions = InitialConditions(spectral_grid, PrimitiveWetModel)
longwave_radiation = OneBandLongwave(spectral_grid, transmissivity = ConstantLongwaveTransmissivity(spectral_grid))
model = PrimitiveWetModel(spectral_grid; spectral_transform = M, feedback = nothing, longwave_radiation, initial_conditions, convection=nothing, output=nothing)
simulation = initialize!(model)

# full model run, stops at the problematic kernel
#run!(simulation, steps=10)

# Just the problematic kernel by itself 
initialize!(simulation, steps=10)
SpeedyWeather.column_parameterizations!(simulation.variables, simulation.model)

