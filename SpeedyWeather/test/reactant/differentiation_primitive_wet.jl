"""
Differentiability of a PrimitiveWetModel time step with Reactant against Enzyme on the CPU, as
`differentation.jl` does for the BarotropicModel in CI (too slow for CI). Uses the functions defined
there. Enzyme on the CPU, the reference, currently needs Julia 1.10 or 1.11. Run from the repository root:

    julia --project=SpeedyWeather/test/reactant SpeedyWeather/test/reactant/differentiation_primitive_wet.jl

Reactant uses its default backend (a GPU if available), set e.g. `REACTANT_BACKEND=cpu` to choose.
"""

import Pkg
Pkg.activate(@__DIR__)

using SpeedyWeather
using Reactant
using Test
using CUDA      # needed to raise the KernelAbstractions kernels
using Enzyme

haskey(ENV, "REACTANT_BACKEND") && Reactant.set_default_backend(ENV["REACTANT_BACKEND"])

const TRUNCATION = 32       # spectral truncation, as in runtests.jl

include("setup.jl")

const DIFFERENTIATION_FUNCTIONS_ONLY = true     # skip the CI tests in differentation.jl
include("differentation.jl")

test_differentiability(PrimitiveWetModel)
