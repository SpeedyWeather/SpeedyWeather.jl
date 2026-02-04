module LowerTriangularArrays

# STRUCTURE
using DocStringExtensions

# GPU
import Adapt: Adapt, adapt
import GPUArrays
import KernelAbstractions
import KernelAbstractions: @kernel, @index, @Const

import SpeedyWeatherInternals.Architectures: Architectures, AbstractArchitecture, on_architecture,
    array_type, ismatching, CPU, GPU, architecture, nonparametric_type

import SpeedyWeatherInternals.Utils: launch!, SpectralWorkOrder

# NUMERICS
import LinearAlgebra: Transpose, tril!

export AbstractSpectrum, Spectrum, resolution, truncation

export LowerTriangularMatrix, LowerTriangularArray
export eachharmonic, eachmatrix, eachorder, orders
export OneBased, ZeroBased
export lta_view
export zero_last_degree!

include("spectrum.jl")
include("lower_triangular_array.jl")
include("rotate_reverse.jl")
include("truncation.jl")

end
