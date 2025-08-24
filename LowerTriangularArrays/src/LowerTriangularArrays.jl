module LowerTriangularArrays

# STRUCTURE
using DocStringExtensions

# GPU
import Adapt: Adapt, adapt
import GPUArrays
import KernelAbstractions
 
import SpeedyInternals.Architectures: Architectures, AbstractArchitecture, on_architecture, 
    array_type, ismatching, CPU, GPU, architecture, nonparametric_type

# NUMERICS
import LinearAlgebra: tril!

export AbstractSpectrum, Spectrum, resolution, truncation
 
export LowerTriangularMatrix, LowerTriangularArray
export eachharmonic, eachmatrix, eachorder, orders
export OneBased, ZeroBased
export lta_view

include("spectrum.jl")
include("lower_triangular_array.jl")
include("rotate_reverse.jl")

end