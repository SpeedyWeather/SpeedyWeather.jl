module LowerTriangularArrays

# STRUCTURE
using DocStringExtensions

# GPU
import Adapt
import GPUArrays
import KernelAbstractions

# NUMERICS
import LinearAlgebra: tril!

export AbstractSpectrum, Spectrum, resolution, truncation
 
export LowerTriangularMatrix, LowerTriangularArray
export eachharmonic, eachmatrix, eachorder
export OneBased, ZeroBased

include("spectrum.jl")
include("lower_triangular_array.jl")
include("rotate_reverse.jl")

end