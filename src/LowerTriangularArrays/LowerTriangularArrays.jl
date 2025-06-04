module LowerTriangularArrays

# STRUCTURE
using DocStringExtensions

# GPU
import Adapt
import Adapt: adapt
import GPUArrays
import KernelAbstractions
import ..Architectures: AbstractArchitecture, on_architecture, array_type, ismatching, CPU, architecture

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