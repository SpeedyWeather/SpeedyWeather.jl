module LowerTriangularMatrices

# STRUCTURE
using DocStringExtensions

# GPU
import Adapt
import GPUArrays
import KernelAbstractions
import ..Architectures: AbstractArchitecture, on_architecture, array_type

# NUMERICS
import LinearAlgebra: tril!

export LowerTriangularMatrix, LowerTriangularArray
export eachharmonic, eachmatrix
export OneBased, ZeroBased

include("lower_triangular_array.jl")
include("rotate_reverse.jl")

end