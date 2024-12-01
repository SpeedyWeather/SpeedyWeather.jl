module LowerTriangularMatrices

# STRUCTURE
using DocStringExtensions

# GPU
import Adapt
import GPUArrays
import KernelAbstractions
import KernelAbstractions: get_backend

# NUMERICS
import LinearAlgebra: tril!

# VISUALISATION
import UnicodePlots

export LowerTriangularMatrix, LowerTriangularArray
export eachharmonic, eachmatrix
export OneBased, ZeroBased

include("lower_triangular_array.jl")
include("plot.jl")

end