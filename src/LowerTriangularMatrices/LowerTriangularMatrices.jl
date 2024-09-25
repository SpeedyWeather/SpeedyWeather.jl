module LowerTriangularMatrices

# STRUCTURE
using DocStringExtensions

# GPU
import Adapt
import GPUArrays

# NUMERICS
import LinearAlgebra: tril!

# VISUALISATION
import UnicodePlots

export LowerTriangularMatrix, LowerTriangularArray
export eachharmonic, eachmatrix
export OneBased, ZeroBased

include("lower_triangular_matrix.jl")
include("plot.jl")

end