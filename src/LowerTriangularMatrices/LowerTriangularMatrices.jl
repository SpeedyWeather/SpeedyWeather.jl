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
export IndexBasis, OneBased, ZeroBased
export eachharmonic, eachmatrix, matrix_size

include("lower_triangular_matrix.jl")
include("plot.jl")

end