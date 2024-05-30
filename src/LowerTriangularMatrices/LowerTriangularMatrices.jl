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
# export plot

export LowerTriangularMatrix, LowerTriangularArray
export eachharmonic, eachmatrix, add!

include("lower_triangular_matrix.jl")
include("plot.jl")

end