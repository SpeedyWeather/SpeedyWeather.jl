module LowerTriangularMatrices

using DocStringExtensions
import Adapt, KernelAbstractions
import UnicodePlots
import LinearAlgebra: tril!
import GPUArrays
import Base: view

export LowerTriangularMatrix, LowerTriangularArray, eachharmonic
# export plot

include("lower_triangular_matrix.jl")
include("plot.jl")

end