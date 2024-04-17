module LowerTriangularMatrices

using DocStringExtensions
import Adapt, KernelAbstractions
import UnicodePlots
import LinearAlgebra: tril!
import GPUArrays

export LowerTriangularMatrix, LowerTriangularArray
export eachharmonic
# export plot

include("lower_triangular_matrix.jl")
include("plot.jl")

end