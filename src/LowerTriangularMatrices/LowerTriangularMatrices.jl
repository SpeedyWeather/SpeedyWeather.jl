module LowerTriangularMatrices

using DocStringExtensions
import Adapt
import UnicodePlots

export LowerTriangularMatrix, LowerTriangularArray, eachharmonic
# export plot

include("lower_triangular_matrix.jl")
include("plot.jl")

end