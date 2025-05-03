module LowerTriangularMatrices

# STRUCTURE
using DocStringExtensions

# GPU
import Adapt
import GPUArrays
import KernelAbstractions

# NUMERICS
import LinearAlgebra: tril!

# VISUALISATION
import UnicodePlots

export LowerTriangularMatrix, LowerTriangularArray
export eachharmonic, eachmatrix
export OneBased, ZeroBased

include("lower_triangular_array.jl")
include("rotate_reverse.jl")
include("plot.jl")

end