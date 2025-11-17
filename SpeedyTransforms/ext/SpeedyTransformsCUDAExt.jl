module SpeedyTransformsCUDAExt
    
import CUDA: CUDA, CUFFT, CuArray
import AbstractFFTs
import LinearAlgebra
using DocStringExtensions

using SpeedyTransforms
using SpeedyTransforms.RingGrids
using SpeedyTransforms.LowerTriangularArrays

# Use CUFFT for CuArrays
SpeedyTransforms.which_FFT_library(::CuArray) = CUFFT 
    
end 
