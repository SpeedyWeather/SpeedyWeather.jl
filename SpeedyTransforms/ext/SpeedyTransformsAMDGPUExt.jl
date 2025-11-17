module SpeedyTransformsAMDGPUExt
    
import AMDGPU: AMDGPU, ROCArray
import AbstractFFTs
import LinearAlgebra
using DocStringExtensions

using SpeedyTransforms
using SpeedyTransforms.RingGrids
using SpeedyTransforms.LowerTriangularArrays

# RocFFT is automatically chosen for ROCArray
SpeedyTransforms.which_FFT_library(::ROCArray) = AbstractFFTs
    
end 
