module SpeedyWeatherCUDAExt

using SpeedyWeather
import CUDA: CUDA, CUDAKernels, CuArray, CuDeviceArray, CUFFT
import AbstractFFTs
using SpeedyWeather.DocStringExtensions
import SpeedyWeather: GPU, CPU, CUDAGPU, array_type, architecture, on_architecture, architecture

# DEVICE SETUP FOR CUDA
# extend functions from main SpeedyWeather 

# for RingGrids and LowerTriangularArrays:
# every Array needs this method to strip away the parameters
RingGrids.nonparametric_type(::Type{<:CuArray}) = CuArray
LowerTriangularArrays.nonparametric_type(::Type{<:CuArray}) = CuArray

array_type(::GPU) = CuArray
array_type(::Type{<:GPU}) = CuArray

CUDAGPU() = GPU(CUDA.CUDABackend(always_inline=true))
GPU() = CUDAGPU() # default to CUDA

architecture(::CuArray) = CUDAGPU()
architecture(::Type{<:CuArray}) = CUDAGPU()
architecture(::Type{<:CuDeviceArray}) = CUDAGPU()

on_architecture(::CPU, a::CuArray) = Array(a)
on_architecture(::GPU, a::CuArray) = a
on_architecture(::CPU, a::SubArray{<:Any, <:Any, <:CuArray}) = Array(a)

on_architecture(::GPU, a::Array) = CuArray(a)
on_architecture(::GPU, a::BitArray) = CuArray(a)
on_architecture(::GPU, a::SubArray{<:Any, <:Any, <:CuArray}) = a
on_architecture(::GPU, a::SubArray{<:Any, <:Any, <:Array}) = CuArray(a)
on_architecture(::GPU, a::StepRangeLen) = a

@inline SpeedyWeather.convert_to_device(::GPU, args) = CUDA.cudaconvert(args)
@inline SpeedyWeather.convert_to_device(::GPU, args::Tuple) = map(CUDA.cudaconvert, args)

include("spectral_transform.jl")
include("fourier.jl")
include("legendre.jl")

end # module