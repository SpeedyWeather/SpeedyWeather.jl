module SpeedyWeatherCUDAExt

using SpeedyWeather
import CUDA: CUDA, CUDAKernels, CuArray, CuDeviceArray, CUFFT
import AbstractFFTs
using DocStringExtensions
import SpeedyInternals.Architectures: Architectures, GPU, CPU, CUDAGPU, array_type, architecture, on_architecture, architecture, compatible_array_types, nonparametric_type

# DEVICE SETUP FOR CUDA
# extend functions from Architectures
Architectures.array_type(::GPU) = CuArray
Architectures.array_type(::Type{<:GPU}) = CuArray
Architectures.array_type(::GPU, NF::Type, N::Int) = CuArray{NF, N, CUDA.DeviceMemory}

Architectures.compatible_array_types(::GPU) = (CuArray, CuDeviceArray)
Architectures.compatible_array_types(::Type{<:GPU}) = (CuArray, CuDeviceArray)

Architectures.nonparametric_type(::Type{<:CuArray}) = CuArray

Architectures.CUDAGPU() = GPU(CUDA.CUDABackend(always_inline=true))
Architectures.GPU() = CUDAGPU() # default to CUDA

Architectures.architecture(::CuArray) = CUDAGPU()
Architectures.architecture(::Type{<:CuArray}) = CUDAGPU()
Architectures.architecture(::Type{<:CuDeviceArray}) = CUDAGPU()

Architectures.on_architecture(::CPU, a::CuArray) = Array(a)
Architectures.on_architecture(::GPU, a::CuArray) = a
Architectures.on_architecture(::CPU, a::SubArray{<:Any, <:Any, <:CuArray}) = Array(a)

Architectures.on_architecture(::GPU, a::Array) = CuArray(a)
Architectures.on_architecture(::GPU, a::BitArray) = CuArray(a)
Architectures.on_architecture(::GPU, a::SubArray{<:Any, <:Any, <:CuArray}) = a
Architectures.on_architecture(::GPU, a::SubArray{<:Any, <:Any, <:Array}) = CuArray(a)
Architectures.on_architecture(::GPU, a::StepRangeLen) = a

@inline Architectures.convert_to_device(::GPU, args) = CUDA.cudaconvert(args)
@inline Architectures.convert_to_device(::GPU, args::Tuple) = map(CUDA.cudaconvert, args)

include("fourier.jl")

end # module