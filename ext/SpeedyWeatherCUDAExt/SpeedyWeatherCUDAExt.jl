module SpeedyWeatherCUDAExt

using SpeedyWeather
import CUDA: CUDA, CUDAKernels, CuArray, CUFFT
import SpeedyWeather.AbstractFFTs
using SpeedyWeather.DocStringExtensions
import SpeedyWeather: GPU, CPU

# DEVICE SETUP FOR CUDA
# extend functions from main SpeedyWeather 
 
# for RingGrids and LowerTriangularMatrices:
# every Array needs this method to strip away the parameters
RingGrids.nonparametric_type(::Type{<:CuArray}) = CuArray
LowerTriangularMatrices.nonparametric_type(::Type{<:CuArray}) = CuArray

SpeedyWeather.array_type(::GPU) = CuArray
SpeedyWeather.array_type(::Type{GPU}) = CuArray

SpeedyWeather.CUDAGPU() = GPU(CUDA.CUDABackend(always_inline=true))

const CUDAGPU = GPU{<:CUDA.CUDABackend}

function SpeedyWeather.GPU()
    if CUDA.has_cuda_gpu()
        return CUDAGPU()
    else
        msg = """We cannot make a GPU with the CUDA backend:
                 a CUDA GPU was not found!"""
        throw(ArgumentError(msg))
    end
end

SpeedyWeather.architecture(::CuArray) = CUDAGPU
SpeedyWeather.architecture(::Type{<:CuArray}) = CUDAGPU

SpeedyWeather.on_architecture(::CPU, a::CuArray) = Array(a)
SpeedyWeather.on_architecture(::CUDAGPU, a::CuArray) = a
SpeedyWeather.on_architecture(::CPU, a::SubArray{<:Any, <:Any, <:CuArray}) = Array(a)

SpeedyWeather.on_architecture(::CUDAGPU, a::Array) = CuArray(a)
SpeedyWeather.on_architecture(::CUDAGPU, a::BitArray) = CuArray(a)
SpeedyWeather.on_architecture(::CUDAGPU, a::SubArray{<:Any, <:Any, <:CuArray}) = a
SpeedyWeather.on_architecture(::CUDAGPU, a::SubArray{<:Any, <:Any, <:Array}) = CuArray(a)
SpeedyWeather.on_architecture(::CUDAGPU, a::StepRangeLen) = a

@inline SpeedyWeather.convert_to_device(::CUDAGPU, args) = CUDA.cudaconvert(args)
@inline SpeedyWeather.convert_to_device(::CUDAGPU, args::Tuple) = map(CUDA.cudaconvert, args)

include("spectral_transform.jl")
include("fourier.jl")
include("legendre.jl")

end # module