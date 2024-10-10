module SpeedyWeatherCUDAExt

using SpeedyWeather
import CUDA: CUDA, CUDAKernels, CuArray

# for RingGrids and LowerTriangularMatrices:
# every Array needs this method to strip away the parameters
RingGrids.nonparametric_type(::Type{<:CuArray}) = CuArray
LowerTriangularMatrices.nonparametric_type(::Type{<:CuArray}) = CuArray

SpeedyWeather.default_array_type(::Type{GPU}) = CuArray

SpeedyWeather.DeviceSetup() = DeviceSetup(Device(), Device_KernelAbstractions(Device()), workgroup_size(Device()))

# DEVICE SETUP FOR CUDA

"""$(TYPEDSIGNATURES)
Return default used device for internal purposes, either `CPU` or `GPU` if a GPU is available."""
SpeedyWeather.Device() = CUDA.functional() ? GPU() : CPU()

"""$(TYPEDSIGNATURES)
Return default used device for KernelAbstractions, either `CPU` or `CUDADevice` if a GPU is available."""
SpeedyWeather.Device_KernelAbstractions() = CUDA.functional() ? KernelAbstractions.CUDADevice : KernelAbstractions.CPU
SpeedyWeather.Device_KernelAbstractions(::GPU) = KernelAbstractions.CUDADevice

SpeedyWeahter.DeviceArray(::GPU, x) = Adapt.adapt(CuArray, x)

"""$(TYPEDSIGNATURES)
Returns a `CuArray` when `device<:GPU` is used. Doesn't uses `adapt`, therefore always returns CuArray."""
SpeedyWeather.DeviceArrayNotAdapt(::GPU, x) = CuArray(x)

end # module