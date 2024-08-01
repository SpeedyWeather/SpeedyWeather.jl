"""
    abstract type AbstractDevice 

Supertype of all devices SpeedyWeather.jl can run on
"""
abstract type AbstractDevice end 

export CPU, GPU

"""
    CPU <: AbstractDevice

Indicates that SpeedyWeather.jl runs on a single CPU 
"""
struct CPU <: AbstractDevice end

"""
    GPU <: AbstractDevice

Indicates that SpeedyWeather.jl runs on a single GPU
"""
struct GPU <: AbstractDevice end 

"""$(TYPEDSIGNATURES)
Return default used device for internal purposes, either `CPU` or `GPU` if a GPU is available.
"""
Device() = CUDA.functional() ? GPU() : CPU()

"""$(TYPEDSIGNATURES)
Default array type on `device`."""
default_array_type(device::AbstractDevice) = default_array_type(typeof(device))
default_array_type(::Type{GPU}) = CuArray
default_array_type(::Type{CPU}) = Array

"""$(TYPEDSIGNATURES)
Return default used device for KernelAbstractions, either `CPU` or `CUDADevice` if a GPU is available
"""
Device_KernelAbstractions() = CUDA.functional() ? KernelAbstractions.CUDADevice : KernelAbstractions.CPU

"""$(TYPEDSIGNATURES)
Return used device for use with KernelAbstractions
"""
Device_KernelAbstractions(::CPU) = KernelAbstractions.CPU
Device_KernelAbstractions(::GPU) = KernelAbstractions.CUDADevice

"""$(TYPEDSIGNATURES)
Holds information about the device the model is running on and workgroup size. 
$(TYPEDFIELDS)"""
struct DeviceSetup{S<:AbstractDevice, T}
    "::AbstractDevice, device the model is running on."
    device::S
    
    "::KernelAbstractions.Device, device for use with KernelAbstractions"
    device_KA::T
    
    "workgroup size" 
    n::Int
end 

DeviceSetup() = DeviceSetup(Device(), Device_KernelAbstractions(Device()), workgroup_size(Device()))
DeviceSetup(device::AbstractDevice) = DeviceSetup(device, Device_KernelAbstractions(device), workgroup_size(device))
DeviceSetup(device::AbstractDevice, n::Integer) = DeviceSetup(device, Device_KernelAbstractions(device), n)

"""$(TYPEDSIGNATURES)
Returns a workgroup size depending on `device`. 
WIP: Will be expanded in the future to also include grid information. 
"""
function workgroup_size(device::AbstractDevice)
    return device isa GPU ? 32 : 4 
end

"""$(TYPEDSIGNATURES)
Adapts `x` to a `CuArray` when `device::GPU` is used, otherwise a regular `Array`.
Uses `adapt`, thus also can return SubArrays etc."""
DeviceArray(::GPU, x) = Adapt.adapt(CuArray, x)
DeviceArray(::CPU, x) = Adapt.adapt(Array, x)
DeviceArray(dev::DeviceSetup, x) = DeviceArray(dev.device, x)

"""$(TYPEDSIGNATURES)
Returns a `CuArray` when `device<:GPU` is used, otherwise a regular `Array`.
Doesn't uses `adapt`, therefore always returns CuArray/Array."""
DeviceArrayNotAdapt(::GPU, x) = CuArray(x)
DeviceArrayNotAdapt(::CPU, x) = Array(x)
DeviceArrayNotAdapt(dev::DeviceSetup, x) = DeviceArrayNotAdapt(dev.device, x)

"""$(TYPEDSIGNATURES)
Launches the `kernel!` on the `device_setup` with `ndrange` computations over the
kernel and arguments `kernel_args`."""
function launch_kernel!(device_setup::DeviceSetup, kernel!, ndrange, kernel_args...)
    device = device_setup.device_KA()
    n = device_setup.n 

    k! = kernel!(device, n)
    k!(kernel_args...; ndrange=ndrange)

    return nothing
end