"""
    abstract type AbstractDevice 

Supertype of all devices SpeedyWeather.jl can ran on
"""
abstract type AbstractDevice end 

"""
    CPUDevice <: AbstractDevice

Indicates that SpeedyWeather.jl runs on a single CPU 
"""
struct CPUDevice <: AbstractDevice end

"""
    GPUDevice <: AbstractDevice

Indicates that SpeedyWeather.jl runs on a single GPU
"""
struct GPUDevice <: AbstractDevice end 

"""
    Device()

Return default used device for internal purposes, either `CPUDevice` or `GPUDevice` if a GPU is available.
"""
Device() = CUDA.functional() ? GPUDevice() : CPUDevice()
      
"""
    Device_KernelAbstractions()

Return default used device for KernelAbstractions, either `CPU` or `CUDADevice` if a GPU is available
"""
Device_KernelAbstractions() = CUDA.functional() ? KernelAbstractions.CUDADevice : KernelAbstractions.CPU

"""
    Device_KernelAbstractions(::AbstractDevice)

Return used device for use with KernelAbstractions
"""
Device_KernelAbstractions(::CPUDevice) = KernelAbstractions.CPU
Device_KernelAbstractions(::GPUDevice) = KernelAbstractions.CUDADevice

"""
    DeviceSetup{S<:AbstractDevice}

Holds information about the device the model is running on and workgroup size. 

* `device::AbstractDevice`: Device the model is running on 
* `device_KA::KernelAbstractions.Device`: Device for use with KernelAbstractions
* `n`: workgroup size  
"""
struct DeviceSetup{S<:AbstractDevice,T}
    device::S # for internal purposes
    device_KA::T # for KernelAbstractions
    n::Int # workgroup size
end 

DeviceSetup() = DeviceSetup(Device(), Device_KernelAbstractions(Device()), workgroup_size(Device()))
DeviceSetup(device::AbstractDevice) = DeviceSetup(device, Device_KernelAbstractions(device), workgroup_size(device))
DeviceSetup(device::AbstractDevice, n::Integer) = DeviceSetup(device, Device_KernelAbstractions(device), n)

"""
    workgroup_size(dev::AbstractDevice)

Returns a workgroup size depending on `dev`. 
WIP: Will be expanded in the future to also include grid information. 
"""
function workgroup_size(device::AbstractDevice)
    return device isa GPUDevice ? 32 : 4 
end

"""
    DeviceArray(device::AbstractDevice, x) 

Adapts `x` to a `CuArray` when `device<:GPUDevice` is used, otherwise a regular `Array`. Uses `adapt`, thus also can return SubArrays etc.
"""
DeviceArray(::GPUDevice, x) = Adapt.adapt(CuArray, x)
DeviceArray(::CPUDevice, x) = Adapt.adapt(Array, x)
DeviceArray(dev::DeviceSetup, x) = DeviceArray(dev.device, x)

"""
    DeviceArrayNotAdapt(device::AbstractDevice, x) 

Returns a `CuArray` when `device<:GPUDevice` is used, otherwise a regular `Array`. Doesn't uses `adapt`, therefore always returns CuArray/Array
"""
DeviceArrayNotAdapt(::GPUDevice, x) = CuArray(x)
DeviceArrayNotAdapt(::CPUDevice, x) = Array(x)
DeviceArrayNotAdapt(dev::DeviceSetup, x) = DeviceArrayNotAdapt(dev.device, x)

"""
    launch_kernel!(device_setup::DeviceSetup, kernel!, ndrange, kernel_args...)

Launches the `kernel!` on the `device_setup` with `ndrange` computations over the kernel and arguments `kernel_args`. Returns an event.
"""
function launch_kernel!(device_setup::DeviceSetup, kernel!, ndrange, kernel_args...)
    device = device_setup.device_KA()
    n = device_setup.n 

    k! = kernel!(device, n)
    event = k!(kernel_args...; ndrange=ndrange)

    return event 
end