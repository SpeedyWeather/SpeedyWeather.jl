"""
    abstract type AbstractDevice

Supertype of all devices SpeedyWeather.jl can run on.
"""
abstract type AbstractDevice end

export CPU, GPU

"""
    CPU <: AbstractDevice

Indicates that SpeedyWeather.jl runs on a single CPU.
"""
struct CPU <: AbstractDevice end

"""
    GPU <: AbstractDevice

Indicates that SpeedyWeather.jl runs on a single GPU.
"""
struct GPU <: AbstractDevice end


"""$(TYPEDSIGNATURES)
Default array type on `device`."""
default_array_type(device::AbstractDevice) = default_array_type(typeof(device))
default_array_type(::Type{CPU}) = Array

"""$(TYPEDSIGNATURES)
Return used device for use with KernelAbstractions.
"""
Device_KernelAbstractions(::CPU) = KernelAbstractions.CPU

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
Adapts `x` to an `Array` when `device::CPU` is used. Define for `CPU` for compatibility with adapt to CuArrays etc.
Uses `adapt`, thus also can return SubArrays etc."""
DeviceArray(::CPU, x) = Adapt.adapt(Array, x)
DeviceArray(dev::DeviceSetup, x) = DeviceArray(dev.device, x)

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