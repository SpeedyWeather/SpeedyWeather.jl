module Architectures

import KernelAbstractions

export AbstractArchitecture, AbstractCPU
export CPU, CPUStatic, GPU, ReactantDevice
export CUDAGPU, MetalGPU, ROCGPU
export array_type, on_architecture, architecture, device
export convert_to_device, ismatching, compatible_array_types, nonparametric_type
export synchronize
export memorysize

"""
    AbstractArchitecture

Abstract supertype for architectures supported by SpeedyWeather.
"""
abstract type AbstractArchitecture end

"""
    AbstractCPU <: AbstractArchitecture

Abstract supertype for CPU architectures supported by SpeedyWeather.
"""
abstract type AbstractCPU <: AbstractArchitecture end

"""
    CPU <: AbstractCPU

Run SpeedyWeather on one CPU node.
"""
struct CPU{D} <: AbstractCPU
    device::D
end

Base.show(io::IO, a::CPU) = print(io, "CPU")

CPU() = CPU(KernelAbstractions.CPU())
CPUStatic() = CPU(KernelAbstractions.CPU(; static = true))

"""
    GPU(device)

Return a GPU architecture using `device`.
`device` defauls to CUDA.CUDABackend(always_inline=true). 
CUDA has to be loaded for this to work. 
"""
struct GPU{D} <: AbstractArchitecture
    device::D
end

Base.show(io::IO, a::GPU) = print(io, "GPU($(device(a)))")

"""
    ReactantDevice(device)

Return a Reactant architecture to run SpeedyWeather using Reactant.
The actual device it's running on depends on the Reactant settings.
"""
struct ReactantDevice{D} <: AbstractArchitecture
    device::D
end

Base.show(io::IO, a::ReactantDevice) = print(io, "ReactantDevice")

# defined here so that it can be extended
function GPU end
function CUDAGPU end
function MetalGPU end
function ROCGPU end
function ReactantDevice end

#####
##### These methods are extended in SpeedyWeatherCUDAExt
#####

device(a::CPU) = a.device
device(a::GPU) = a.device

"""
    architecture(x)

Return the default architecture that's associated with the input `x`.
"""
architecture() = nothing
architecture(::Number) = nothing
architecture(::Array) = CPU()
architecture(::Type{<:Array}) = CPU()
architecture(a::SubArray) = architecture(parent(a))

"""
    architecture(a::AbstractArchitecture)

Return the architecture of `a`.
"""
architecture(a::AbstractArchitecture) = a

"""
    array_type(arch::AbstractArchitecture)

Return the array type that's used with the architecture `arch`.
"""
array_type(::CPU) = Array
array_type(::Type{<:CPU}) = Array

"""
    array_type(::AbstractArchitecture, NF::Type, N::Int)

Return the concrete array type that's used with the architecture 
for a number type `NF` and dimension `N`.
"""
array_type(::CPU, NF::Type, N::Int) = Array{NF, N}

"""
    compatible_array_types(arch::AbstractArchitecture)

Return the array types that are compatible with the architecture `arch`.
This includes the regular `array_type(arch)` as well as any other array types
that are compatible with the architecture, e.g. device arrays on GPU.
"""
compatible_array_types(arch::Type{<:AbstractArchitecture}) = (array_type(arch),) # fallback to array_type
compatible_array_types(arch::AbstractArchitecture) = (array_type(arch),) # fallback to array_type

ismatching(arch::Type{<:AbstractArchitecture}, array_T::Type{<:AbstractArray}) = any(arch_array -> array_T <: arch_array, compatible_array_types(arch))
ismatching(arch::AbstractArchitecture, array_T::Type{<:AbstractArray}) = ismatching(typeof(arch), array_T)
ismatching(arch::AbstractArchitecture, array::AbstractArray) = ismatching(arch, typeof(array))
ismatching(arch::AbstractArchitecture, a::SubArray) = ismatching(arch, parent(a))
ismatching(arch::AbstractArchitecture, a::Type{<:SubArray{T, N, A}}) where {T, N, A} = ismatching(arch, A)

# TODO: currently we just chech matching array types, sufficient?
ismatching(arch_1::AbstractArchitecture, arch_2::AbstractArchitecture) = array_type(arch_1) == array_type(arch_2)

"""
    nonparametric_type(array_T::Type{<:AbstractArray})

Strip away all the parameters of `array_T` and return the nonparametric type.
"""
nonparametric_type(::Type{<:Array}) = Array

# nonparametric_type for a SubArray is the arraytype it is viewing. Needed to construct new arrays from SubArrays!
nonparametric_type(::Type{<:SubArray{T, N, A}}) where {T, N, A} = nonparametric_type(A)
nonparametric_type(::Type{<:SubArray}) = SubArray   # if ArrayType A is not specified, return SubArray

# Tupled implementation
on_architecture(arch::AbstractArchitecture, t::Tuple) = map(elem -> on_architecture(arch, elem), t)   # map (not a generator) keeps the element types inferred
on_architecture(arch::AbstractArchitecture, nt::NamedTuple) = NamedTuple{keys(nt)}(on_architecture(arch, Tuple(nt)))

# pass-through for types that don't need architecture transfer
on_architecture(::AbstractArchitecture, x::Number) = x
on_architecture(::AbstractArchitecture, x::Base.RefValue) = x

# On architecture for array types
"""
    on_architecture(arch::AbstractArchitecture, a)

Return `a`, but on the architecture `arch`. 
"""
on_architecture(::CPU, a::Array) = a
on_architecture(::CPU, a::BitArray) = a
on_architecture(::CPU, a::SubArray{<:Any, <:Any, <:Array}) = a
on_architecture(::CPU, a::StepRangeLen) = a

# Convert arguments to GPU-compatible types
@inline convert_to_device(arch, args) = args
@inline convert_to_device(::CPU, args) = args

KernelAbstractions.synchronize(arch::AbstractArchitecture) = KernelAbstractions.synchronize(arch.device)
"""
    device_bytes(x) -> Int

Bytes of `x` held in device (e.g. GPU) memory if `x` is a device array, 0 otherwise.
Extended for `CuArray`, `ROCArray`, `MtlArray` and Reactant's concrete arrays in the
respective extensions."""
device_bytes(x) = 0

"""
    memorysize(x) -> Int

Like `Base.summarysize(x)` (memory reachable from `x` in bytes) but also counting the memory
of device arrays (e.g. `CuArray`), of which `Base.summarysize` only counts the host-side
wrapper. Every mutable object (incl. arrays) is counted only once."""
memorysize(x) = Base.summarysize(x) + _device_bytes(x, Base.IdSet{Any}())

# recursively sum `device_bytes` over all objects reachable from `x`
function _device_bytes(x, seen::Base.IdSet{Any})
    x isa Union{Module, Function, Type, Symbol, String} && return 0
    if ismutable(x)
        x in seen && return 0
        push!(seen, x)
    end
    nbytes = device_bytes(x)
    nbytes > 0 && return nbytes
    isbits(x) && return 0
    if x isa Array                  # host arrays: recurse into elements only if they're not bits
        isbitstype(eltype(x)) && return 0
        return sum((_device_bytes(x[i], seen) for i in eachindex(x) if isassigned(x, i)); init = 0)
    end
    # other objects incl. wrappers (views, LowerTriangularArray, Field, ...): recurse into fields
    return sum((_device_bytes(getfield(x, i), seen) for i in 1:nfields(x) if isdefined(x, i)); init = 0)
end

end
