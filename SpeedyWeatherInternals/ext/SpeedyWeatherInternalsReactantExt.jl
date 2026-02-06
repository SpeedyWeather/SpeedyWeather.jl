module SpeedyWeatherInternalsReactantExt

using KernelAbstractions
using SpeedyWeatherInternals
using Reactant

import SpeedyWeatherInternals.Architectures: Architectures, CPU, ReactantDevice, array_type, architecture, on_architecture, compatible_array_types, nonparametric_type, device
import SpeedyWeatherInternals.Utils: _jit

# grab the proper KernelAbstractions backend for Reactant
const ReactantKernelAbstractionsExt = Base.get_extension(
    Reactant, :ReactantKernelAbstractionsExt
)
const ReactantBackend = ReactantKernelAbstractionsExt.ReactantBackend
const AnyConcreteReactantArray = Union{Reactant.AnyConcretePJRTArray, Reactant.AnyConcreteIFRTArray}

Architectures.ReactantDevice() = ReactantDevice(ReactantBackend())

Architectures.array_type(::ReactantDevice) = ConcreteRArray
Architectures.array_type(::Type{<:ReactantDevice}) = ConcreteRArray
Architectures.array_type(::ReactantDevice, NF::Type, N::Int) = ConcretePJRTArray{NF, N, 1}

Architectures.compatible_array_types(::ReactantDevice) = (ConcreteRArray, AnyConcreteReactantArray)
Architectures.compatible_array_types(::Type{<:ReactantDevice}) = (ConcreteRArray, AnyConcreteReactantArray)

Architectures.nonparametric_type(::Type{<:ConcreteRArray}) = ConcreteRArray

device(dev::ReactantDevice) = dev.device

architecture(::AnyConcreteReactantArray) = ReactantDevice()
architecture(::Reactant.AnyTracedRArray) = ReactantDevice()
architecture(::Type{<:AnyConcreteReactantArray}) = ReactantDevice()
architecture(::Type{<:Reactant.AnyTracedRArray}) = ReactantDevice()

on_architecture(::ReactantDevice, a::Reactant.AnyTracedRArray) = a
on_architecture(::CPU, a::AnyConcreteReactantArray) = Array(a)
on_architecture(::CPU, a::SubArray{<:Any, <:Any, <:AnyConcreteReactantArray}) = Array(a)

# which arrays should converted?
const ArraysToRArray = Union{
    Array,
    Reactant.AnyConcretePJRTArray,
    Reactant.AnyConcreteIFRTArray,
    BitArray,
    SubArray{<:Any, <:Any, <:Array},
}

on_architecture(::ReactantDevice, a::ArraysToRArray) = Reactant.to_rarray(a)

# defined for automatic conversion during struct construction
# could also generalize those, but currently we are very conservative here to only define exactly those we need
Base.convert(::Type{ConcretePJRTArray{T, 1, 1}}, a::AbstractVector{T}) where {T <: Number} = Reactant.to_rarray(a)
Base.convert(::Type{ConcretePJRTArray{T, 1, 1}}, a::AbstractVector{S}) where {T <: Number, S <: Number} = Reactant.to_rarray(T.(a))
Base.convert(::Type{ConcretePJRTArray{T, 2}}, a::AbstractMatrix{T}) where {T <: Number} = Reactant.to_rarray(a)
Base.convert(::Type{ConcretePJRTArray{T, 2}}, a::AbstractMatrix{S}) where {T <: Number, S <: Number} = Reactant.to_rarray(T.(a))
Base.convert(::Type{ConcretePJRTArray{T, 2, 1}}, a::AbstractMatrix{T}) where {T <: Number} = Reactant.to_rarray(a)
Base.convert(::Type{ConcretePJRTArray{T, 2, 1}}, a::AbstractMatrix{S}) where {T <: Number, S <: Number} = Reactant.to_rarray(T.(a))
Base.convert(::Type{ConcretePJRTArray{T, 3, 1}}, a::AbstractArray{T, 3}) where {T <: Number} = Reactant.to_rarray(a)
Base.convert(::Type{ConcretePJRTArray{T, 3, 1}}, a::AbstractArray{S, 3}) where {T <: Number, S <: Number} = Reactant.to_rarray(T.(a))

Reactant.ConcretePJRTArray{T, N, D}(a::AbstractArray{T, N}) where {T, N, D} = Reactant.to_rarray(a)

# For @maybe_jit macro - extend _jit for ReactantDevice
_jit(::ReactantDevice, f, args...; kwargs...) = Reactant.@jit f(args...; kwargs...)

end
