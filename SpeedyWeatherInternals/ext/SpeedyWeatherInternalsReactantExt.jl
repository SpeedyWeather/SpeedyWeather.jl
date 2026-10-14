module SpeedyWeatherInternalsReactantExt

using KernelAbstractions
using SpeedyWeatherInternals
using Reactant

import SpeedyWeatherInternals.Architectures: Architectures, CPU, ReactantDevice, array_type, architecture, on_architecture, compatible_array_types, nonparametric_type, device
import SpeedyWeatherInternals.Utils: _jit, _set_scalar!

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

Architectures.compatible_array_types(::ReactantDevice) = (ConcreteRArray, AnyConcreteReactantArray, Reactant.AnyTracedRArray)
Architectures.compatible_array_types(::Type{<:ReactantDevice}) = (ConcreteRArray, AnyConcreteReactantArray, Reactant.AnyTracedRArray)

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

# TODO: move upstream to Reactant
# Broadcasting into views of concrete Reactant arrays (e.g. a variable that is a view into a fused
# parent array) has no Reactant broadcast style, falls back to Base's elementwise scalar indexing
# loop and silently produces wrong results. Use Reactant's broadcast style and compile instead.
_broadcast_into!(dest, f, args...) = (dest .= f.(args...); nothing)

for ArrayType in (:ConcretePJRTArray, :ConcreteIFRTArray)
    @eval begin
        Base.Broadcast.BroadcastStyle(::Type{<:SubArray{<:Any, <:Any, <:$ArrayType}}) =
            Base.Broadcast.ArrayStyle{$ArrayType}()

        function Base.copyto!(
                dest::SubArray{<:Any, <:Any, <:$ArrayType},
                bc::Base.Broadcast.Broadcasted{Base.Broadcast.ArrayStyle{$ArrayType}},
            )
            bc = Base.Broadcast.flatten(bc)
            Reactant.@jit _broadcast_into!(dest, bc.f, bc.args...)
            return dest
        end
    end
end

# TODO: move upstream to Reactant
# `deepcopy` of a concrete Reactant array does not register the copy in deepcopy's `IdDict`, so every
# reference to the same array gets its own independent copy and aliasing is lost, e.g. between a
# fused parent array and the variables that are views into it. Register the copy as `deepcopy` of an
# `Array` (or Reactant's `TracedRArray`) does.
for ArrayType in (:ConcretePJRTArray, :ConcreteIFRTArray)
    @eval function Base.deepcopy_internal(x::$ArrayType, stackdict::IdDict)
        haskey(stackdict, x) && return stackdict[x]::typeof(x)
        y = copy(x)
        stackdict[x] = y
        return y
    end
end

# For @maybe_jit macro - extend _jit for ReactantDevice.
# If we are already inside a Reactant compile/trace context, call f directly (no nested @jit).
# Otherwise JIT-compile on the fly, raising KernelAbstractions kernels to StableHLO as in the compiled
# time step (see `time_stepping!` in SpeedyWeatherReactantExt): kernels that aren't raised can be handed
# operands with a non-default layout on CUDA and read scrambled data (EnzymeAD/Reactant.jl#3377)
function _jit(::ReactantDevice, f, args...; kwargs...)
    if Reactant.ReactantCore.within_compile()
        return f(args...; kwargs...)
    else
        return Reactant.@jit raise = true f(args...; kwargs...)
    end
end

# In Reactant actually use @allowscalar (in contrast to regular model which needs a workaround for Enzyme)
function _set_scalar!(::ReactantDevice, x::AbstractArray, i::Integer, v::Number)
    Reactant.@allowscalar x[i] = v
    return x
end

end
