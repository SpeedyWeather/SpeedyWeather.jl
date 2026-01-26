module SpeedyWeatherInternalsReactantExt 

using Reactant: AnyConcreteRArray

import SpeedyWeatherInternals.Architectures: Architectures, ReactantDevice, array_type, architecture, on_architecture, compatible_array_types, nonparametric_type

# grab the proper KernelAbstractions backend for Reactant
const ReactantKernelAbstractionsExt = Base.get_extension(
    Reactant, :ReactantKernelAbstractionsExt
)
const ReactantBackend = ReactantKernelAbstractionsExt.ReactantBackend
const AnyConcreteReactantArray = Union{Reactant.AnyConcretePJRTArray,Reactant.AnyConcreteIFRTArray}

Architectures.ReactantDevice() = ReactantDevice(ReactantBackend())

Architectures.array_type(::ReactantDevice) = ConcreteRArray
Architectures.array_type(::Type{<:ReactantDevice}) = ConcreteRArray
Architectures.array_type(::ReactantDevice, NF::Type, N::Int) = ConcreteRArray{NF, N, Runtime.Mem.HIPBuffer}

Architectures.compatible_array_types(::ReactantDevice) = (ConcreteRArray, AnyConcreteReactantArray)
Architectures.compatible_array_types(::Type{<:ReactantDevice}) = (ConcreteRArray, AnyConcreteReactantArray)

Architectures.nonparametric_type(::Type{<:ConcreteRArray}) = ConcreteRArray

device(::ReactantDevice) = ReactantBackend()

architecture(::AnyConcreteReactantArray) = ReactantDevice()
architecture(::Reactant.AnyTracedRArray) = ReactantDevice()

on_architecture(::ReactantDevice, a::Reactant.AnyTracedRArray) = a
on_architecture(::CPU, a::AnyConcreteReactantArray) = Array(a)
on_architecture(::CPU, a::SubArray{<:Any,<:Any,<:AnyConcreteReactantArray}) = Array(a)

# which arrays should converted? 
const ArraysToRArray = Union{Array,
    Reactant.AnyConcretePJRTArray,
    # Reactant.AnyConcreteIFRTArray, # needed?
    BitArray,
    SubArray{<:Any,<:Any,<:Array}}

on_architecture(::ReactantDevice, a::ArraysToRArray) = Reactant.to_rarray(a)

end