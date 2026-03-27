module LowerTriangularArraysReactantExt

using LowerTriangularArrays
using Reactant

const AnyReactantArray = Union{Reactant.AnyConcretePJRTArray, Reactant.AnyConcreteIFRTArray, Reactant.AnyTracedRArray}

Base.Broadcast.BroadcastStyle(::Type{LowerTriangularArray{T, N, ArrayType, S}}) where {T, N, ArrayType <: AnyReactantArray, S} = Base.Broadcast.BroadcastStyle(ArrayType)

end
