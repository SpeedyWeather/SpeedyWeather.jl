module RingGridsReactantExt

using RingGrids
using Reactant

const AnyReactantArray = Union{Reactant.AnyConcretePJRTArray, Reactant.AnyConcreteIFRTArray, Reactant.AnyTracedRArray}

Base.Broadcast.BroadcastStyle(::Type{F}) where {F <: AbstractField{T, N, ArrayType, Grid}} where {T, N, ArrayType <: AnyReactantArray, Grid} = Base.Broadcast.BroadcastStyle(ArrayType)

end
