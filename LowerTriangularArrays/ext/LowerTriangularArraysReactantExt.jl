module LowerTriangularArraysReactantExt

using LowerTriangularArrays
using Reactant

const AnyReactantArray = Union{Reactant.AnyConcretePJRTArray, Reactant.AnyConcreteIFRTArray, Reactant.AnyTracedRArray}

Base.Broadcast.BroadcastStyle(::Type{LowerTriangularArray{T, N, ArrayType, S}}) where {T, N, ArrayType <: AnyReactantArray, S} = Base.Broadcast.BroadcastStyle(ArrayType)

# Unwrap LowerTriangularArrays with Reactant-backed data to their underlying `.data` so
# that Reactant's broadcast JIT compiler sees plain ConcretePJRTArrays (avoids scalar
# `lta[i]` fallback via LowerTriangularArray's AbstractArray interface).
Base.Broadcast.broadcastable(lta::LowerTriangularArray{T, N, ArrayType}) where {T, N, ArrayType <: AnyReactantArray} = lta.data

# Forward `copyto!(lta, bc)` to `copyto!(lta.data, bc)` so Reactant's broadcast copyto!
# can write directly into the ConcretePJRTArray.
for AT in (:ConcretePJRTArray, :ConcreteIFRTArray)
    @eval function Base.copyto!(
            dest::LowerTriangularArray{T, N, ArrayType},
            bc::Base.Broadcast.Broadcasted{Base.Broadcast.ArrayStyle{Reactant.$AT}},
        ) where {T, N, ArrayType <: AnyReactantArray}
        Base.copyto!(dest.data, bc)
        return dest
    end
end

end
