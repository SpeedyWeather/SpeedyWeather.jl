"""
    AssociatedLegendrePolArray{T,N,M,V} <: AbstractArray{T,N}

Type that wraps around a `LowerTriangularArray{T,M,V}` but is a subtype of `AbstractArray{T,M+1}`. 
This enables easier use with AssociatedLegendrePolynomials.jl which otherwise couldn't use the 
"matrix-style" (l, m) indexing of `LowerTriangularArray`. This type however doesn't support any
other operations than indexing and is purerly intended for internal purposes. 
"""
struct AssociatedLegendrePolArray{T,N,M,V} <: AbstractArray{T,N}
    data::LowerTriangularArray{T,M,V}
end 

"""2-dimensional `AssociatedLegendrePolArray` of type `T`` with its non-zero entries unravelled into a `Vector{T}`"""
const AssociatedLegendrePolMatrix{T} = AssociatedLegendrePolArray{T, 2, 1, Vector{T}}

Base.size(A::AssociatedLegendrePolArray) = size(A.data; as=Matrix)
@inline Base.getindex(A::AssociatedLegendrePolArray, I...) = getindex(A.data, I...)
@inline Base.setindex!(A::AssociatedLegendrePolArray, v, I...) = setindex!(A.data, v, I...)

# AssociatedLegendrePolynomials.jl does indexing with an empty cartesian index that would otherwise lead to an error
@inline Base.getindex(A::AssociatedLegendrePolArray, i::CartesianIndex{0}, I...) = getindex(A.data, I...)
@inline Base.setindex!(A::AssociatedLegendrePolArray, v, i::CartesianIndex{0}, I...) = setindex!(A.data, v, I...)
