"""A lower triangular matrix implementation that only stores the
non-zero entries explicitly."""
struct LowerTriangularMatrix{T} <: AbstractMatrix{T}
    v::Vector{T}    # non-zero elements unravelled into a vector
    m::Int          # number of rows
    n::Int          # number of columns
end

triangle_number(n::Integer) = n*(n-1)÷2
ij2k(i::Integer,j::Integer,m::Integer) = i+(j-1)*m-triangle_number(j)
Base.size(L::LowerTriangularMatrix) = (L.m,L.n)
Base.sizeof(L::LowerTriangularMatrix) = sizeof(L.v)

@inline function Base.getindex(L::LowerTriangularMatrix,k::Integer)
    @boundscheck 0 < k <= length(L.v) || throw(BoundsError(L,k))
    @inbounds r = L.v[k]
    return r
end

@inline function Base.getindex(L::LowerTriangularMatrix{T},i::Integer,j::Integer) where T
    @boundscheck (0 < i <= L.m && 0 < j <= L.n) || throw(BoundsError(L,(i,j)))
    j > i && return zero(T)
    k = ij2k(i,j,L.m)
    @inbounds r = L.v[k]
    return r
end

function LowerTriangularMatrix(M::AbstractMatrix{T}) where T
    m,n = size(M)
    v = Vector{T}(undef,m*n-triangle_number(n))
    @inbounds for j in 1:n
        for i in j:m
            v[ij2k(i,j,m)] = M[i,j]
        end
    end
    return LowerTriangularMatrix(v,m,n)
end