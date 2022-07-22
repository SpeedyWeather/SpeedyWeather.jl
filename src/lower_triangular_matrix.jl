"""A lower triangular matrix implementation that only stores the
non-zero entries explicitly."""
struct LowerTriangularMatrix{T} <: AbstractMatrix{T}
    v::Vector{T}    # non-zero elements unravelled into a vector
    m::Int          # number of rows
    n::Int          # number of columns
end

function fibonacci(j::Integer)
    ∑ = zero(j)
    for jj in 1:j-1
        ∑ += jj
    end
    return ∑
end

ij2k(i::Integer,j::Integer,m::Integer) = i+(j-1)*m-fibonacci(j)
Base.size(L::LowerTriangularMatrix) = (L.m,L.n)
Base.sizeof(L::LowerTriangularMatrix) = sizeof(L.v)

Base.getindex(L::LowerTriangularMatrix,k::Integer) = getindex(@inbounds L.v[k])

function Base.getindex(L::LowerTriangularMatrix{T},i::Integer,j::Integer) where T
    @boundscheck (i > L.m || j > L.n) && throw(BoundsError(L,(i,j)))
    j > i && return zero(T)
    return getindex(L.v,ij2k(i,j,L.m))
end

function LowerTriangularMatrix(M::AbstractMatrix{T}) where T
    m,n = size(M)
    v = Vector{T}(undef,m*n-fibonacci(n))
    @inbounds for j in 1:n
        for i in j:m
            v[ij2k(i,j,m)] = M[i,j]
        end
    end
    return LowerTriangularMatrix{T}(v,m,n)
end