struct LowerTriangularMatrix{T} <: AbstractMatrix{T}
    v::Vector{T}
    m::Int
    n::Int
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

function Base.getindex(L::LowerTriangularMatrix{T},i::Integer,j::Integer) where T
    j > i && return zero(T)
    @boundscheck (i > L.m || j > L.n) && throw(
        BoundsError("attempt to access $(L.m)x$(L.n) LowerTriangularMatrix{$T} at index [$i,$j]"))
    k = ij2k(i,j,L.m)
    return getindex(L.v,k)
end

function LowerTriangularMatrix(M::AbstractMatrix{T}) where T
    m,n = size(M)
    v = Vector{T}(undef,m*n-fibonacci(n))
    @inbounds for j in 1:n
        for i in 1:m
            v[ij2k(i,j,m)] = M[i,j]
        end
    end
    return LowerTriangularMatrix{T}(v,m,n)
end