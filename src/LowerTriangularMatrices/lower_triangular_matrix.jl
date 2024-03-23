import LinearAlgebra: tril!

"""
    L = LowerTriangularArray{T, N, ArrayType}(v::AbstractArray{T,N}, m::Int, n::Int)

A lower triangular matrix implementation that only stores the non-zero entries explicitly.
`L<:AbstractMatrix` although in general we have `L[i] != Matrix(L)[i]`, the former skips
zero entries, tha latter includes them."""
struct LowerTriangularArray{T, N, ArrayType <: AbstractArray{T}} <: AbstractArray{T,N}
    data::ArrayType     # non-zero elements unravelled into an array with one dimension less N-1
    m::Int              # number of rows
    n::Int              # number of columns

    LowerTriangularArray{T, N, ArrayType}(data, m, n) where {T, N, ArrayType<:AbstractArray} = length(data) == prod(size(data)[2:end]) * nonzeros(m, n) ? new(data, m, n) : error("$(size(data))-sized Array{$(eltype(data))} cannot be used to create a $(m)x$(n)x$(size(data)[2:end]) LowerTriangularArray{$T,$N,$ArrayType} with $(prod(size(data)[2:end]))x$(nonzeros(m, n)) non-zero entries.")  #TODO: adjust error string
end

const LowerTriangularMatrix{T, ArrayType} = LowerTriangularArray{T, 2, ArrayType}

LowerTriangularArray(data::ArrayType, m::Integer, n::Integer) where {T, N, ArrayType<:AbstractArray{T,N}} = LowerTriangularArray{T, N+1, ArrayType}(data, m, n)
LowerTriangularMatrix(data::ArrayType, m::Integer, n::Integer) where {T, ArrayType<:AbstractVector{T}} = LowerTriangularMatrix{T, ArrayType}(data, m, n)

# SIZE ETC
Base.length(L::LowerTriangularArray) = length(L.data)  # define length as number of non-zero elements
Base.size(L::LowerTriangularArray) = (L.m, L.n, size(L.data)[2:end]...)  # define size as matrix size (doesnt give an error, even for LowerTriangularMatrix due to tuple indexing logic)
Base.sizeof(L::LowerTriangularArray) = sizeof(L.data)  # sizeof the underlying data vector

# CREATE INSTANCES (ZEROS, ONES, UNDEF)
for zeros_or_ones in (:zeros, :ones)
    @eval begin
        function Base.$zeros_or_ones(::Type{LowerTriangularArray{T,N,ArrayType}}, m::Integer, n::Integer, I::Vararg{Integer,M}) where {T,N,M,ArrayType}
            return LowerTriangularArray($zeros_or_ones(T, nonzeros(m, n), I...), m, n)
        end
        
        # use Float64 and Vector as default if type T and ArrayType not provided
        Base.$zeros_or_ones(::Type{LowerTriangularMatrix}, m::Integer, n::Integer) = $zeros_or_ones(LowerTriangularMatrix{Float64, Vector{Float64}}, m, n)
    end
end

Base.zero(L::LowerTriangularArray{T,N,AT}) where {T,N,AT} = zeros(LowerTriangularArray{T,N,AT}, size(L)...)
Base.one(L::LowerTriangularArray{T,N,AT}) where {T,N,AT} = zeros(LowerTriangularArray{T,N,AT}, size(L)...)

function LowerTriangularArray{T,N,ArrayType}(::UndefInitializer, I::Vararg{Integer,N}) where {T,N,ArrayType<:AbstractArray{T}}
    return LowerTriangularArray(ArrayType(undef, nonzeros(I[1], I[2]), I[3:end]...), I[1], I[2])
end

function LowerTriangularMatrix{T,ArrayType}(::UndefInitializer, m::Integer, n::Integer) where {T,ArrayType<:AbstractArray{T}}
    return LowerTriangularMatrix(ArrayType(undef, nonzeros(m, n)), m, n)
end

Base.randn(::Type{LowerTriangularArray{T,N,ArrayType}}, m::Integer, n::Integer, I::Vararg{Integer,N}) where {T,N,ArrayType<:AbstractArray{T}} = LowerTriangularArray(ArrayType(randn(T, prod(I)*nonzeros(m, n), I...)), m, n)
Base.rand(::Type{LowerTriangularArray{T,N,ArrayType}}, m::Integer, n::Integer, I::Vararg{Integer,N}) where {T,N,ArrayType<:AbstractArray{T}} = LowerTriangularArray(ArrayType(rand(T, prod(I)*nonzeros(m, n), I...)), m, n)

Base.randn(::Type{LowerTriangularArray{T}}, m::Integer, n::Integer, I::Vararg{Integer,N}) where {T,N} = LowerTriangularArray(randn(T, prod(I)*nonzeros(m, n), I...), m, n)
Base.rand(::Type{LowerTriangularArray{T}}, m::Integer, n::Integer, I::Vararg{Integer,N}) where {T,N} = LowerTriangularArray(randn(T, prod(I)*nonzeros(m, n), I...), m, n)

Base.randn(::Type{LowerTriangularMatrix}, m::Integer, n::Integer) = LowerTriangularMatrix(randn(nonzeros(m, n)), m, n)
Base.rand(::Type{LowerTriangularMatrix}, m::Integer, n::Integer) = LowerTriangularMatrix(rand(nonzeros(m, n)), m, n)

# INDEXING
"""
    k = ij2k(   i::Integer,     # row index of matrix
                j::Integer,     # column index of matrix
                m::Integer)     # number of rows in matrix

Converts the index pair `i, j` of an `m`x`n` LowerTriangularMatrix `L` to a single
index `k` that indexes the same element in the corresponding vector that stores
only the lower triangle (the non-zero entries) of `L`."""
@inline ij2k(i::Integer, j::Integer, m::Integer) = i+(j-1)*m-j*(j-1)÷2
triangle_number(n::Integer) = n*(n+1)÷2
nonzeros(m::Integer, n::Integer) = m*n-triangle_number(n-1)

# direct indexing, no. indices have to be one less than `N` for the correct dimensionality, so N!=M
@inline function Base.getindex(L::LowerTriangularArray{T,N}, I::Vararg{Integer,M}) where {T,N,M}
    @boundscheck M == N-1 || throw(BoundsError(L, I))
    @boundscheck 0 < I[1] <= size(L.data,1) || throw(BoundsError(L, I[1])) 
    @inbounds r = getindex(L.data, I...) # TODO: Question: We could also just let the regular getindex handle the boundscheck
    return r
end

# integer + other indices (:, ranges, etc...)
@inline function Base.getindex(L::LowerTriangularArray{T,N}, k::Integer, I::Vararg{R,M}) where {T,N,R,M}
    @boundscheck M == N-2 || throw(BoundsError(L, I))
    @boundscheck 0 < k <= size(L.data,1) || throw(BoundsError(L, k)) 
    @inbounds r = getindex(L.data, k, I...) # TODO: Question: We could also just let the regular getindex handle the boundscheck
    return r
end

# l,m sph indexing, no. indices has to be equal to N 
@inline function Base.getindex(L::LowerTriangularArray{T,N}, I::Vararg{Integer,N}) where {T,N}
    i, j = I[1:2]
    @boundscheck (0 < i <= L.m && 0 < j <= L.n) || throw(BoundsError(L, (i, j)))
    j > i && return zero(getindex(L.data, 1, I[3:end]...)) # to get a zero element in the correct shape, we just take the zero element of some valid element, there are probably faster ways to do this, but I don't know how (espacially when ":" are used), and this is just a fallback anyway 
    k = ij2k(i, j, L.m)
    @inbounds r = getindex(L.data, k, I[3:end]...)
    return r
end

# l,m sph indexing with integer + other indices
@inline function Base.getindex(L::LowerTriangularArray{T,N}, i::Integer, j::Integer, I::Vararg{R,M}) where {T,N,R,M}
    @boundscheck M == N-2 || throw(BoundsError(L, I))
    @boundscheck (0 < i <= L.m && 0 < j <= L.n) || throw(BoundsError(L, (i, j)))
    j > i && return zero(getindex(L.data, 1, I...)) # to get a zero element in the correct shape, we just take the zero element of some valid element, there are probably faster ways to do this, but I don't know how, and this is just a fallback anyway 
    k = ij2k(i, j, L.m)
    @inbounds r = getindex(L.data, k, I...)
    return r
end

@inline Base.getindex(L::LowerTriangularArray, r::AbstractRange) = getindex(L.data,r)
@inline Base.getindex(L::LowerTriangularArray, r::AbstractRange, I...) = getindex(L.data, r, I...)

Base.@propagate_inbounds function Base.setindex!(L::LowerTriangularArray{T,N}, x, I::Vararg{Any, M}) where {T,N,M} 
    @boundscheck M == N-1 || throw(BoundsError(L, I))
    setindex!(L.data, x, I...)
end 

Base.@propagate_inbounds function Base.setindex!(L::LowerTriangularArray{T,N}, x, I::Vararg{Any, N}) where {T,N}
    i, j = I[1:2] # TODO: check if i,j are Int?
    @boundscheck i >= j || throw(BoundsError(L, (i, j)))
    k = ij2k(i, j, L.m)
    setindex!(L.data, x, k, I[3:end]...)
end

@inline Base.setindex!(L::LowerTriangularMatrix, x::AbstractVector, r::AbstractRange, I...) = setindex!(L.data, x, r, I...)

# propagate index to data vector
Base.eachindex(L ::LowerTriangularMatrix)    = eachindex(L.data)
Base.eachindex(Ls::LowerTriangularMatrix...) = eachindex((L.data for L in Ls)...)

"""
    unit_range = eachharmonic(L::LowerTriangularArray)

creates `unit_range::UnitRange` to loop over all non-zeros/spherical harmonics numbers in a LowerTriangularArray `L`.
Like `eachindex` but skips the upper triangle with zeros in `L`."""
eachharmonic(L::LowerTriangularArray) = axes(L.data, 1) # TODO: do want to keep this like that for n-dim arrays or change it? also adjust docstring when this is decided 

"""
    unit_range = eachharmonic(Ls::LowerTriangularMatrix...)

creates `unit_range::UnitRange` to loop over all non-zeros in the LowerTriangularMatrices
provided as arguments. Checks bounds first. All LowerTriangularMatrix's need to be of the same size.
Like `eachindex` but skips the upper triangle with zeros in `L`."""
eachharmonic(Ls::LowerTriangularArray...) = eachindex(Ls...)

# CONVERSIONS
""" 
    L = LowerTriangularMatrix(M)

Create a LowerTriangularMatrix `L` from Matrix `M` by copying over the non-zero elements in `M`."""
function LowerTriangularArray(M::Matrix{T}) where T # CPU version 
    m, n = size(M)
    L = LowerTriangularArray{T,2,Array{T,1}}(undef, size(M)...)
    
    k = 0
    @inbounds for j in 1:n      # only loop over lower triangle
        for i in j:m
            k += 1              # next element in lower triangle
            L[k] = M[i, j]       # and copy data into vector
        end
    end
    return L
end

function LowerTriangularArray(M::Array{T,3}) where T # CPU version 
    m, n = size(M)[1:2]
    L = LowerTriangularArray{T,3,Array{T,2}}(undef, size(M)...)
    
    k = 0
    @inbounds for j in 1:n      # only loop over lower triangle
        for i in j:m
            k += 1              # next element in lower triangle
            L[k, :] = M[i, j, :]       # and copy data into vector 
        end
    end
    return L
end

# helper function for conversion etc on GPU, returns indices of the lower triangle
lowertriangle_indices(M::AbstractMatrix) = tril!(trues(size(M)))
lowertriangle_indices(m::Integer, n::Integer) = tril!(trues((m,n)))

function lowertriangle_indices(M::AbstractArray{T,N}) where {T,N}
    @boundscheck N >= 3 || throw(BoundsError)

    indices = lowertriangle_indices(getindex(M,:,:,[1 for i=1:(N-2)]...)) # TODO: this assumes 1-indexing and won't work for some fancy indexed arrays
    indices = reshape(indices, size(indices,1), size(indices,2), [1 for i=1:(N-2)]...)
    
    repeat(indices, 1, 1, size(M)[3:end]...)
end  

# GPU version and fallback for higher dime  nsions
function LowerTriangularArray(M::ArrayType) where {T, N, ArrayType<:AbstractArray{T,N}} 
    m, n = size(M)[1:2]
    LowerTriangularArray(M[lowertriangle_indices(m,n),[Colon() for i=1:N-2]...], m, n)
end

LowerTriangularMatrix(M::AbstractMatrix) = LowerTriangularArray(M)

function Base.Matrix(L::LowerTriangularMatrix{T}) where T
    m, n = size(L)
    M = zeros(T, m, n)
    copyto!(M, L)
end
            
Base.copy(L::LowerTriangularArray{T}) where T = LowerTriangularArray(L)

function Base.copyto!(L1::LowerTriangularArray{T,N,ArrayType}, L2::LowerTriangularArray{S,N,ArrayType}) where {T, S, N, ArrayType}
    # if sizes don't match copy over the largest subset of indices
    size(L1) != size(L2) && return copyto!(L1, L2,  Base.OneTo(minimum(size.((L1, L2), 1))),
                                                    Base.OneTo(minimum(size.((L1, L2), 2))))

    L1.data .= convert.(T, L2.data)
    L1
end

function Base.copyto!(  L1::LowerTriangularArray{T,N,ArrayType},   # copy to L1
                        L2::LowerTriangularArray{S,N,ArrayType},      # copy from L2
                        ls::AbstractUnitRange,          # range of indices in 1st dim
                        ms::AbstractUnitRange) where {T, S, N, ArrayType}  # range of indices in 2nd dim

    lmax, mmax = size(L2)        # use the size of L2 for boundscheck
    @boundscheck maximum(ls) <= lmax || throw(BoundsError)
    @boundscheck maximum(ms) <= mmax || throw(BoundsError)

    lmax, mmax = size(L1)        # but the size of L1 to loop
    lm = 0
    @inbounds for m in 1:mmax
        for l in m:lmax
            lm += 1
            L1[lm] = (l in ls) && (m in ms) ? convert(T, L2[l, m]) : L1[lm]
        end
    end
    L1
end # do with intersection of the index arrays instead

function Base.copyto!(  L::LowerTriangularArray{T},    # copy to L
                        M::AbstractMatrix) where T      # copy from M

    @boundscheck size(L) == size(M) || throw(BoundsError)
    L.data .= convert.(T, M[lowertriangle_indices(M)])

    L
end

function Base.copyto!(  M::AbstractArray{T},               # copy to M
                        L::LowerTriangularArray) where T   # copy from L

    @boundscheck size(L) == size(M) || throw(BoundsError)

    lower_triangle_indices = lowertriangle_indices(M)
    upper_triangle_indices = @. ~lower_triangle_indices

    M[upper_triangle_indices] .= zero(T)
    M[lower_triangle_indices] = convert.(T, L.data)

    M
end

function LowerTriangularMatrix{T}(M::LowerTriangularMatrix{T2,ArrayType}) where {T,T2,ArrayType}
    L = LowerTriangularMatrix{T,ArrayType}(undef, size(M)...)
    copyto!(L, M)
    return L
end

function Base.convert(::Type{LowerTriangularArray{T1,N,ArrayType}}, L::LowerTriangularArray{T2,N,ArrayType}) where {T1, T2, N, ArrayType}
    return LowerTriangularArray{T1,N,ArrayType}(L.data, L.m, L.n)
end

function Base.similar(::LowerTriangularArray{T,N,ArrayType}, size::NTuple{N, Integer}) where {T, N, ArrayType}
    return LowerTriangularArray{T,N,ArrayType}(undef, size...)
end

function Base.similar(L::LowerTriangularArray{S,N,ArrayType}, ::Type{T}) where {T, S, N,ArrayType}
    return LowerTriangularArray{T,N,ArrayType}(undef, size(L)...)
end

Base.similar(L::LowerTriangularArray{T}) where T = similar(L, T)
 
# ARITHMETIC
# only mul/div with scalar and addition/subtraction, others are converted to Matrix
function Base.:(*)(L::LowerTriangularArray{T}, s::Number) where T
    sT = convert(T, s)
    LowerTriangularArray(L.data .* sT, L.m, L.n)
end

Base.:(*)(s::Number, L::LowerTriangularMatrix) = L*s         # commutative
Base.:(/)(L::LowerTriangularMatrix, s::Number) = L*inv(s)
Base.:(/)(s::Number, L::LowerTriangularMatrix) = L/s

# TODO: what if m / n unequal 
Base.:(+)(L1::LowerTriangularArray{T,N,ArrayType}, L2::LowerTriangularArray{T,N,ArrayType}) where {T,N,ArrayType} = LowerTriangularArray{T,N,ArrayType}(L1.data + L2.data, L1.m, L1.n)
function Base.:(+)(L1::LowerTriangularArray{T,N,ArrayType}, L2::LowerTriangularArray{S,N,ArrayType}) where {T,S,N,ArrayType}
    R = promote_type(T, S)
    LowerTriangularArray{R,N,ArrayType}(R.(L1.data) + R.(L2.data), L1.m, L1.n)
end

Base.:(-)(L1::LowerTriangularArray{T,N,ArrayType}, L2::LowerTriangularArray{T,N,ArrayType}) where {T,N,ArrayType} = LowerTriangularArray{T,N,ArrayType}(L1.data - L2.data, L1.m, L1.n)
function Base.:(-)(L1::LowerTriangularArray{T,N,ArrayType}, L2::LowerTriangularArray{S,N,ArrayType}) where {T,S,N,ArrayType}
    R = promote_type(T, S)
    LowerTriangularArray{R,N,ArrayType}(R.(L1.data) - R.(L2.data), L1.m, L1.n)
end

Base.:(-)(L::LowerTriangularArray) = LowerTriangularMatrix(-L.data, size(L)...)
Base.prod(L::LowerTriangularArray{NF}) where NF = zero(NF)

"""
    fill!(L::LowerTriangularArray, x)

Fills the elements of `L` with `x`. Faster than fill!(::AbstractArray, x)
as only the non-zero elements in `L` are assigned with x."""
function Base.fill!(L::LowerTriangularArray{T}, x) where T
    xT = convert(T, x)
    fill!(L.data, xT)

    L
end

function scale!(L::LowerTriangularArray{T}, s::Number) where T
    sT = convert(T, s)
    L.data .*= sT

    L
end

# Broadcast (more or less copied and adjusted from LinearAlgebra.jl)
import Base: similar, copyto!
import Base.Broadcast: BroadcastStyle, Broadcasted, DefaultArrayStyle
import LinearAlgebra: isstructurepreserving, fzeropreserving

struct LowerTriangularStyle <: Broadcast.AbstractArrayStyle{2} end

Base.BroadcastStyle(::Type{<:LowerTriangularMatrix}) = LowerTriangularStyle()

function Base.similar(bc::Broadcasted{LowerTriangularStyle}, ::Type{NF}) where NF
    inds = axes(bc)
    if isstructurepreserving(bc) || fzeropreserving(bc) 
        return LowerTriangularMatrix{NF}(undef, inds[1][end], inds[2][end])
    end
    return similar(convert(Broadcasted{DefaultArrayStyle{ndims(bc)}}, bc), NF)
end

LowerTriangularStyle(::Val{0}) = LowerTriangularStyle()
LowerTriangularStyle(::Val{1}) = LowerTriangularStyle()
LowerTriangularStyle(::Val{2}) = LowerTriangularStyle()

function Base.copyto!(dest::LowerTriangularMatrix{T}, bc::Broadcasted{<:LowerTriangularStyle}) where T
    axs = axes(dest)
    axes(bc) == axs || Broadcast.throwdm(axes(bc), axs)

    lmax, mmax = size(dest)
    lm = 0
    for m in 1:mmax
        for l in m:lmax
            lm += 1
            dest.data[lm] = Broadcast._broadcast_getindex(bc, CartesianIndex(l, m))
        end
    end
    return dest
end

# GPU methods (could be moved to an extension, in case this becomes a standalone package)
Adapt.adapt_structure(to, x::LowerTriangularArray) = LowerTriangularArray(Adapt.adapt(to, x.data), x.m, x.n)

