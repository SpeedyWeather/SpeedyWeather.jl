"""
    L = LowerTriangularMatrix{T}(v::Vector{T},m::Int,n::Int)

A lower triangular matrix implementation that only stores the non-zero entries explicitly.
`L<:AbstractMatrix` although in general we have `L[i] != Matrix(L)[i]`, the former skips
zero entries, tha latter includes them."""
struct LowerTriangularMatrix{T} <: AbstractMatrix{T}
    data::Vector{T}     # non-zero elements unravelled into a vector
    m::Int              # number of rows
    n::Int              # number of columns

    function LowerTriangularMatrix{T}(data, m, n) where {T}
        length(data) == nonzeros(m, n) ?
        new(data, m, n) :
        error("$(length(data))-element Vector{$(eltype(data))} cannot be used to create a " *
              "$(m)x$(n) LowerTriangularMatrix{$T} with $(nonzeros(m,n)) non-zero entries.")
    end
end

function LowerTriangularMatrix(data::AbstractVector{T}, m::Integer, n::Integer) where {T}
    LowerTriangularMatrix{T}(data, m, n)
end

# SIZE ETC
Base.length(L::LowerTriangularMatrix) = length(L.data)  # define length as number of non-zero elements
Base.size(L::LowerTriangularMatrix) = (L.m, L.n)         # define size as matrix size
Base.sizeof(L::LowerTriangularMatrix) = sizeof(L.data)  # sizeof the underlying data vector

# CREATE INSTANCES (ZEROS, ONES, UNDEF)
for zeros_or_ones in (:zeros, :ones)
    @eval begin
        function Base.$zeros_or_ones(::Type{LowerTriangularMatrix{T}}, m::Integer,
                                     n::Integer) where {T}
            return LowerTriangularMatrix($zeros_or_ones(T, nonzeros(m, n)), m, n)
        end

        # use Float64 as default if type T not provided
        function Base.$zeros_or_ones(::Type{LowerTriangularMatrix}, m::Integer, n::Integer)
            $zeros_or_ones(LowerTriangularMatrix{Float64}, m, n)
        end
    end
end

function Base.zero(L::LowerTriangularMatrix{T}) where {T}
    zeros(LowerTriangularMatrix{T}, size(L)...)
end
Base.one(L::LowerTriangularMatrix{T}) where {T} = ones(LowerTriangularMatrix{T}, size(L)...)

function LowerTriangularMatrix{T}(::UndefInitializer, m::Integer, n::Integer) where {T}
    return LowerTriangularMatrix(Vector{T}(undef, nonzeros(m, n)), m, n)
end

function Base.randn(::Type{LowerTriangularMatrix{T}}, m::Integer, n::Integer) where {T}
    LowerTriangularMatrix(randn(T, nonzeros(m, n)), m, n)
end
function Base.rand(::Type{LowerTriangularMatrix{T}}, m::Integer, n::Integer) where {T}
    LowerTriangularMatrix(rand(T, nonzeros(m, n)), m, n)
end

function Base.randn(::Type{LowerTriangularMatrix}, m::Integer, n::Integer)
    LowerTriangularMatrix(randn(nonzeros(m, n)), m, n)
end
function Base.rand(::Type{LowerTriangularMatrix}, m::Integer, n::Integer)
    LowerTriangularMatrix(rand(nonzeros(m, n)), m, n)
end

# INDEXING
"""
    k = ij2k(   i::Integer,     # row index of matrix
                j::Integer,     # column index of matrix
                m::Integer)     # number of rows in matrix

Converts the index pair `i,j` of an `m`x`n` LowerTriangularMatrix `L` to a single
index `k` that indexes the same element in the corresponding vector that stores
only the lower triangle (the non-zero entries) of `L`."""
@inline ij2k(i::Integer, j::Integer, m::Integer) = i + (j - 1) * m - j * (j - 1) ÷ 2
triangle_number(n::Integer) = n * (n + 1) ÷ 2
nonzeros(m::Integer, n::Integer) = m * n - triangle_number(n - 1)

@inline function Base.getindex(L::LowerTriangularMatrix, k::Integer)
    @boundscheck 0 < k <= length(L.data) || throw(BoundsError(L, k))
    @inbounds r = L.data[k]
    return r
end

# for L[i,j] pull corresponding entry in data vector via index k or return zero
@inline function Base.getindex(L::LowerTriangularMatrix{T}, i::Integer,
                               j::Integer) where {T}
    @boundscheck (0 < i <= L.m && 0 < j <= L.n) || throw(BoundsError(L, (i, j)))
    j > i && return zero(T)
    k = ij2k(i, j, L.m)
    @inbounds r = L.data[k]
    return r
end

@inline Base.getindex(L::LowerTriangularMatrix, r::AbstractRange) = L.data[r]

Base.@propagate_inbounds function Base.setindex!(L::LowerTriangularMatrix, x, k::Integer)
    setindex!(L.data, x, k)
end
Base.@propagate_inbounds function Base.setindex!(L::LowerTriangularMatrix{T}, x, i::Integer,
                                                 j::Integer) where {T}
    @boundscheck i >= j || throw(BoundsError(L, (i, j)))
    k = ij2k(i, j, L.m)
    setindex!(L.data, x, k)
end

@inline function Base.setindex!(L::LowerTriangularMatrix, x::AbstractVector,
                                r::AbstractRange)
    setindex!(L.data, x, r)
end

# propagate index to data vector
Base.eachindex(L::LowerTriangularMatrix) = eachindex(L.data)
Base.eachindex(Ls::LowerTriangularMatrix...) = eachindex((L.data for L in Ls)...)

"""
    unit_range = eachharmonic(L::LowerTriangular)

creates `unit_range::UnitRange` to loop over all non-zeros in a LowerTriangularMatrix `L`.
Like `eachindex` but skips the upper triangle with zeros in `L`."""
eachharmonic(L::LowerTriangularMatrix) = eachindex(L.data)

"""
    unit_range = eachharmonic(Ls::LowerTriangularMatrix...)

creates `unit_range::UnitRange` to loop over all non-zeros in the LowerTriangularMatrices
provided as arguments. Checks bounds first. All LowerTriangularMatrix's need to be of the same size.
Like `eachindex` but skips the upper triangle with zeros in `L`."""
eachharmonic(Ls::LowerTriangularMatrix...) = eachindex(Ls...)

# CONVERSIONS
""" 
    L = LowerTriangularMatrix(M)

Create a LowerTriangularMatrix `L` from Matrix `M` by copying over the non-zero elements in `M`."""
function LowerTriangularMatrix(M::AbstractMatrix{T}) where {T}
    m, n = size(M)
    L = LowerTriangularMatrix{T}(undef, m, n)
    k = 0
    @inbounds for j in 1:n      # only loop over lower triangle
        for i in j:m
            k += 1              # next element in lower triangle
            L[k] = M[i, j]       # and copy data into vector
        end
    end
    return L
end

function Base.Matrix(L::LowerTriangularMatrix{T}) where {T}
    m, n = size(L)
    M = zeros(T, m, n)
    k = 0
    @inbounds for j in 1:n
        for i in j:m
            k += 1
            M[i, j] = L[k]
        end
    end
    return M
end

Base.copy(L::LowerTriangularMatrix{T}) where {T} = LowerTriangularMatrix(L)

function Base.copyto!(L1::LowerTriangularMatrix{T}, L2::LowerTriangularMatrix) where {T}
    # if sizes don't match copy over the largest subset of indices
    size(L1) != size(L2) && return copyto!(L1, L2, Base.OneTo(minimum(size.((L1, L2), 1))),
                   Base.OneTo(minimum(size.((L1, L2), 2))))

    @inbounds for i in eachindex(L1, L2)
        L1[i] = convert(T, L2[i])
    end
    L1
end

function Base.copyto!(L1::LowerTriangularMatrix{T},   # copy to L1
                      L2::LowerTriangularMatrix,      # copy from L2
                      ls::AbstractUnitRange,          # range of indices in 1st dim
                      ms::AbstractUnitRange) where {T}  # range of indices in 2nd dim
    lmax, mmax = size(L2)        # use the size of L2 for boundscheck
    @boundscheck maximum(ls) <= lmax || throw(BoundsError)
    @boundscheck maximum(ms) <= mmax || throw(BoundsError)

    lmax, mmax = size(L1)        # but the size of L1 to loop
    lm = 0
    @inbounds for m in 1:maximum(ms)
        for l in m:lmax
            lm += 1
            L1[lm] = (l in ls) && (m in ms) ? convert(T, L2[l, m]) : L1[lm]
        end
    end
    L1
end

function Base.copyto!(L::LowerTriangularMatrix{T},    # copy to L
                      M::AbstractMatrix) where {T}      # copy from M
    @boundscheck size(L) == size(M) || throw(BoundsError)
    lmax, mmax = size(L)

    lm = 0
    @inbounds for m in 1:mmax
        for l in m:lmax
            lm += 1
            L[lm] = convert(T, M[l, m])
        end
    end
    L
end

function Base.copyto!(M::AbstractMatrix{T},               # copy to M
                      L::LowerTriangularMatrix) where {T}   # copy from L
    @boundscheck size(L) == size(M) || throw(BoundsError)
    lmax, mmax = size(L)

    lm = 0
    @inbounds for m in 1:mmax
        for l in 1:(m - 1)          # zero for upper triangle (excl diagonal)
            M[l, m] = zero(T)
        end

        for l in m:lmax         # convert and copy for lower triangle
            lm += 1
            M[l, m] = convert(T, M[lm])
        end
    end
    M
end

function LowerTriangularMatrix{T}(M::LowerTriangularMatrix) where {T}
    L = LowerTriangularMatrix{T}(undef, size(M)...)
    copyto!(L, M)
    return L
end

LowerTriangularMatrix(M::LowerTriangularMatrix{T}) where {T} = LowerTriangularMatrix{T}(M)

function Base.convert(::Type{LowerTriangularMatrix{T1}},
                      L::LowerTriangularMatrix{T2}) where {T1, T2}
    return LowerTriangularMatrix{T1}(L.data, L.m, L.n)
end

function Base.similar(::LowerTriangularMatrix{T}, size::Integer...) where {T}
    return LowerTriangularMatrix{T}(undef, size...)
end

function Base.similar(::LowerTriangularMatrix{T}, size::NTuple{2, Integer}) where {T}
    return LowerTriangularMatrix{T}(undef, size...)
end

function Base.similar(L::LowerTriangularMatrix, ::Type{T}) where {T}
    return LowerTriangularMatrix{T}(undef, size(L)...)
end

Base.similar(L::LowerTriangularMatrix{T}) where {T} = similar(L, T)

# ARITHMETIC
# only mul/div with scalar and addition/subtraction, others are converted to Matrix
function Base.:(*)(L::LowerTriangularMatrix{T}, s::Number) where {T}
    M = similar(L)
    sT = convert(T, s)
    @inbounds for i in eachindex(L, M)
        M[i] = L[i] * sT
    end
    M
end

Base.:(*)(s::Number, L::LowerTriangularMatrix) = L * s         # commutative
Base.:(/)(L::LowerTriangularMatrix, s::Number) = L * inv(s)
Base.:(/)(s::Number, L::LowerTriangularMatrix) = L / s

function Base.:(+)(L1::LowerTriangularMatrix, L2::LowerTriangularMatrix)
    T = promote_type(eltype(L1), eltype(L2))
    M = similar(L1, T)
    @inbounds for i in eachindex(M, L1, L2)
        M[i] = L1[i] + L2[i]
    end
    M
end

function Base.:(-)(L1::LowerTriangularMatrix, L2::LowerTriangularMatrix)
    T = promote_type(eltype(L1), eltype(L2))
    M = similar(L1, T)
    @inbounds for i in eachindex(M, L1, L2)
        M[i] = L1[i] - L2[i]
    end
    M
end

Base.:(-)(L::LowerTriangularMatrix) = LowerTriangularMatrix(-L.data, size(L)...)
Base.prod(L::LowerTriangularMatrix{NF}) where {NF} = zero(NF)

"""
    fill!(L::LowerTriangularMatrix,x)

Fills the elements of `L` with `x`. Faster than fill!(::AbstractArray,x)
as only the non-zero elements in `L` are assigned with x."""
function Base.fill!(L::LowerTriangularMatrix{T}, x) where {T}
    xT = convert(T, x)
    @inbounds for i in eachindex(L)
        L[i] = xT
    end
end

function scale!(L::LowerTriangularMatrix{T}, s::Number) where {T}
    sT = convert(T, s)
    @inbounds for i in eachindex(L)
        L[i] *= sT
    end
    L
end

# Broadcast (more or less copied and adjusted from LinearAlgebra.jl)
import Base: similar, copyto!
import Base.Broadcast: BroadcastStyle, Broadcasted, DefaultArrayStyle
import LinearAlgebra: isstructurepreserving, fzeropreserving

struct LowerTriangularStyle <: Broadcast.AbstractArrayStyle{2} end

Base.BroadcastStyle(::Type{<:LowerTriangularMatrix}) = LowerTriangularStyle()

function Base.similar(bc::Broadcasted{LowerTriangularStyle}, ::Type{NF}) where {NF}
    inds = axes(bc)
    if isstructurepreserving(bc) || fzeropreserving(bc)
        return LowerTriangularMatrix{NF}(undef, inds[1][end], inds[2][end])
    end
    return similar(convert(Broadcasted{DefaultArrayStyle{ndims(bc)}}, bc), NF)
end

LowerTriangularStyle(::Val{0}) = LowerTriangularStyle()
LowerTriangularStyle(::Val{1}) = LowerTriangularStyle()
LowerTriangularStyle(::Val{2}) = LowerTriangularStyle()

function Base.copyto!(dest::LowerTriangularMatrix{T},
                      bc::Broadcasted{<:LowerTriangularStyle}) where {T}
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

# GPU methods
function Adapt.adapt_structure(to, x::LowerTriangularMatrix{T}) where {T}
    LowerTriangularMatrix(Adapt.adapt(to, x.data), x.m, x.n)
end
