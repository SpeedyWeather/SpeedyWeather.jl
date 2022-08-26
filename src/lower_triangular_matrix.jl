"""
    L = LowerTriangularMatrix{T}(v::Vector{T},m::Int,n::Int)

A lower triangular matrix implementation that only stores the non-zero entries explicitly.
`L<:AbstractMatrix` although in general we have `L[i] != Matrix(L)[i]`, the former skips
zero entries, tha latter includes them."""
struct LowerTriangularMatrix{T} <: AbstractMatrix{T}
    v::Vector{T}    # non-zero elements unravelled into a vector
    m::Int          # number of rows
    n::Int          # number of columns

    LowerTriangularMatrix{T}(v,m,n) where T = length(v) == nonzeros(m,n) ?
        new(v,m,n) : error("$(length(v))-element Vector{$(eltype(v))} cannot be used to create a "*
            "$(m)x$(n) LowerTriangularMatrix{$T} with $(nonzeros(m,n)) non-zero entries.")
end

LowerTriangularMatrix(v::AbstractVector{T},m::Integer,n::Integer) where T = LowerTriangularMatrix{T}(v,m,n)

# SIZE ETC
Base.length(L::LowerTriangularMatrix) = length(L.v)     # define length as number of non-zero elements
Base.size(L::LowerTriangularMatrix) = (L.m,L.n)         # define size as matrix size
Base.sizeof(L::LowerTriangularMatrix) = sizeof(L.v)     # sizeof the underlying data vector

# CREATE INSTANCES (ZEROS, ONES, UNDEF)
for zeros_or_ones in (:zeros,:ones)
    @eval begin
        function Base.$zeros_or_ones(::Type{LowerTriangularMatrix{T}},m::Integer,n::Integer) where T
            return LowerTriangularMatrix($zeros_or_ones(T,nonzeros(m,n)),m,n)
        end
        
        # use Float64 as default if type T not provided
        Base.$zeros_or_ones(::Type{LowerTriangularMatrix},m::Integer,n::Integer) = $zeros_or_ones(LowerTriangularMatrix{Float64},m,n)
    end
end

Base.zero(L::LowerTriangularMatrix{T}) where T = zeros(LowerTriangularMatrix{T},size(L)...)
Base.one(L::LowerTriangularMatrix{T}) where T = ones(LowerTriangularMatrix{T},size(L)...)

function LowerTriangularMatrix{T}(::UndefInitializer,m::Integer,n::Integer) where T
    return LowerTriangularMatrix(Vector{T}(undef,nonzeros(m,n)),m,n)
end

Base.randn(::Type{LowerTriangularMatrix{T}},m::Integer,n::Integer) where T = LowerTriangularMatrix(randn(T,nonzeros(m,n)),m,n)
Base.rand(::Type{LowerTriangularMatrix{T}},m::Integer,n::Integer) where T = LowerTriangularMatrix(rand(T,nonzeros(m,n)),m,n)

Base.randn(::Type{LowerTriangularMatrix},m::Integer,n::Integer) = LowerTriangularMatrix(randn(nonzeros(m,n)),m,n)
Base.rand(::Type{LowerTriangularMatrix},m::Integer,n::Integer) = LowerTriangularMatrix(rand(nonzeros(m,n)),m,n)

# INDEXING
"""
    k = ij2k(   i::Integer,     # row index of matrix
                j::Integer,     # column index of matrix
                m::Integer)     # number of rows in matrix

Converts the index pair `i,j` of an `m`x`n` LowerTriangularMatrix `L` to a single
index `k` that indexes the same element in the corresponding vector that stores
only the lower triangle (the non-zero entries) of `L`."""
@inline ij2k(i::Integer, j::Integer,m::Integer) = i+(j-1)*m-j*(j-1)÷2
triangle_number(n::Integer) = n*(n+1)÷2
nonzeros(m::Integer,n::Integer) = m*n-triangle_number(n-1)

@inline function Base.getindex(L::LowerTriangularMatrix,k::Integer)
    @boundscheck 0 < k <= length(L.v) || throw(BoundsError(L,k))
    @inbounds r = L.v[k]
    return r
end

# for L[i,j] pull corresponding entry in data vector via index k or return zero
@inline function Base.getindex(L::LowerTriangularMatrix{T},i::Integer,j::Integer) where T
    @boundscheck (0 < i <= L.m && 0 < j <= L.n) || throw(BoundsError(L,(i,j)))
    j > i && return zero(T)
    k = ij2k(i,j,L.m)
    @inbounds r = L.v[k]
    return r
end

@inline Base.getindex(L::LowerTriangularMatrix,r::AbstractRange) = L.v[r]

@inline Base.setindex!(L::LowerTriangularMatrix,x,k::Integer) = setindex!(L.v,x,k)
@inline function Base.setindex!(L::LowerTriangularMatrix{T},x,i::Integer,j::Integer) where T
    j > i && return zero(T)
    k = ij2k(i,j,L.m)
    setindex!(L.v,x,k)
end

# @inline Base.setindex!(L::LowerTriangularMatrix,x::AbstractVector,r::AbstractRange) = setindex!(L.v,x,r)

"""
    unit_range = eachharmonic(L::LowerTriangular)

creates `unit_range::UnitRange` to loop over all non-zeros in a LowerTriangularMatrix `L`.
Like `eachindex` but skips the upper triangle with zeros in `L`."""
eachharmonic(L::LowerTriangularMatrix) = Base.OneTo(length(L))

"""
    unit_range = eachharmonic(Ls::LowerTriangularMatrix...)

creates `unit_range::UnitRange` to loop over all non-zeros in the LowerTriangularMatrices
provided as arguments. Checks bounds first. All LowerTriangularMatrixs need to be of the same size.
Like `eachindex` but skips the upper triangle with zeros in `L`."""
function eachharmonic(Ls::LowerTriangularMatrix...)
    for L in Ls
        @boundscheck size(L) == size(Ls[1]) || throw(BoundsError)
    end
    return eachharmonic(Ls[1])
end

# CONVERSIONS
""" 
    L = LowerTriangularMatrix(M)

Create a LowerTriangularMatrix `L` from Matrix `M` by copying over the non-zero elements in `M`."""
function LowerTriangularMatrix(M::AbstractMatrix{T}) where T
    m,n = size(M)
    L = LowerTriangularMatrix{T}(undef,m,n)
    k = 0
    @inbounds for j in 1:n      # only loop over lower triangle
        for i in j:m
            k += 1              # next element in lower triangle
            L[k] = M[i,j]       # and copy data into vector
        end
    end
    return L
end

Base.copy(L::LowerTriangularMatrix{T}) where T = LowerTriangularMatrix(L)

function LowerTriangularMatrix{T}(M::LowerTriangularMatrix) where T
    L = LowerTriangularMatrix{T}(undef,size(M)...)
    for i in eachindex(L,M)
        L[i] = convert(T,M[i])  # copy data over
    end
    return L
end

LowerTriangularMatrix(M::LowerTriangularMatrix{T}) where T = LowerTriangularMatrix{T}(M)

function Base.convert(::Type{LowerTriangularMatrix{T1}},L::LowerTriangularMatrix{T2}) where {T1,T2}
    return LowerTriangularMatrix{T1}(L.v,L.m,L.n)
end

"""
    fill!(L::LowerTriangularMatrix,x)

Fills the elements of `L` with `x`. Faster than fill!(::AbstractArray,x)
as only the non-zero elements in `L` are assigned with x."""
function Base.fill!(L::LowerTriangularMatrix{T}, x) where T
    xT = convert(T, x)
    @inbounds for i in eachindex(L.v)
        L[i] = xT
    end
end

# GPU methods
Adapt.adapt_structure(to, x::LowerTriangularMatrix{T}) where T = LowerTriangularMatrix(Adapt.adapt(to, x.v), x.m, x.n)
