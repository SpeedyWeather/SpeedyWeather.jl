"""
$(TYPEDSIGNATURES)
A lower triangular array implementation that only stores the non-zero entries explicitly.
`L<:AbstractArray{T,N-1}` although we do allow both "flat" `N-1`-dimensional indexing and 
additional `N`-dimensional or "matrix-style" indexing.

Supports n-dimensional lower triangular arrays, so that for all trailing dimensions `L[:, :, ..]`
is a matrix in lower triangular form, e.g. a (5x5x3)-LowerTriangularArray would hold 3 lower
triangular matrices."""
struct LowerTriangularArray{T, N, ArrayType <: AbstractArray{T,N}} <: AbstractArray{T,N}
    data::ArrayType     # non-zero elements unravelled into an in which the lower triangle is flattened
    m::Int              # number of rows
    n::Int              # number of columns

    LowerTriangularArray{T, N, ArrayType}(data, m, n) where {T, N, ArrayType<:AbstractArray{T}} =
        check_lta_input_array(data, m, n, N) ? 
        new(data, m, n)  :
        error(lta_error_message(data, m, n, T, N, ArrayType))
end

check_lta_input_array(data, m, n, N) =
    (ndims(data) == N) & (length(data) == prod(size(data)[2:end]) * nonzeros(m, n)) 

function lta_error_message(data, m, n, T, N, ArrayType) 
    return "$(Base.dims2string(size(data))) $(typeof(data)) cannot be used to create "*
            "a $(Base.dims2string(matrix_size(data, m, n))) LowerTriangularArray{$T, $N, $ArrayType}"
end 

"""2-dimensional `LowerTriangularArray` of type `T`` with its non-zero entries unravelled into a `Vector{T}`"""
const LowerTriangularMatrix{T} = LowerTriangularArray{T, 1, Vector{T}}

LowerTriangularArray(data::ArrayType, m::Integer, n::Integer) where {T, N, ArrayType<:AbstractArray{T,N}} = LowerTriangularArray{T, N, ArrayType}(data, m, n)

LowerTriangularMatrix(data::Vector{T}, m::Integer, n::Integer) where T =
    LowerTriangularMatrix{T}(data, m, n)

# SIZE ETC
"""$(TYPEDSIGNATURES)
Length of a `LowerTriangularArray` defined as number of non-zero elements."""
Base.length(L::LowerTriangularArray) = length(L.data)

# super type to be used for ZeroBased (0-based), OneBased (1-based) indexing of the spherical harmonics
abstract type IndexBasis end

"""Abstract type to dispatch for 0-based indexing of the spherical harmonic
degree l and order m, i.e. l=m=0 is the mean, the zonal modes are m=0 etc.
This indexing is more common in mathematics."""
abstract type ZeroBased <: IndexBasis end

"""Abstract type to dispatch for 1-based indexing of the spherical harmonic
degree l and order m, i.e. l=m=1 is the mean, the zonal modes are m=1 etc.
This indexing matches Julia's 1-based indexing for arrays."""
abstract type OneBased <: IndexBasis end

# get matrix size of LTA from its data array and m, n (number of rows and columns)
matrix_size(data::AbstractArray, m::Integer, n::Integer) = (m, n, size(data)[2:end]...)

# extend to get the size of the i-th dimension, with 1 returned for any additional dimension
# as it is also defined for Array
function matrix_size(data::AbstractArray, m::Integer, n::Integer, i::Integer)
    i == 1 && return m      # first dimension is the number of rows m
    i == 2 && return n      # second dimension is the number of columns n
    return size(data, i-1)  # -1 as m, n are collapsed into a vector in the data array
end 

"""$(TYPEDSIGNATURES)
Size of a `LowerTriangularArray` defined as size of the flattened array if `as <: AbstractVector`
and as if it were a full matrix when `as <: AbstractMatrix`` ."""
Base.size(L::LowerTriangularArray, base::Type{<:IndexBasis}=OneBased; as=Vector) = size(L, base, as)
Base.size(L::LowerTriangularArray, i::Integer, base::Type{<:IndexBasis}=OneBased; as=Vector) = size(L, i, base, as)

# use multiple dispatch to chose the right options of basis and vector/flat vs matrix indexing
# matrix indexing can be zero based (natural for spherical harmonics) or one-based,
# vector/flat indexing has only one-based indexing
Base.size(L::LowerTriangularArray, base::Type{OneBased}, as::Type{Matrix}) = matrix_size(L.data, L.m, L.n)
Base.size(L::LowerTriangularArray, base::Type{ZeroBased}, as::Type{Matrix}) = matrix_size(L.data, L.m-1, L.n-1)
Base.size(L::LowerTriangularArray, base::Type{OneBased}, as::Type{Vector}) = size(L.data)

# size(L, i, ...) to get the size of the i-th dimension, with 1 returned for any additional dimension 
Base.size(L::LowerTriangularArray, i::Integer, base::Type{OneBased}, as::Type{Matrix}) = matrix_size(L.data, L.m, L.n, i)
Base.size(L::LowerTriangularArray, i::Integer, base::Type{ZeroBased}, as::Type{Matrix}) = matrix_size(L.data, L.m-1, L.n-1, i)
Base.size(L::LowerTriangularArray, i::Integer, base::Type{OneBased}, as::Type{Vector}) = size(L.data, i)

# sizeof the underlying data vector
Base.sizeof(L::LowerTriangularArray) = sizeof(L.data)

function Base.show(io::IO, ::MIME"text/plain", L::LowerTriangularMatrix)
    Base.array_summary(io, L, axes(L))
    
    # the following is copied over from base/arrayshow.jl
    X = Matrix(L)
    # 2) compute new IOContext
    if !haskey(io, :compact) && length(axes(X, 2)) > 1
        io = IOContext(io, :compact => true)
    end
    if get(io, :limit, false)::Bool && eltype(X) === Method
        # override usual show method for Vector{Method}: don't abbreviate long lists
        io = IOContext(io, :limit => false)
    end

    if get(io, :limit, false)::Bool && displaysize(io)[1]-4 <= 0
        return print(io, " …")
    else
        println(io)
    end

    # 3) update typeinfo
    #
    # it must come after printing the summary, which can exploit :typeinfo itself
    # (e.g. views)
    # we assume this function is always called from top-level, i.e. that it's not nested
    # within another "show" method; hence we always print the summary, without
    # checking for current :typeinfo (this could be changed in the future)
    io = IOContext(io, :typeinfo => eltype(X))

    # 4) show actual content
    recur_io = IOContext(io, :SHOWN_SET => X)
    Base.print_array(recur_io, X)
end

function Base.array_summary(io::IO, L::LowerTriangularMatrix{T}, inds::Tuple{Vararg{Base.OneTo}}) where T
    mn = size(L; as=Matrix)
    print(io, Base.dims2string(length.(inds)), ", $(mn[1])x$(mn[2]) LowerTriangularMatrix{$T}")
end

# CREATE INSTANCES (ZEROS, ONES, UNDEF)
for f in (:zeros, :ones, :rand, :randn)
    @eval begin
        # use ArrayType from LowerTriangularArray parameter
        function Base.$f(
            ::Type{LowerTriangularArray{T, N, ArrayType}}, 
            m::Integer,
            n::Integer,
            I::Vararg{Integer, M},
        ) where {T, N, M, ArrayType}
            ArrayType_ = nonparametric_type(ArrayType)
            return LowerTriangularArray(ArrayType_($f(T, nonzeros(m, n), I...)), m, n)
        end
        
        # default CPU, use Array
        function Base.$f(
            ::Type{LowerTriangularArray{T}},
            m::Integer,
            n::Integer, 
            I::Vararg{Integer, M},
        ) where {T, M}
            return LowerTriangularArray($f(T, nonzeros(m, n), I...), m, n)
        end
        
        # use Float64 as default if type T is not provided
        Base.$f(::Type{LowerTriangularArray}, m::Integer, nk::Integer...) =
            $f(LowerTriangularArray{Float64}, m, nk...)
        Base.$f(::Type{LowerTriangularMatrix}, m::Integer, n::Integer) =
            $f(LowerTriangularArray{Float64}, m, n)
    end
end

Base.zero(L::LTA) where {LTA <: LowerTriangularArray} = zeros(LTA, size(L; as=Matrix)...)
Base.one(L::LTA) where {LTA <: LowerTriangularArray} = ones(LTA, size(L; as=Matrix)...)

function LowerTriangularArray{T, N, ArrayType}(
    ::UndefInitializer,
    I::Vararg{Integer,M},
) where {T, N, M, ArrayType<:AbstractArray{T}}
    return LowerTriangularArray(ArrayType(undef, nonzeros(I[1], I[2]), I[3:end]...), I[1], I[2])
end

function LowerTriangularArray{T, N, ArrayType}(
    ::UndefInitializer,
    size::S,
) where {T, N, ArrayType<:AbstractArray{T}, S<:Tuple}
    return LowerTriangularArray{T, N, ArrayType}(undef, size...)
end
   
function LowerTriangularMatrix{T}(::UndefInitializer, m::Integer, n::Integer) where T
    return LowerTriangularMatrix(Vector{T}(undef, nonzeros(m, n)), m, n)
end

# INDEXING
"""
$(TYPEDSIGNATURES)
Converts the index pair `i, j` of an `m`x`n` LowerTriangularMatrix `L` to a single
index `k` that indexes the same element in the corresponding vector that stores
only the lower triangle (the non-zero entries) of `L`."""
@inline ij2k(i::Integer, j::Integer, m::Integer) = i + (j-1)*m - j*(j-1)÷2
@inline triangle_number(n::Integer) = n*(n+1)÷2
nonzeros(m::Integer, n::Integer) = m*n - triangle_number(n-1)

"""
$(TYPEDSIGNATURES)
Converts the linear index `k` in the lower triangle into a pair `(i, j)` of indices 
of the matrix in column-major form. (Formula taken from 
Angeletti et al, 2019, https://hal.science/hal-02047514/document)
"""
@inline function k2ij(k::Integer, m::Integer) 
    kp = triangle_number(m) - k 
    p = Int(floor((sqrt(1 + 8*kp) - 1)/2))
    i = k - m*(m-1)÷2 + p*(p+1)÷2
    j = m - p
    return i, j
end 
k2ij(I::CartesianIndex, m::Int) = CartesianIndex(k2ij(I[1], m)...,I.I[2:end]...) 
# ij2l(i::Int, j::Int, m::Int) = m*(m-1)÷2 - (m-j)*(m-j-1)÷2 + i - j

# direct indexing, no. indices have to be equal to `N` for the correct dimensionality
@inline Base.getindex(L::LowerTriangularArray{T, N}, I::Vararg{Any, N}) where {T, N} = getindex(L.data, I...) 

# indexing with : + other indices, returns a LowerTriangularArray
@inline function Base.getindex(L::LowerTriangularArray{T,N}, col::Colon, I...) where {T,N}
    return LowerTriangularArray(getindex(L.data, col, I...), L.m, L.n)
end

# l,m sph "matrix-style"  indexing with integer + other indices
@inline function Base.getindex(L::LowerTriangularArray{T,N}, I::Vararg{Any, M}) where {T, N, M}
    i, j = I[1:2]
    @boundscheck (0 < i <= L.m && 0 < j <= L.n) || throw(BoundsError(L, (i, j)))
    # to get a zero element in the correct shape, we just take the zero element of some valid element,
    # there are probably faster ways to do this, but I don't know how, and this is just a fallback anyway 
    @boundscheck j > i && return zero(getindex(L.data, 1, I[3:end]...)) 
    k = ij2k(i, j, L.m)
    return getindex(L.data, k, I[3:end]...)
end

@inline function Base.getindex(L::LowerTriangularMatrix{T}, col::Colon, i::Integer) where T
    if i==1
        return L.data[1:size(L, 1, as=Matrix)]
    else 
        return Matrix(L)[:,i]
    end
end

@inline function Base.getindex(L::LowerTriangularMatrix{T}, i::Integer, col::Colon) where T
    return Matrix(L)[i,:] 
end

# important to do Tuple(I) here for the j > i case as one of the getindex methods above is called
Base.@propagate_inbounds Base.getindex(L::LowerTriangularArray{T,N}, I::CartesianIndex{M}) where {T,N,M} = getindex(L, Tuple(I)...)
Base.@propagate_inbounds Base.getindex(L::LowerTriangularArray{T,N}, I::CartesianIndex{N}) where {T,N} = getindex(L.data, I)

Base.@propagate_inbounds Base.getindex(L::LowerTriangularArray{T,N}, i::Integer) where {T,N} = getindex(L.data, i)

# needed to remove ambigouities in 1D case 
Base.@propagate_inbounds Base.getindex(L::LowerTriangularArray{T,1,V}, i::Integer) where {T,V<:AbstractVector{T}} = getindex(L.data, i)
Base.@propagate_inbounds Base.getindex(L::LowerTriangularArray{T,1,V}, I::CartesianIndex{M}) where {T,V<:AbstractVector{T},M} = getindex(L, Tuple(I)...)
Base.@propagate_inbounds Base.getindex(L::LowerTriangularArray{T,1,V}, i::Integer, I::CartesianIndex{0}) where {T,V<:AbstractVector{T}} = getindex(L, i)
Base.@propagate_inbounds Base.getindex(L::LowerTriangularArray{T,1,V}, i::Integer, I::CartesianIndices{0}) where {T,V<:AbstractVector{T}} = getindex(L, i)

# setindex with lm, ..
@inline Base.setindex!(L::LowerTriangularArray{T,N}, x, I::Vararg{Any, N}) where {T, N} = setindex!(L.data, x, I...)

# setindex with il, im, ..
@inline function Base.setindex!(L::LowerTriangularArray{T,N}, x, I::Vararg{Any, M}) where {T, N, M}
    @boundscheck N+1==M || throw(BoundsError(L, I))
    i, j = I[1:2] 
    @boundscheck i >= j || throw(BoundsError(L, I))
    k = ij2k(i, j, L.m)
    setindex!(L.data, x, k, I[3:end]...)
end

# this is specifically here for the weird way in which AssociatedLegendrePolynomials.jl indexes with an empty CartesianIndex
@inline Base.setindex!(L::LowerTriangularArray{T,N}, x, i::CartesianIndex{0}, I::Vararg{Any, M}) where {T,N, M} = setindex!(L, x, I)


@inline Base.setindex!(L::LowerTriangularArray, x, I::CartesianIndex) = setindex!(L, x, Tuple(I)...)
@inline Base.setindex!(L::LowerTriangularArray{T,N}, x, I::CartesianIndex{N}) where {T,N} = setindex!(L.data, x, I)
@inline Base.setindex!(L::LowerTriangularArray{T,N}, x, i::Integer) where {T,N} = setindex!(L.data, x, i)

# this is to avoid ambigouities
@inline Base.setindex!(L::LowerTriangularArray{T,1,V}, x, i::Integer) where {T,V<:AbstractVector{T}} = setindex!(L.data, x, i)
@inline Base.setindex!(L::LowerTriangularArray{T,1,V}, x, I::CartesianIndex{1}) where {T,V<:AbstractVector{T}} = setindex!(L.data, x, I)
@inline Base.setindex!(L::LowerTriangularArray{T,1,V}, x, I::CartesianIndex{2}) where {T,V<:AbstractVector{T}} = setindex!(L, x, Tuple(I)...)
@inline Base.setindex!(L::LowerTriangularArray{T,1,V}, x, i::Integer, I::CartesianIndex{0}) where {T,V<:AbstractVector{T}} = setindex!(L, x, i)

# propagate index to data vector
Base.eachindex(L ::LowerTriangularArray)    = eachindex(L.data)
Base.eachindex(Ls::LowerTriangularArray...) = eachindex((L.data for L in Ls)...)

"""
$(TYPEDSIGNATURES)
creates `unit_range::UnitRange` to loop over all non-zeros/spherical harmonics numbers in a LowerTriangularArray `L`.
Like `eachindex` but skips the upper triangle with zeros in `L`."""
eachharmonic(L::LowerTriangularArray) = axes(L.data, 1) 

"""
$(TYPEDSIGNATURES)
creates `unit_range::UnitRange` to loop over all non-zeros in the LowerTriangularMatrices
provided as arguments. Checks bounds first. All LowerTriangularMatrix's need to be of the same size.
Like `eachindex` but skips the upper triangle with zeros in `L`."""
function eachharmonic(L1::LowerTriangularArray, Ls::LowerTriangularArray...) 
    lowertriangular_match(L1, Ls...; horizontal_only=true) || throw(DimensionMismatch(L1, Ls...))
    return eachharmonic(L1) 
end

"""$(TYPEDSIGNATURES) Iterator for the non-horizontal dimensions in
LowerTriangularArrays. To be used like

    for k in eachmatrix(L)
        L[1, k]
    
to loop over every non-horizontal dimension of L."""
eachmatrix(L::LowerTriangularArray) = CartesianIndices(size(L)[2:end])

"""$(TYPEDSIGNATURES) Iterator for the non-horizontal dimensions in
LowerTriangularArrays. Checks that the LowerTriangularArrays match according to
`lowertriangular_match`."""
function eachmatrix(L1::LowerTriangularArray, Ls::LowerTriangularArray...)
    lowertriangular_match(L1, Ls...) || throw(DimensionMismatch(L1, Ls...))
    return eachmatrix(L1)
end

"""$(TYPEDSIGNATURES) True if both `L1` and `L2` are of the same size (as matrix),
but ignores singleton dimensions, e.g. 5x5 and 5x5x1 would match.
With `horizontal_only=true` (default `false`) ignore the non-horizontal dimensions,
e.g. 5x5, 5x5x1, 5x5x2 would all match."""
function lowertriangular_match(
    L1::LowerTriangularArray,
    L2::LowerTriangularArray;
    horizontal_only::Bool=false,
)
    horizontal_match = size(L1, as=Matrix)[1:2] == size(L2, as=Matrix)[1:2]
    horizontal_only && return horizontal_match
    return horizontal_match && length(L1) == length(L2)     # ignores singleton dimensions
end


"""$(TYPEDSIGNATURES) True if all lower triangular matrices provided as arguments
match according to `lowertriangular_match` wrt to `L1` (and therefore all)."""
function lowertriangular_match(L1::LowerTriangularArray, Ls::LowerTriangularArray...; kwargs...)
    length(Ls) == 0 && return true   # single L1 always matches itself
    # cut the Ls tuple short on every iteration of the recursion
    return lowertriangular_match(L1, Ls[1]; kwargs...) && lowertriangular_match(L1, Ls[2:end]...; kwargs...)
end

function Base.DimensionMismatch(L1::LowerTriangularArray, Ls::LowerTriangularArray...)
    s = "LowerTriangularArrays do not match; $(Base.dims2string(size(L1, as=Matrix)))"
    for L in Ls
        s *= ", $(Base.dims2string(size(L, as=Matrix)))"
    end
    return DimensionMismatch(s)
end

# CONVERSIONS
""" 
$(TYPEDSIGNATURES)
Create a LowerTriangularArray `L` from Matrix `M` by copying over the non-zero elements in `M`."""
function LowerTriangularMatrix(M::Matrix{T}) where T # CPU version 
    m, n = size(M)
    L = LowerTriangularMatrix{T}(undef, size(M)...)
    
    k = 0
    @inbounds for j in 1:n      # only loop over lower triangle
        for i in j:m
            k += 1              # next element in lower triangle
            L[k] = M[i, j]       # and copy data into vector
        end
    end
    return L
end

# helper function for conversion etc on GPU, returns indices of the lower triangle
lowertriangle_indices(M::AbstractMatrix) = lowertriangle_indices(size(M)...)
lowertriangle_indices(m::Integer, n::Integer) = tril!(trues((m,n)))

function lowertriangle_indices(M::AbstractArray{T, N}) where {T, N}
    @boundscheck N >= 3 || throw(BoundsError)

    indices = lowertriangle_indices(getindex(M,:,:,[axes(M,3)[1] for i=3:N]...)) 
    indices = reshape(indices, size(indices,1), size(indices,2), [1 for i=1:(N-2)]...)
    
    repeat(indices, 1, 1, size(M)[3:end]...)
end  

# GPU version and fallback for higher dimensions
"""
$(TYPEDSIGNATURES)
Create a LowerTriangularArray `L` from Array `M` by copying over the non-zero elements in `M`.
"""
function LowerTriangularArray(M::ArrayType) where {T, N, ArrayType<:AbstractArray{T,N}} 
    m, n = size(M)[1:2]
    LowerTriangularArray(M[lowertriangle_indices(m,n),[Colon() for i=1:N-2]...], m, n)
end

LowerTriangularArray(L::LowerTriangularArray) = LowerTriangularArray(L.data, L.m, L.n)
LowerTriangularMatrix(M::AbstractMatrix) = LowerTriangularArray(M)

function Base.Matrix(L::LowerTriangularMatrix{T}) where T
    m, n = size(L, as=Matrix)
    M = zeros(T, m, n)
    copyto!(M, L)
end

function Base.Array(L::LowerTriangularArray{T}) where T
    A = zeros(T, size(L, as=Matrix))
    copyto!(A, L)
end 
            
function Base.copyto!(
    L1::LowerTriangularArray{T},
    L2::LowerTriangularArray
) where T
    # if sizes don't match copy over the largest subset of indices
    size(L1) != size(L2) && return copyto!(L1, L2,  Base.OneTo(minimum(size.((L1, L2), 1; as=Matrix))),
                                                    Base.OneTo(minimum(size.((L1, L2), 2; as=Matrix))))

    L1.data .= convert.(T, L2.data)
    return L1
end

# CPU version
function Base.copyto!(  
    L1::LowerTriangularArray{T,N,ArrayType1},   # copy to L1
    L2::LowerTriangularArray{S,N,ArrayType2},   # copy from L2
    ls::AbstractUnitRange,                      # range of indices in 1st dim
    ms::AbstractUnitRange,                      # range of indices in 2nd dim
) where {T, S, N, ArrayType1<:Array{T,N}, ArrayType2<:Array{S,N}}

    lmax, mmax = size(L2; as=Matrix)            # use the size of L2 for boundscheck
    @boundscheck maximum(ls) <= lmax || throw(BoundsError)
    @boundscheck maximum(ms) <= mmax || throw(BoundsError)

    lmax, mmax = size(L1; as=Matrix)            # but the size of L1 to loop
    lm = 0
    @inbounds for m in 1:mmax
        for l in m:lmax
            lm += 1
            L1[lm,[Colon() for i=1:(N-1)]...] = (l in ls) && (m in ms) ?
                                                T.(L2[l, m, [Colon() for i=1:(N-1)]...]) :
                                                L1[lm, [Colon() for i=1:(N-1)]...]
        end
    end

    return L1
end 

# Fallback / GPU version (the two versions _copyto! and copyto! are there to enable tests of this function with regular Arrays)
function Base.copyto!(
    L1::LowerTriangularArray{T, N, ArrayType1},
    L2::LowerTriangularArray{S, N, ArrayType2},
    ls::AbstractUnitRange,
    ms::AbstractUnitRange
) where {T, S, N, ArrayType1<:AbstractArray{T}, ArrayType2<:AbstractArray{S}}
    return _copyto_core!(L1, L2, ls, ms)
end

function _copyto_core!( 
    L1::LowerTriangularArray{T,N,ArrayType1},   # copy to L1
    L2::LowerTriangularArray{S,N,ArrayType2},   # copy from L2
    ls::AbstractUnitRange,                      # range of indices in 1st dim
    ms::AbstractUnitRange,                      # range of indices in 2nd dim
) where {T, S, N, ArrayType1<:AbstractArray{T}, ArrayType2<:AbstractArray{S}}

    lmax_2, mmax_2 = size(L2, as=Matrix)     # use the size of L2 for boundscheck
    @boundscheck maximum(ls) <= lmax_2 || throw(BoundsError)
    @boundscheck maximum(ms) <= mmax_2 || throw(BoundsError)

    lmax_1, mmax_1 = size(L1, as=Matrix)

    # we'll setup the 1D-indices that correspond to the given range here
    ind_L1 = lowertriangle_indices(lmax_1, mmax_1)
    ind_L1_smaller = deepcopy(ind_L1) # threshold the lower triangular based on the ls and ms
    if maximum(ls) < lmax_1 # threshold L1 based on ls 
        ind_L1_smaller[ls[end]+1:end, :] .= 0
    end 
    if maximum(ms) < mmax_1 # threshold L1 based on ms
        ind_L1_smaller[:, ms[end]+1:end] .= 0 
    end 
    ind_L1 = ind_L1_smaller[ind_L1] # also flattens the indices into indices for the L1.data array

    ind_L2 = lowertriangle_indices(lmax_2, mmax_2)
    ind_L2_smaller = deepcopy(ind_L2) # threshold the lower triangular based on the ls and ms and L1,L2
    if maximum(ls) > lmax_1 # so that L2 can fit into L1 
        ind_L2_smaller[lmax_1+1:end, :] .= 0 
    end 
    if maximum(ls) < lmax_2 # threshold L2 based on ls
        ind_L2_smaller[ls[end]+1:end, :] .= 0
    end 
    if maximum(ms) > mmax_1 # so that L2 can fit into L1 
        ind_L2_smaller[:, mmax_1+1:end] .= 0 
    end
    if maximum(ms) < mmax_2 # threshold L2 based on ms
        ind_L2_smaller[:, ms[end]+1:end] .= 0 
    end 
    ind_L2 = ind_L2_smaller[ind_L2] # also flattens the indices into indices for the L2.data array

    L1.data[ind_L1,[Colon() for i=1:(N-1)]...] = T.(L2.data[ind_L2,[Colon() for i=1:(N-1)]...])

    return L1
end 


# copyto! using matrix indexing from Matrix/Array
function Base.copyto!(  L::LowerTriangularArray{T},  # copy to L
                        M::AbstractArray) where T    # copy from M
    @boundscheck size(L, as=Matrix) == size(M) || throw(BoundsError)
    L.data .= convert.(T, M[lowertriangle_indices(M)])
    return L
end

function Base.copyto!(  M::AbstractArray{T},               # copy to M
                        L::LowerTriangularArray) where T   # copy from L
  
    @boundscheck size(L, as=Matrix) == size(M) || throw(BoundsError)

    lower_triangle_indices = lowertriangle_indices(M)
    upper_triangle_indices = @. ~lower_triangle_indices # upper triangle without main diagonal 

    M[upper_triangle_indices] .= zero(T)
    M[lower_triangle_indices] = convert.(T, L.data)

    return M
end

# copyto! from Vector/Array to using vector indexing
function Base.copyto!(  L::LowerTriangularArray{T,N},       # copy to L
                        V::AbstractArray{S,N}) where {T,S,N}# copy from V
    @boundscheck size(L, as=Vector) == size(V) || throw(BoundsError)
    L.data .= convert.(T, V)
    return L 
end 

function LowerTriangularMatrix{T}(M::LowerTriangularMatrix{T2}) where {T,T2}
    L = LowerTriangularMatrix{T}(undef, size(M)...)
    copyto!(L, M)
    return L
end

function Base.convert(
    ::Type{LowerTriangularArray{T1, N, ArrayTypeT1}},
    L::LowerTriangularArray{T2, N, ArrayTypeT2},
) where {T1, T2, N, ArrayTypeT1<:AbstractArray{T1}, ArrayTypeT2<:AbstractArray{T2}}
    return LowerTriangularArray{T1, N, ArrayTypeT1}(L.data, L.m, L.n)
end

function Base.convert(::Type{LowerTriangularMatrix{T}}, L::LowerTriangularMatrix) where T
    return LowerTriangularMatrix{T}(L.data, L.m, L.n)
end

function Base.similar(::LowerTriangularArray{T, N, ArrayType}, I::Integer...) where {T, N, ArrayType}
    return LowerTriangularArray{T,N,ArrayType}(undef, I...)
end

function Base.similar(::LowerTriangularArray{T, N, ArrayType}, size::S) where {T, N, ArrayType, S<:Tuple}
    return LowerTriangularArray{T,N,ArrayType}(undef, size...)
end

function Base.similar(L::LowerTriangularArray{S, N, ArrayType}, ::Type{T}) where {T, S, N, ArrayType}
    ArrayType_ = nonparametric_type(ArrayType) # TODO: not sure how else to infer this type 
    return LowerTriangularArray{T, N, ArrayType_{T, N}}(undef, size(L; as=Matrix)...)
end

Base.similar(L::LowerTriangularArray{T, N, ArrayType}, ::Type{T}) where {T, N, ArrayType} =
    LowerTriangularArray{T, N, ArrayType}(undef, size(L; as=Matrix)...)
Base.similar(L::LowerTriangularArray{T}) where T = similar(L, T)
 
Base.prod(L::LowerTriangularArray{NF}) where NF = zero(NF)

"""
$(TYPEDSIGNATURES)
Fills the elements of `L` with `x`. Faster than fill!(::AbstractArray, x)
as only the non-zero elements in `L` are assigned with x."""
function Base.fill!(L::LowerTriangularArray, x)
    fill!(L.data, x)
    return L
end

Base.:(==)(L1::LowerTriangularArray, L2::LowerTriangularArray) = 
    L1.m == L2.m && L1.n == L2.n && L1.data == L2.data
Base.isapprox(L1::LowerTriangularArray, L2::LowerTriangularArray; kwargs...) =
    isapprox(L1.data, L2.data; kwargs...)
Base.all(L::LowerTriangularArray) = all(L.data)
Base.any(L::LowerTriangularArray) = any(L.data)

Base.repeat(L::LowerTriangularArray, counts...) = LowerTriangularArray(repeat(L.data, counts...), L.m, L.n)

# Broadcast CPU/GPU
import Base.Broadcast: BroadcastStyle, Broadcasted, DefaultArrayStyle
import LinearAlgebra: isstructurepreserving, fzeropreserving

# CPU with scalar indexing
struct LowerTriangularStyle{N, ArrayType} <: Broadcast.AbstractArrayStyle{N} end

# GPU without scalar indexing
struct LowerTriangularGPUStyle{N, ArrayType} <: GPUArrays.AbstractGPUArrayStyle{N} end

function BroadcastStyle(::Type{LowerTriangularArray{T, N, ArrayType}}) where {T, N, ArrayType}
    # remove number format parameter for broadcasting with type promotion
    ArrayType_ = nonparametric_type(ArrayType)
    return LowerTriangularStyle{N, ArrayType_}()
end

function BroadcastStyle(
    ::Type{LowerTriangularArray{T, N, ArrayType}},
) where {T, N, ArrayType <: GPUArrays.AbstractGPUArray}
    # remove number format parameter for broadcasting with type promotion
    ArrayType_ = nonparametric_type(ArrayType)
    return LowerTriangularGPUStyle{N, ArrayType_}()
end

# ::Val{0} for broadcasting with 0-dimensional, ::Val{1} for broadcasting with vectors, etc
LowerTriangularStyle{N, ArrayType}(::Val{M}) where {N, ArrayType, M} =
    LowerTriangularStyle{N, ArrayType}()
LowerTriangularGPUStyle{N, ArrayType}(::Val{M}) where {N, ArrayType, M} =
    LowerTriangularGPUStyle{N, ArrayType}()

# also needed for other array types
nonparametric_type(::Type{<:Array}) = Array

"`L = find_L(Ls)` returns the first LowerTriangularArray among the arguments. 
Adapted from Julia documentation of Broadcast interface"
find_L(bc::Base.Broadcast.Broadcasted) = find_L(bc.args)
find_L(args::Tuple) = find_L(find_L(args[1]), Base.tail(args))
find_L(x) = x
find_L(::Tuple{}) = nothing
find_L(a::LowerTriangularArray, rest) = a
find_L(::Any, rest) = find_L(rest)

function Base.similar(
    bc::Broadcasted{LowerTriangularStyle{N, ArrayType}},
    ::Type{T},
) where {N, ArrayType, T}
    L = find_L(bc)
    return LowerTriangularArray{T, N, ArrayType{T,N}}(undef, size(L; as=Matrix))
end

# same function as above, but needs to be defined for both CPU and GPU style
function Base.similar(
    bc::Broadcasted{LowerTriangularGPUStyle{N, ArrayType}},
    ::Type{T},
) where {N, ArrayType, T}
    L = find_L(bc)
    return LowerTriangularArray{T, N, ArrayType{T,N}}(undef, size(L; as=Matrix))
end

function GPUArrays.backend(
    ::Type{LowerTriangularArray{T, N, ArrayType}}
) where {T, N, ArrayType <: GPUArrays.AbstractGPUArray}
    return GPUArrays.backend(ArrayType)
end

Adapt.adapt_structure(to, L::LowerTriangularArray) =
    LowerTriangularArray(Adapt.adapt(to, L.data), L.m, L.n)