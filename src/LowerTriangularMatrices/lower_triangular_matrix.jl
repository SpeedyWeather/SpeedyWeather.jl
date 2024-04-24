"""
$(TYPEDSIGNATURES)
A lower triangular array implementation that only stores the non-zero entries explicitly.
`L<:AbstractMatrix` although in general we have `L[i] != Matrix(L)[i]`, the former skips
zero entries, tha latter includes them.

Supports n-dimensional lower triangular arrays, so that for all trailing dimensions `L[:, :, ..]`
is a matrix in lower triangular form, e.g. a (5x5x3)-LowerTriangularArray would hold 3 lower
triangular matrices."""
struct LowerTriangularArray{T, N, ArrayType <: AbstractArray{T}} <: AbstractArray{T,N}
    data::ArrayType     # non-zero elements unravelled into an array with one dimension less N-1
    m::Int              # number of rows
    n::Int              # number of columns

    LowerTriangularArray{T, N, ArrayType}(data, m, n) where {T, N, ArrayType<:AbstractArray{T}} =
        check_lta_input_array(data, m, n, N) ?
        new(data, m, n) :
        error(lta_error_message(data, m, n, T, N, ArrayType))
end

check_lta_input_array(data, m, n, N) =
    (ndims(data) == N-1) & (length(data) == prod(size(data)[2:end]) * nonzeros(m, n)) 

function lta_error_message(data, m, n, T, N, ArrayType) 
    size_tuple = (m, n, size(data[2:end])...)
    return "$(size(data))-sized $(typeof(data)) cannot be used to create "*
                "a $size_tuple LowerTriangularArray{$T,$N,$ArrayType}"
end 

"""2-dimensional `LowerTriangularArray` of type `T`` with its non-zero entries unravelled into a `Vector{T}`"""
const LowerTriangularMatrix{T} = LowerTriangularArray{T, 2, Vector{T}}

LowerTriangularArray(data::ArrayType, m::Integer, n::Integer) where {T, N, ArrayType<:AbstractArray{T,N}} =
    LowerTriangularArray{T, N+1, ArrayType}(data, m, n)
LowerTriangularMatrix(data::Vector{T}, m::Integer, n::Integer) where T =
    LowerTriangularMatrix{T}(data, m, n)

# SIZE ETC
"""$(TYPEDSIGNATURES)
Length of a `LowerTriangularArray` defined as number of non-zero elements."""
Base.length(L::LowerTriangularArray) = length(L.data)

"""$(TYPEDSIGNATURES)
Size of a `LowerTriangularArray` defined as matrix size including zero upper triangle."""
Base.size(L::LowerTriangularArray) = (L.m, L.n, size(L.data)[2:end]...)

# sizeof the underlying data vector
Base.sizeof(L::LowerTriangularArray) = sizeof(L.data)

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
            return LowerTriangularArray(ArrayType($f(T, nonzeros(m, n), I...)), m, n)
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

Base.zero(L::LTA) where {LTA <: LowerTriangularArray} = zeros(LTA, size(L)...)
Base.one(L::LTA) where {LTA <: LowerTriangularArray} = ones(LTA, size(L)...)

function LowerTriangularArray{T, N, ArrayType}(
    ::UndefInitializer,
    I::Vararg{Integer,N},
) where {T, N, ArrayType<:AbstractArray{T}}
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
    (k - m*(m-1)÷2 + p*(p+1)÷2, m - p)
end 
k2ij(I::CartesianIndex, m::Int) = CartesianIndex(k2ij(I[1], m)...,I.I[2:end]...) 

# direct indexing, no. indices have to be one less than `N` for the correct dimensionality, so N!=M
@inline function Base.getindex(L::LowerTriangularArray{T, N}, I::Vararg{Integer, M}) where {T, N, M}
    @boundscheck M == N-1 || throw(BoundsError(L, I))
    return getindex(L.data, I...) 
end

# integer + other indices (:, ranges, etc...)
@inline function Base.getindex(L::LowerTriangularArray{T, N}, k::Integer, I::Vararg{R, M}) where {T,N,R,M}
    @boundscheck M == N-2 || throw(BoundsError(L, I))
    return getindex(L.data, k, I...)
end

# double index i,j, i.e. spherical harmoinc indexing (1-based), no. indices has to be equal to N 
Base.@propagate_inbounds function Base.getindex(L::LowerTriangularArray{T, N}, I::Vararg{Integer, N}) where {T, N}
    i, j = I[1:2]
    @boundscheck (0 < i <= L.m && 0 < j <= L.n) || throw(BoundsError(L, (i, j)))
    @boundscheck j > i && return zero(T)
    k = ij2k(i, j, L.m)
    return getindex(L.data, k, I[3:end]...)
end

# l,m sph indexing with integer + other indices
@inline function Base.getindex(L::LowerTriangularArray{T,N}, i::Integer, j::Integer, I::Vararg{R, M}) where {T, N, R, M}
    @boundscheck M == N-2 || throw(BoundsError(L, I))
    @boundscheck (0 < i <= L.m && 0 < j <= L.n) || throw(BoundsError(L, (i, j)))
    # to get a zero element in the correct shape, we just take the zero element of some valid element,
    # there are probably faster ways to do this, but I don't know how, and this is just a fallback anyway 
    @boundscheck j > i && return zero(getindex(L.data, 1, I...)) 
    k = ij2k(i, j, L.m)
    return getindex(L.data, k, I...)
end

# indexing with : + other indices, returns a LowerTriangularArray
@inline function Base.getindex(L::LowerTriangularArray{T,N}, col::Colon, I...) where {T,N}
    return LowerTriangularArray(getindex(L.data, col, I...), L.m, L.n)
end

@inline function Base.getindex(L::LowerTriangularMatrix{T}, col::Colon, i::Integer) where T
    return Matrix(L)[:,i]
end

# Base.@propagate_inbounds Base.getindex(L::LowerTriangularArray, r::AbstractRange) = getindex(L.data, r)
Base.@propagate_inbounds Base.getindex(L::LowerTriangularArray, r::AbstractRange, I...) = getindex(L.data, r, I...)
Base.@propagate_inbounds Base.getindex(L::LowerTriangularArray, i::Integer) = getindex(L.data, i)

# important to do Tuple(I) here for the j > i case as one of the getindex methods above is called
Base.@propagate_inbounds Base.getindex(L::LowerTriangularArray, I::CartesianIndex) = getindex(L, Tuple(I)...)

# setindex with lm, ..
@inline function Base.setindex!(L::LowerTriangularArray{T,N}, x, I::Vararg{Any, M}) where {T, N, M} 
    @boundscheck M == N-1 || throw(BoundsError(L, I))
    setindex!(L.data, x, I...)
end 

# setindex with il, im, ..
@inline function Base.setindex!(L::LowerTriangularArray{T,N}, x, I::Vararg{Any, N}) where {T, N}
    i, j = I[1:2] # TODO: check if i,j are Int?
    @boundscheck i >= j || throw(BoundsError(L, (i, j)))
    k = ij2k(i, j, L.m)
    setindex!(L.data, x, k, I[3:end]...)
end

# this is specifically here for the weird way in which AssociatedLegendrePolynomials.jl indexes with an empty CartesianIndex
@inline function Base.setindex!(L::LowerTriangularArray{T,N}, x, i::CartesianIndex, I::Vararg{Any, N}) where {T,N}
    i, j = I[1:2] # TODO: check if i,j are Int?
    @boundscheck i >= j || throw(BoundsError(L, (i, j)))
    k = ij2k(i, j, L.m)
    setindex!(L.data, x, k, I[3:end]...)
end

@inline Base.setindex!(L::LowerTriangularArray, x::AbstractVector, r::AbstractRange, I...) = setindex!(L.data, x, r, I...)
@inline Base.setindex!(L::LowerTriangularArray, x, i::Integer) = setindex!(L.data, x, i)
@inline Base.setindex!(L::LowerTriangularArray, x, I::CartesianIndex) = setindex!(L, x, Tuple(I)...)

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
    n = size(L1.data,1) 
    Base._all_match_first(L->size(L.data,1), n, L1, Ls...) || throw(BoundsError)
    return eachharmonic(L1) 
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
lowertriangle_indices(M::AbstractMatrix) = tril!(trues(size(M)))
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
    m, n = size(L)
    M = zeros(T, m, n)
    copyto!(M, L)
end
            
function Base.copyto!(
    L1::LowerTriangularArray{T, N, ArrayType1},
    L2::LowerTriangularArray
) where {T, N, ArrayType1<:AbstractArray{T}}
    # if sizes don't match copy over the largest subset of indices
    size(L1) != size(L2) && return copyto!(L1, L2,  Base.OneTo(minimum(size.((L1, L2), 1))),
                                                    Base.OneTo(minimum(size.((L1, L2), 2))))

    L1.data .= convert.(T, L2.data)
    L1
end

# CPU version
function Base.copyto!(  
    L1::LowerTriangularArray{T,N,ArrayType1},   # copy to L1
    L2::LowerTriangularArray{S,N,ArrayType2},   # copy from L2
    ls::AbstractUnitRange,                      # range of indices in 1st dim
    ms::AbstractUnitRange,                      # range of indices in 2nd dim
) where {T, S, N, ArrayType1<:Array{T}, ArrayType2<:Array{S}}

    lmax, mmax = size(L2)        # use the size of L2 for boundscheck
    @boundscheck maximum(ls) <= lmax || throw(BoundsError)
    @boundscheck maximum(ms) <= mmax || throw(BoundsError)

    lmax, mmax = size(L1)        # but the size of L1 to loop
    lm = 0
    @inbounds for m in 1:mmax
        for l in m:lmax
            lm += 1
            L1[lm,[Colon() for i=1:(N-2)]...] = (l in ls) && (m in ms) ?
                                                T.(L2[l, m, [Colon() for i=1:(N-2)]...]) :
                                                L1[lm, [Colon() for i=1:(N-2)]...]
        end
    end

    L1
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

    lmax_2, mmax_2 = size(L2)[1:2]       # use the size of L2 for boundscheck
    @boundscheck maximum(ls) <= lmax_2 || throw(BoundsError)
    @boundscheck maximum(ms) <= mmax_2 || throw(BoundsError)

    lmax_1, mmax_1 = size(L1)[1:2]

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

    L1.data[ind_L1,[Colon() for i=1:(N-2)]...] = T.(L2.data[ind_L2,[Colon() for i=1:(N-2)]...])

    L1
end 

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
    upper_triangle_indices = @. ~lower_triangle_indices # upper triangle without main diagonal 

    M[upper_triangle_indices] .= zero(T)
    M[lower_triangle_indices] = convert.(T, L.data)

    M
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
    return LowerTriangularArray{T1,N,ArrayTypeT1}(L.data, L.m, L.n)
end

function Base.convert(::Type{LowerTriangularMatrix{T}}, L::LowerTriangularMatrix) where T
    return LowerTriangularMatrix{T}(L.data, L.m, L.n)
end

function Base.similar(::LowerTriangularArray{T,N,ArrayType}, size::S) where {T, N, ArrayType, S<:Tuple}
    return LowerTriangularArray{T,N,ArrayType}(undef, size...)
end

function Base.similar(L::LowerTriangularArray{S,N,ArrayType}, ::Type{T}) where {T, S, N, ArrayType}
    new_array_type = typeof(T.(L.data)) # TODO: not sure how else to infer this type 
    return LowerTriangularArray{T,N,new_array_type}(undef, size(L)...)
end

Base.similar(L::LowerTriangularArray{T,N,ArrayType}, ::Type{T}) where {T, N, ArrayType} =
    LowerTriangularArray{T,N,ArrayType}(undef, size(L)...)
Base.similar(L::LowerTriangularArray{T}) where T = similar(L, T)
 
# ARITHMETIC
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
    typeof(L1) == typeof(L2) && L1.data == L2.data
Base.isapprox(L1::LowerTriangularArray, L2::LowerTriangularArray; kwargs...) =
    isapprox(L1.data, L2.data; kwargs...)
Base.all(L::LowerTriangularArray) = all(L.data)
Base.any(L::LowerTriangularArray) = any(L.data)

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
# nonparametric_type(::Type{<:CUDA.CuArray}) = CuArray

function Base.similar(
    bc::Broadcasted{LowerTriangularStyle{N, ArrayType}},
    ::Type{T},
) where {N, ArrayType, T}
    ArrayType_ = nonparametric_type(ArrayType)
    #TODO branch currently not hit because these are always false, although they should be true for
    # operations like +, -, ... then returning a `LowerTriangularArray` again
    # if isstructurepreserving(bc) || fzeropreserving(bc)
    #     return LowerTriangularArray{T, N, ArrayType_{T}}(undef, size(bc)...)
    # end
    # TODO should return the similar for operations that escape the structure, e.g. .==
    # but given the TODO above return a `LowerTriangularArray` in both cases
    # return similar(convert(Broadcasted{DefaultArrayStyle{ndims(bc)}}, bc), T)
    return LowerTriangularArray{T, N, ArrayType_{T}}(undef, size(bc)...)
end

# same function as above, but needs to be defined for both CPU and GPU style
function Base.similar(
    bc::Broadcasted{LowerTriangularGPUStyle{N, ArrayType}},
    ::Type{T},
) where {N, ArrayType, T}
    ArrayType_ = nonparametric_type(ArrayType)
    #TODO branch currently not hit because these are always false, although they should be true for
    # operations like +, -, ... then returning a `LowerTriangularArray` again
    # if isstructurepreserving(bc) || fzeropreserving(bc)
    #     return LowerTriangularArray{T, N, ArrayType_{T}}(undef, size(bc)...)
    # end
    # TODO should return the similar for operations that escape the structure, e.g. .==
    # but given the TODO above return a `LowerTriangularArray` in both cases
    # return similar(convert(Broadcasted{DefaultArrayStyle{ndims(bc)}}, bc), T)
    return LowerTriangularArray{T, N, ArrayType_{T}}(undef, size(bc)...)
end

function Base.copyto!(
    dest::LowerTriangularArray{T, N, ArrayTypeP},
    bc::Broadcasted{LowerTriangularStyle{N, ArrayType}},
) where {T, N, ArrayTypeP <: Array, ArrayType <: Array}
    axs = axes(dest)
    axes(bc) == axs || Broadcast.throwdm(axes(bc), axs)

    lmax, mmax = size(dest)
    lm = 0

    for m in 1:mmax
        for l in m:lmax
            lm += 1
            for I in Iterators.product(axs[3:end]...)
                dest.data[lm, I...] = Broadcast._broadcast_getindex(bc, CartesianIndex(l, m, I...))
            end 
        end
    end
    return dest
end

# GPU methods (should be moved to an extension, in case this becomes a standalone package)
import GPUArrays._copyto!

# copied and adjusted from GPUArrays.jl
@inline function GPUArrays._copyto!(
    dest::LowerTriangularArray{T, N, ArrayType},
    bc::Broadcasted
) where {T, N, ArrayType}
    
    axes(dest) == axes(bc) || Broadcast.throwdm(axes(dest), axes(bc))
    isempty(dest) && return dest
    bc = Broadcast.preprocess(dest, bc)

    broadcast_kernel = function (ctx, dest, bc, nelem)
            i = 0
            while i < nelem
                i += 1
                I = GPUArrays.@cartesianidx(dest.data, i) 
                # TODO: the k2ij is costly, if we intend to accelarate this, this needs to be precomputed
                @inbounds dest.data[I] = bc[k2ij(I, dest.m)]
            end                                                                     
            return
        end

    elements = length(dest.data)
    elements_per_thread = typemax(Int)
    heuristic = GPUArrays.launch_heuristic(GPUArrays.backend(dest), broadcast_kernel, dest, bc, 1;
                                 elements, elements_per_thread)
    config = GPUArrays.launch_configuration(GPUArrays.backend(dest), heuristic;
                                  elements, elements_per_thread)
    GPUArrays.gpu_call(broadcast_kernel, dest, bc, config.elements_per_thread;
             threads=config.threads, blocks=config.blocks)

    return dest
end

function GPUArrays.backend(
    ::Type{LowerTriangularArray{T, N, ArrayType}}
) where {T, N, ArrayType <: GPUArrays.AbstractGPUArray}
    return GPUArrays.backend(ArrayType)
end

Adapt.adapt_structure(to, L::LowerTriangularArray) =
    LowerTriangularArray(Adapt.adapt(to, L.data), L.m, L.n)