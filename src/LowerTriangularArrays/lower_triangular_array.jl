"""
$(TYPEDSIGNATURES)
A lower triangular array implementation that only stores the non-zero entries explicitly.
`L<:AbstractArray{T,N-1}` although we do allow both "flat" `N-1`-dimensional indexing and 
additional `N`-dimensional or "matrix-style" indexing.

Supports n-dimensional lower triangular arrays, so that for all trailing dimensions `L[:, :, ..]`
is a matrix in lower triangular form, e.g. a (5x5x3)-LowerTriangularArray would hold 3 lower
triangular matrices."""
struct LowerTriangularArray{T, N, ArrayType <: AbstractArray{T, N}, S<:AbstractSpectrum} <: AbstractArray{T, N}
    data::ArrayType     # non-zero elements unravelled into an in which the lower triangle is flattened
    spectrum::S         # spectrum, that holds all spectral discretization information and the architecture the array is on

    LowerTriangularArray{T, N, ArrayType, S}(data::AbstractArray, spectrum::S) where {T, N, ArrayType <: AbstractArray{T, N}, S <: AbstractSpectrum} =
        check_lta_input_array(data, spectrum, N) ? 
        new(data, spectrum)  :
        error(lta_error_message(data, spectrum, T, N, ArrayType, S))
end

check_lta_input_array(data, spectrum, N) =
    (ndims(data) == N) & (length(data) == prod(size(data)[2:end]) * nonzeros(spectrum)) # & ismatching(spectrum, typeof(data)) TODO: reactivate this? problem for constructors

function lta_error_message(data, spectrum, T, N, ArrayType, S) 
    return "$(Base.dims2string(size(data))) $(typeof(data)) cannot be used to create "*
            "a $(Base.dims2string(matrix_size(data, spectrum))) LowerTriangularArray{$T, $N, $ArrayType, $S}"
end 

"""2D `LowerTriangularArray` of type `T`"""
const LowerTriangularMatrix = LowerTriangularArray{T, 1} where T

# construct LTA from data and spectrum
LowerTriangularArray( data::ArrayType, spectrum::S) where {T, N, ArrayType <: AbstractArray{T, N}, S <: AbstractSpectrum} =
    LowerTriangularArray{T, N, ArrayType, S}(data, spectrum)

# or construct using LowerTriangularMatrix for 1D arrays
LowerTriangularMatrix(data::ArrayType, spectrum::S) where {T,    ArrayType <: AbstractArray{T, 1}, S <: AbstractSpectrum} =
    LowerTriangularArray{T, 1, ArrayType, S}(data, spectrum)

function LowerTriangularArray(data::ArrayType, lmax::Integer, mmax::Integer) where {T, N, ArrayType <: AbstractArray{T, N}}
    spectrum = Spectrum(lmax, mmax, architecture=architecture(ArrayType))
    return LowerTriangularArray{T, N, ArrayType, typeof(spectrum)}(data, spectrum)
end

function LowerTriangularMatrix(data::ArrayType, lmax::Integer, mmax::Integer) where {T, ArrayType <: AbstractArray{T, 1}}
    spectrum = Spectrum(lmax, mmax) 
    return LowerTriangularArray{T, 1, ArrayType, typeof(spectrum)}(data, spectrum)
end 

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

# get matrix size of LTA from its data array and lmax, mmax (number of rows and columns)
matrix_size(data::AbstractArray, spectrum::AbstractSpectrum) = matrix_size(data, spectrum.lmax, spectrum.mmax)
matrix_size(data::AbstractArray, spectrum::AbstractSpectrum, i::Integer) = matrix_size(data, spectrum.lmax, spectrum.mmax, i)
matrix_size(data::AbstractArray, lmax::Integer, mmax::Integer) = (lmax, mmax, size(data)[2:end]...)

# extend to get the size of the i-th dimension, with 1 returned for any additional dimension
# as it is also defined for Array
function matrix_size(data::AbstractArray, lmax::Integer, mmax::Integer, i::Integer)
    i == 1 && return lmax      # first dimension is the number of rows m
    i == 2 && return mmax      # second dimension is the number of columns n
    return size(data, i-1)  # -1 as lmax, mmax are collapsed into a vector in the data array
end 

"""$(TYPEDSIGNATURES)
Size of a `LowerTriangularArray` defined as size of the flattened array if `as <: AbstractVector`
and as if it were a full matrix when `as <: AbstractMatrix`` ."""
Base.size(L::LowerTriangularArray, base::Type{<:IndexBasis}=OneBased; as=Vector) = size(L, base, as)
Base.size(L::LowerTriangularArray, i::Integer, base::Type{<:IndexBasis}=OneBased; as=Vector) = size(L, i, base, as)

# use multiple dispatch to chose the right options of basis and vector/flat vs matrix indexing
# matrix indexing can be zero based (natural for spherical harmonics) or one-based,
# vector/flat indexing has only one-based indexing
Base.size(L::LowerTriangularArray, base::Type{OneBased}, as::Type{Matrix}) = matrix_size(L.data, L.spectrum.lmax, L.spectrum.mmax)
Base.size(L::LowerTriangularArray, base::Type{ZeroBased}, as::Type{Matrix}) = matrix_size(L.data, L.spectrum.lmax-1, L.spectrum.mmax-1)
Base.size(L::LowerTriangularArray, base::Type{OneBased}, as::Type{Vector}) = size(L.data)

# size(L, i, ...) to get the size of the i-th dimension, with 1 returned for any additional dimension 
Base.size(L::LowerTriangularArray, i::Integer, base::Type{OneBased}, as::Type{Matrix}) = matrix_size(L.data, L.spectrum.lmax, L.spectrum.mmax, i)
Base.size(L::LowerTriangularArray, i::Integer, base::Type{ZeroBased}, as::Type{Matrix}) = matrix_size(L.data, L.spectrum.lmax-1, L.spectrum.mmax-1, i)
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

function Base.show(io::IO, ::MIME"text/plain", L::LowerTriangularArray)
    Base.array_summary(io, L, axes(L))

    if get(io, :limit, false)::Bool && displaysize(io)[1]-4 <= 0
        return print(io, " …")
    else
        println(io)
    end
    
    Base.print_array(io, L.data)
end

@inline Base.dataids(L::LowerTriangularArray) = Base.dataids(L.data)

# CREATE INSTANCES (ZEROS, ONES, UNDEF)
for f in (:zeros, :ones, :rand, :randn)
    @eval begin
        # use ArrayType from LowerTriangularArray parameter
        function Base.$f(
            ::Type{LowerTriangularArray{T, N, ArrayType, S}}, 
            lmax::Integer,
            mmax::Integer,
            I::Vararg{Integer, M},
        ) where {T, N, M, ArrayType, S<:AbstractSpectrum}
            ArrayType_ = nonparametric_type(ArrayType)
            return LowerTriangularArray(ArrayType_($f(T, nonzeros(lmax, mmax), I...)), Spectrum(lmax, mmax, architecture=architecture(ArrayType_)))
        end

        function Base.$f(
            ::Type{LowerTriangularArray{T, N, ArrayType, S}}, 
            spectrum::AbstractSpectrum,
            I::Vararg{Integer, M},
        ) where {T, N, M, ArrayType, S<:AbstractSpectrum}
            ArrayType_ = nonparametric_type(ArrayType)
            return LowerTriangularArray(ArrayType_($f(T, nonzeros(spectrum), I...)), spectrum)
        end
        
        # default CPU, use Array
        function Base.$f(
            ::Type{LowerTriangularArray{T}},
            lmax::Integer,
            mmax::Integer, 
            I::Vararg{Integer, M},
        ) where {T, M}
            return LowerTriangularArray($f(T, nonzeros(lmax, mmax), I...), Spectrum(lmax, mmax, architecture=architecture(Array)))
        end

        function Base.$f(
            ::Type{LowerTriangularMatrix{T}},
            lmax::Integer,
            mmax::Integer, 
        ) where T
            return LowerTriangularMatrix($f(T, nonzeros(lmax, mmax)), Spectrum(lmax, mmax, architecture=architecture(Array)))
        end
        
        function Base.$f(
            ::Type{LowerTriangularArray{T}},
            spectrum::AbstractSpectrum,
            I::Vararg{Integer, M},
        ) where {T, M}
            return LowerTriangularArray($f(T, nonzeros(spectrum), I...), spectrum)
        end

        function Base.$f(
            ::Type{LowerTriangularMatrix{T}},
            spectrum::AbstractSpectrum,
        ) where T
            return LowerTriangularMatrix($f(T, nonzeros(spectrum)), spectrum)
        end

        function Base.$f(
            ::Type{T},
            spectrum::AbstractSpectrum,
            I::Vararg{Integer, M},
        ) where {T <: Number, M}
            return LowerTriangularArray($f(T, nonzeros(spectrum), I...), spectrum)
        end

        # use Float64 as default if type T is not provided
        Base.$f(::Type{LowerTriangularArray}, lmax::Integer, mk::Integer...) =
            $f(LowerTriangularArray{Float64}, lmax, mk...)
        Base.$f(::Type{LowerTriangularMatrix}, lmax::Integer, mmax::Integer) =
            $f(LowerTriangularArray{Float64}, lmax, mmax)

        Base.$f(::Type{LowerTriangularArray}, spectrum::AbstractSpectrum, I::Vararg{Integer, M}) where {M} =
            $f(LowerTriangularArray{Float64}, spectrum, I...)
        Base.$f(::Type{LowerTriangularMatrix}, spectrum::AbstractSpectrum) =
            $f(LowerTriangularArray{Float64}, spectrum)
        Base.$f(spectrum::AbstractSpectrum, I::Vararg{Integer, M}) where M =
            $f(LowerTriangularArray{Float64}, spectrum, I...)
    end
end

Base.zero(L::LTA) where {LTA <: LowerTriangularArray} = zeros(LTA, L.spectrum, size(L)[2:end]...)
Base.one(L::LTA) where {LTA <: LowerTriangularArray} = ones(LTA, L.spectrum, size(L)[2:end]...)

function LowerTriangularArray{T, N, ArrayType, S}(
    ::UndefInitializer,
    I::Vararg{Integer,M},
) where {T, N, M, ArrayType<:AbstractArray{T}, S<:AbstractSpectrum}
    return LowerTriangularArray(ArrayType(undef, nonzeros(I[1], I[2]), I[3:end]...), Spectrum(I[1], I[2], architecture=architecture(ArrayType)))
end

function LowerTriangularArray{T, N, ArrayType, S}(
    ::UndefInitializer,
    spectrum::AbstractSpectrum,
    I::Vararg{Integer,M},
) where {T, N, M, ArrayType<:AbstractArray{T}, S<:AbstractSpectrum}
    return LowerTriangularArray(ArrayType(undef, nonzeros(spectrum), I...), spectrum)
end

function LowerTriangularArray{T, N, ArrayType, SP}(
    ::UndefInitializer,
    size::S,
) where {T, N, SP<:AbstractSpectrum, ArrayType<:AbstractArray{T}, S<:Tuple}
    return LowerTriangularArray{T, N, ArrayType, SP}(undef, size...)
end

# note the following constructors hardcode Vector TODO generalise ?
function LowerTriangularMatrix{T}(::UndefInitializer, lmax::Integer, mmax::Integer) where T
    return LowerTriangularMatrix(Vector{T}(undef, nonzeros(lmax, mmax)), lmax, mmax)
end

function LowerTriangularMatrix{T}(::UndefInitializer, spectrum::AbstractSpectrum) where T
    return LowerTriangularMatrix(Vector{T}(undef, nonzeros(spectrum)), spectrum)
end

Base.eltype(L::LowerTriangularArray) = eltype(L.data)

# INDEXING
"""
$(TYPEDSIGNATURES)
Converts the index pair `l, m` of an `lmax`x`mmax` LowerTriangularMatrix `L` to a single
index `i` that indexes the same element in the corresponding vector that stores
only the lower triangle (the non-zero entries) of `L`."""
@inline lm2i(l::Integer, m::Integer, lmax::Integer) = l + (m-1)*lmax - m*(m-1)÷2

"""
$(TYPEDSIGNATURES)
range of the running indices lm in a l-column (degrees of spherical harmonics)
given the column index m (order of harmonics) 
"""
get_lm_range(m, lmax) = lm2i(2*m - 1, m, lmax):lm2i(lmax+m, m, lmax)

"""
$(TYPEDSIGNATURES)
range of the doubled running indices 2lm in a l-column (degrees of spherical harmonics)
given the column index m (order of harmonics) 
"""
get_2lm_range(m, lmax) = 2*lm2i(2*m - 1, m, lmax)-1:2*lm2i(lmax+m, m, lmax)
 
"""
$(TYPEDSIGNATURES)
Converts the linear index `i` in the lower triangle into a pair `(l, m)` of indices 
of the matrix in column-major form. (Formula taken from 
Angeletti et al, 2019, https://hal.science/hal-02047514/document)
"""
@inline function i2lm(k::Integer, mmax::Integer) 
    kp = triangle_number(mmax) - k 
    p = Int(floor((sqrt(1 + 8*kp) - 1)/2))
    l = k - mmax*(mmax-1)÷2 + p*(p+1)÷2
    m = mmax - p
    return l, m
end 
i2lm(I::CartesianIndex, mmax::Int) = CartesianIndex(i2lm(I[1], mmax)...,I.I[2:end]...) 

# direct indexing, no. indices have to be equal to `N` for the correct dimensionality
@inline Base.getindex(L::LowerTriangularArray{T, N}, I::Vararg{Any, N}) where {T, N} = getindex(L.data, I...) 

# indexing with : + other indices, returns a LowerTriangularArray
@inline function Base.getindex(L::LowerTriangularArray{T,N}, col::Colon, I...) where {T,N}
    return LowerTriangularArray(getindex(L.data, col, I...), L.spectrum.lmax, L.spectrum.mmax)
end

# l,m sph "matrix-style"  indexing with integer + other indices
@inline function Base.getindex(L::LowerTriangularArray{T,N}, I::Vararg{Any, M}) where {T, N, M}
    l, m = I[1:2]
    @boundscheck (0 < l <= L.spectrum.lmax && 0 < m <= L.spectrum.mmax) || throw(BoundsError(L, (l, m)))
    # to get a zero element in the correct shape, we just take the zero element of some valid element,
    # there are probably faster ways to do this, but I don't know how, and this is just a fallback anyway 
    @boundscheck m > l && return zero(getindex(L.data, 1, I[3:end]...)) 
    k = lm2i(l, m, L.spectrum.lmax)
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
Base.@propagate_inbounds Base.getindex(L::LowerTriangularArray{T,1,V,S}, i::Integer) where {S<:AbstractSpectrum,V<:AbstractVector{T}} where T = getindex(L.data, i)
Base.@propagate_inbounds Base.getindex(L::LowerTriangularArray{T,1,V,S}, I::CartesianIndex{M}) where {T,S<:AbstractSpectrum,V<:AbstractVector{T},M} = getindex(L, Tuple(I)...)
Base.@propagate_inbounds Base.getindex(L::LowerTriangularArray{T,1,V,S}, i::Integer, I::CartesianIndex{0}) where {T,S<:AbstractSpectrum,V<:AbstractVector{T}} = getindex(L, i)
Base.@propagate_inbounds Base.getindex(L::LowerTriangularArray{T,1,V,S}, i::Integer, I::CartesianIndices{0}) where {T,S<:AbstractSpectrum,V<:AbstractVector{T}} = getindex(L, i)

# setindex with lm, ..
@inline Base.setindex!(L::LowerTriangularArray{T,N}, x, I::Vararg{Any, N}) where {T, N} = setindex!(L.data, x, I...)

# setindex with il, im, ..
@inline function Base.setindex!(L::LowerTriangularArray{T,N}, x, I::Vararg{Any, M}) where {T, N, M}
    @boundscheck N+1==M || throw(BoundsError(L, I))
    i, j = I[1:2] 
    @boundscheck i >= j || throw(BoundsError(L, I))
    k = lm2i(i, j, L.spectrum.lmax)
    setindex!(L.data, x, k, I[3:end]...)
end

# this is specifically here for the weird way in which AssociatedLegendrePolynomials.jl indexes with an empty CartesianIndex
@inline Base.setindex!(L::LowerTriangularArray{T,N}, x, i::CartesianIndex{0}, I::Vararg{Any, M}) where {T,N, M} = setindex!(L, x, I)

@inline Base.setindex!(L::LowerTriangularArray, x, I::CartesianIndex) = setindex!(L, x, Tuple(I)...)
@inline Base.setindex!(L::LowerTriangularArray{T,N}, x, I::CartesianIndex{N}) where {T,N} = setindex!(L.data, x, I)
@inline Base.setindex!(L::LowerTriangularArray{T,N}, x, i::Integer) where {T,N} = setindex!(L.data, x, i)

# this is to avoid ambigouities
@inline Base.setindex!(L::LowerTriangularArray{T,1,V,S}, x, i::Integer) where {T,S<:AbstractSpectrum,V<:AbstractVector{T}} = setindex!(L.data, x, i)
@inline Base.setindex!(L::LowerTriangularArray{T,1,V,S}, x, I::CartesianIndex{1}) where {T,S<:AbstractSpectrum,V<:AbstractVector{T}} = setindex!(L.data, x, I)
@inline Base.setindex!(L::LowerTriangularArray{T,1,V,S}, x, I::CartesianIndex{2}) where {T,S<:AbstractSpectrum,V<:AbstractVector{T}} = setindex!(L, x, Tuple(I)...)
@inline Base.setindex!(L::LowerTriangularArray{T,1,V,S}, x, i::Integer, I::CartesianIndex{0}) where {T,S<:AbstractSpectrum,V<:AbstractVector{T}} = setindex!(L, x, i)

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
creates `unit_range::UnitRange` to loop over all non-zeros in the LowerTriangularArrays
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

"""$(TYPEDSIGNATURES) Iterator for the order m, for each m return all ls, 
therefore the columns in the lower triangular matrix.

    for lms in eachorder(L)
        for lm in lms
            L[lm] 
        end
    end
    
to loop over every order of L."""
eachorder(L1::LowerTriangularArray) = eachorder(L1.spectrum)

function eachorder(L1::LowerTriangularArray, Ls::LowerTriangularArray...)
    lowertriangular_match(L1, Ls...) || throw(DimensionMismatch(L1, Ls...))
    return eachorder(L1)
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
function LowerTriangularMatrix(M::Matrix{T}, spectrum::AbstractSpectrum) where T # CPU version 
    (; lmax, mmax) = spectrum
    L = LowerTriangularMatrix{T}(undef, spectrum)
    
    k = 0
    @inbounds for j in 1:mmax      # only loop over lower triangle
        for i in j:lmax
            k += 1              # next element in lower triangle
            L[k] = M[i, j]       # and copy data into vector
        end
    end
    return L
end

""" 
$(TYPEDSIGNATURES)
Create a LowerTriangularArray `L` from Matrix `M` by copying over the non-zero elements in `M`."""
function LowerTriangularMatrix(M::Matrix{T}) where T # GPU version 
    lmax, mmax = size(M)
    spectrum = Spectrum(lmax, mmax)
    return LowerTriangularMatrix(M, spectrum)
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
    lmax, mmax = size(M)[1:2]
    LowerTriangularArray(M[lowertriangle_indices(lmax,mmax),[Colon() for i=1:N-2]...], lmax, mmax)
end

LowerTriangularArray(L::LowerTriangularArray) = LowerTriangularArray(L.data, L.spectrum)
LowerTriangularMatrix(M::AbstractMatrix) = LowerTriangularArray(M)

function Base.Matrix(L::LowerTriangularMatrix{T}) where T
    lmax, mmax = size(L, as=Matrix)
    M = zeros(T, lmax, mmax)
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
    L1::LowerTriangularArray{T,N,ArrayType1,SP1},   # copy to L1
    L2::LowerTriangularArray{S,N,ArrayType2,SP2},   # copy from L2
    ls::AbstractUnitRange,                      # range of indices in 1st dim
    ms::AbstractUnitRange,                      # range of indices in 2nd dim
) where {T, S, N, SP1 <: AbstractSpectrum, SP2 <: AbstractSpectrum, ArrayType1<:Array{T,N}, ArrayType2<:Array{S,N}}

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
    L1::LowerTriangularArray{T, N, ArrayType1, SP1},
    L2::LowerTriangularArray{S, N, ArrayType2, SP2},
    ls::AbstractUnitRange,
    ms::AbstractUnitRange
) where {T, S, N, SP1 <: AbstractSpectrum, SP2 <: AbstractSpectrum, ArrayType1<:AbstractArray{T}, ArrayType2<:AbstractArray{S}}
    return _copyto_core!(L1, L2, ls, ms)
end

function _copyto_core!( 
    L1::LowerTriangularArray{T,N,ArrayType1,SP1},   # copy to L1
    L2::LowerTriangularArray{S,N,ArrayType2,SP2},   # copy from L2
    ls::AbstractUnitRange,                      # range of indices in 1st dim
    ms::AbstractUnitRange,                      # range of indices in 2nd dim
) where {T, S, N, SP1 <: AbstractSpectrum, SP2 <: AbstractSpectrum, ArrayType1<:AbstractArray{T}, ArrayType2<:AbstractArray{S}}

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
    ::Type{LowerTriangularArray{T1, N, ArrayTypeT1, S}},
    L::LowerTriangularArray{T2, N, ArrayTypeT2, S},
) where {T1, T2, N, S <: AbstractSpectrum, ArrayTypeT1<:AbstractArray{T1}, ArrayTypeT2<:AbstractArray{T2}}
    return LowerTriangularArray{T1, N, ArrayTypeT1, S}(L.data, L.spectrum)
end

function Base.convert(::Type{LowerTriangularMatrix{T}}, L::LowerTriangularMatrix) where T
    return LowerTriangularArray(T.(L.data), L.spectrum)
end

function Base.similar(L::LowerTriangularArray{T, N, ArrayType, SP}, I::Integer...) where {T, N, ArrayType, SP}
    if resolution(L.spectrum) != I[1:2]
        new_spectrum = Spectrum(I[1], I[2], architecture=L.spectrum.architecture)
        return LowerTriangularArray{T,N,ArrayType,SP}(undef, new_spectrum, I[3:end]...)
    else
        return LowerTriangularArray{T,N,ArrayType,SP}(undef, L.spectrum, I[3:end]...)
    end 
end

function Base.similar(L::LowerTriangularArray{T, N, ArrayType, SP}, size::S) where {T, N, ArrayType, SP, S <: Tuple}
    if resolution(L.spectrum) != size[1:2]
        new_spectrum = Spectrum(size[1], size[2], architecture=L.spectrum.architecture)
        return LowerTriangularArray{T,N,ArrayType,SP}(undef, new_spectrum, size[3:end]...)
    else
        return LowerTriangularArray{T,N,ArrayType,SP}(undef, L.spectrum, size[3:end]...)
    end
end

function Base.similar(L::LowerTriangularArray{S, N, ArrayType, SP}, ::Type{T}) where {T, S, N, SP, ArrayType}
    ArrayType_ = nonparametric_type(ArrayType) 
    return LowerTriangularArray{T, N, ArrayType_{T,N}, SP}(similar(L.data, T), L.spectrum)
end

Base.similar(L::LowerTriangularArray{T, N, ArrayType, SP}, ::Type{T}) where {T, N, ArrayType, SP} =
    LowerTriangularArray{T, N, ArrayType, SP}(similar(L.data, T), L.spectrum)
Base.similar(L::LowerTriangularArray{T}) where T = similar(L, T)
 
array_type(::Type{<:LowerTriangularArray{T, N, ArrayType}}) where {T, N, ArrayType} = ArrayType
array_type(L::LowerTriangularArray) = array_type(typeof(L))

Base.prod(L::LowerTriangularArray{NF}) where NF = zero(NF)
@inline Base.sum(L::LowerTriangularArray; dims=:, kw...) = sum(L.data; dims, kw...)

"""
$(TYPEDSIGNATURES)
Fills the elements of `L` with `x`. Faster than fill!(::AbstractArray, x)
as only the non-zero elements in `L` are assigned with x."""
function Base.fill!(L::LowerTriangularArray, x)
    fill!(L.data, x)
    return L
end

Base.:(==)(L1::LowerTriangularArray, L2::LowerTriangularArray) = 
    L1.spectrum == L2.spectrum && L1.data == L2.data
Base.isapprox(L1::LowerTriangularArray, L2::LowerTriangularArray; kwargs...) =
    isapprox(L1.data, L2.data; kwargs...)
Base.all(L::LowerTriangularArray) = all(L.data)
Base.any(L::LowerTriangularArray) = any(L.data)

Base.repeat(L::LowerTriangularArray, counts...) = LowerTriangularArray(repeat(L.data, counts...), L.spectrum)

# Views that return a LowerTriangularArray again (need to retain all horizontal grid points, hence `:, 1` for example)
# view(array, :) unravels like array[:] does hence "::Colon, i, args..." used to enforce one argument after :
# exception is view(vector, :) which preserves the vector structure, equivalent here is the LowerTriangularMatrix
# TODO extend Base.view?
lta_view(L::LowerTriangularArray,  c::Colon, i, args...) = LowerTriangularArray(view(L.data, c, i, args...), L.spectrum)
lta_view(L::LowerTriangularMatrix, c::Colon) = LowerTriangularArray(view(L.data, c), L.spectrum)
lta_view(L::LowerTriangularArray, args...) = view(L, args...)   # fallback to normal view

# Broadcast CPU/GPU
import Base.Broadcast: BroadcastStyle, Broadcasted, DefaultArrayStyle
import LinearAlgebra: isstructurepreserving, fzeropreserving

# CPU with scalar indexing
struct LowerTriangularStyle{N, ArrayType, S} <: Broadcast.AbstractArrayStyle{N} end

# GPU without scalar indexing
struct LowerTriangularGPUStyle{N, ArrayType, S} <: GPUArrays.AbstractGPUArrayStyle{N} end

function BroadcastStyle(::Type{LowerTriangularArray{T, N, ArrayType, S}}) where {T, N, ArrayType, S}
    # remove number format parameter for broadcasting with type promotion
    ArrayType_ = nonparametric_type(ArrayType)
    return LowerTriangularStyle{N, ArrayType_, S}()
end

function BroadcastStyle(
    ::Type{LowerTriangularArray{T, N, ArrayType, S}},
) where {T, N, ArrayType <: GPUArrays.AbstractGPUArray, S}
    # remove number format parameter for broadcasting with type promotion
    ArrayType_ = nonparametric_type(ArrayType)
    return LowerTriangularGPUStyle{N, ArrayType_, S}()
end

# ::Val{0} for broadcasting with 0-dimensional, ::Val{1} for broadcasting with vectors, etc
LowerTriangularStyle{N, ArrayType, S}(::Val{M}) where {N, ArrayType, S, M} =
    LowerTriangularStyle{N, ArrayType, S}()
LowerTriangularGPUStyle{N, ArrayType, S}(::Val{M}) where {N, ArrayType, S, M} =
    LowerTriangularGPUStyle{N, ArrayType, S}()

# also needed for other array types
nonparametric_type(::Type{<:Array}) = Array

# nonparametric_type for a SubArray is the arraytype it is viewing. Needed to construct new arrays from SubArrays!
nonparametric_type(::Type{<:SubArray{T, N, A}}) where {T, N, A} = nonparametric_type(A)
nonparametric_type(::Type{<:SubArray}) = SubArray   # if ArrayType A is not specified, return SubArray

"`L = find_L(Ls)` returns the first LowerTriangularArray among the arguments. 
Adapted from Julia documentation of Broadcast interface"
find_L(bc::Base.Broadcast.Broadcasted) = find_L(bc.args)
find_L(args::Tuple) = find_L(find_L(args[1]), Base.tail(args))
find_L(x) = x
find_L(::Tuple{}) = nothing
find_L(a::LowerTriangularArray, rest) = a
find_L(::Any, rest) = find_L(rest)

function Base.similar(
    bc::Broadcasted{LowerTriangularStyle{N, ArrayType, S}},
    ::Type{T},
) where {N, ArrayType, S, T}
    L = find_L(bc)
    return LowerTriangularArray{T, N, ArrayType{T,N}, S}(similar(L.data, T), L.spectrum)
end

# same function as above, but needs to be defined for both CPU and GPU style
function Base.similar(
    bc::Broadcasted{LowerTriangularGPUStyle{N, ArrayType, S}},
    ::Type{T},
) where {N, ArrayType, S, T}
    L = find_L(bc)
    return LowerTriangularArray{T, N, ArrayType{T,N}, S}(similar(L.data, T), L.spectrum)
end

function KernelAbstractions.get_backend( 
    a::LowerTriangularArray{T, N, ArrayType, S} 
) where {T, N, ArrayType, S} 
    return KernelAbstractions.get_backend(a.data) 
end 

function Adapt.adapt_structure(to, L::LowerTriangularArray) 
    adapted_data = adapt(to, L.data)
    if ismatching(L.spectrum, typeof(adapted_data)) # if matching, adapt the same spectrum 
        return LowerTriangularArray(adapted_data, adapt(to, L.spectrum))
    else # if not matching, create new spectrum with other architecture
        #@warn "Adapting LowerTriangularArray to new architecture with $(typeof(adapted_data))"
        return LowerTriangularArray(adapted_data, adapt(to, Spectrum(L.spectrum, architecture=architecture(typeof(adapted_data)))))
    end
end

Architectures.architecture(L::LowerTriangularArray) = architecture(L.spectrum)
on_architecture(arch::AbstractArchitecture, a::LowerTriangularArray) = Adapt.adapt(array_type(arch), a)
