struct Field{T, N, ArrayType <: AbstractArray{T, N}, Grid <: AbstractGrid} <: AbstractField{T, N, ArrayType, Grid}
    data::ArrayType
    grid::Grid
end

# for fields, also add info about the number of rings
function Base.array_summary(io::IO, field::AbstractField, inds::Tuple{Vararg{Base.OneTo}})
    print(io, Base.dims2string(length.(inds)), ", $(get_nlat(field))-ring ")
    Base.showarg(io, field, true)
end

# TYPE: reduced (fewer points around poles) or full (constant number of longitudes)
isreduced(::Type{<:AbstractField}) = true
isreduced(::Type{<:AbstractField{T, N, A, <:AbstractFullGrid}}) where {T, N, A} = false
isreduced(field::AbstractField) = isreduced(typeof(field))

# SIZE
Base.size(f::AbstractField, args...) = size(f.data, args...)

for f in (:get_nlat, :get_nlat_half, :get_npoints)
    @eval begin
        $f(field::AbstractField) = $f(field.grid)
    end
end

# needed for unalias
@inline Base.dataids(field::AbstractField) = Base.dataids(field.data)

## INDEXING
# simply propagate all indices forward
Base.@propagate_inbounds Base.getindex(field::AbstractField, ijk...) = getindex(field.data, ijk...)
Base.@propagate_inbounds Base.setindex!(field::AbstractField, x, ijk...) = setindex!(field.data, x, ijk...)
Base.fill!(field::AbstractField, x) = fill!(field.data, x)

# @inline function Base.getindex(
#     G::GridArray,
#     col::Colon,
#     k...,
# ) where {GridArray<:AbstractGridArray}
#     GridArray_ = nonparametric_type(GridArray)  # obtain parameters from G.data
#     return GridArray_(getindex(G.data, col, k...), G.nlat_half, G.rings)
# end

# # simply propagate all indices forward

# ## CONSTRUCTORS
# # if no nlat_half provided calculate it
# (::Type{Grid})(M::AbstractArray; input_as=Vector) where Grid<:AbstractGridArray = Grid(M, input_as)

# function (::Type{Grid})(data::AbstractArray, input_as::Type{Vector}) where Grid<:AbstractGridArray
#     npoints2D = size(data, 1)                   # from 1st dim of data
#     nlat_half = get_nlat_half(Grid, npoints2D)  # get nlat_half of Grid
#     return Grid(data, nlat_half)
# end

# function (::Type{Grid})(data::AbstractArray, input_as::Type{Matrix})  where Grid<:AbstractGridArray
#     error("Only full grids can be created from matrix input")
# end


# """$(TYPEDSIGNATURES) True for `data`, `nlat_half` and `rings` that all match in size
# to construct a grid of type `Grid`."""
# function check_inputs(data, nlat_half, rings, Grid)
#     check = true
#     check &= size(data, 1) == get_npoints2D(Grid, nlat_half)    # test number of 2D grid points
#     check &= length(rings) == get_nlat(Grid, nlat_half)         # test number of rings == nlat
#     # TODO also check that rings map to all and only valid grid points?
#     return check
# end

# function error_message(data, nlat_half, rings, G, T, N, A)
#     nlat = get_nlat(G, nlat_half)
#     nrings = length(rings)
#     if nlat != nrings
#         return error("$nrings-element ring indices "*
#             "cannot be used to create a $nlat-ring $G{$T, $N, $A}.")
#     else
#         return error("$(summary(data)) cannot be used to create a $nlat-ring $G{$T, $N, $A}")
#     end
# end

# ## CONVERSION
# # convert an AbstractMatrix to the full grids, and vice versa
# """
# ($TYPEDSIGNATURES)
# Initialize an instance of the grid from an Array. For keyword argument `input_as=Vector` (default)
# the leading dimension is interpreted as a flat vector of all horizontal entries in one layer.
# For `input_as==Matrx` the first two leading dimensions are interpreted as longitute and latitude.
# This is only possible for full grids that are a subtype of `AbstractFullGridArray`.
# """
# (Grid::Type{<:AbstractFullGridArray})(M::AbstractArray; input_as=Vector) = Grid(M, input_as)

# function (Grid::Type{<:AbstractFullGridArray})(M::AbstractArray, input_as::Type{Matrix})
#     # flatten the two horizontal dimensions into one, identical to vec(M) for M <: AbstractMatrix
#     M_flat = reshape(M, :, size(M)[3:end]...)
#     Grid(M_flat)
# end 

# Base.Array(grid::AbstractFullGridArray) = Array(reshape(grid.data, :, get_nlat(grid), size(grid.data)[2:end]...))
# Base.Matrix(grid::AbstractFullGridArray) = Array(grid)

# ## INDEXING
# # simply propagate all indices forward
# Base.@propagate_inbounds Base.getindex(G::AbstractGridArray, ijk...) = getindex(G.data, ijk...)

# @inline function Base.getindex(
#     G::GridArray,
#     col::Colon,
#     k...,
# ) where {GridArray<:AbstractGridArray}
#     GridArray_ = nonparametric_type(GridArray)  # obtain parameters from G.data
#     return GridArray_(getindex(G.data, col, k...), G.nlat_half, G.rings)
# end

# # simply propagate all indices forward
# Base.@propagate_inbounds Base.setindex!(G::AbstractGridArray, x, ijk...) = setindex!(G.data, x, ijk...)
# Base.fill!(G::AbstractGridArray, x) = fill!(G.data, x)

# Base.@propagate_inbounds Base.getindex(f::Field, args...) = getindex(f.data, args...)

# @inline eachring(field::AbstractField) = eachring(field.grid)



# for f in (:zeros, :ones, :rand, :randn)
#     @eval begin
#         # general version with ArrayType(zeros(...)) conversion
#         function Base.$f(
#             ::Type{Grid},
#             nlat_half::Integer,
#             k::Integer...,
#         ) where {Grid<:AbstractGridArray{T, N, ArrayType}} where {T, N, ArrayType}
#             return Grid(ArrayType($f(T, get_npoints2D(Grid, nlat_half), k...)), nlat_half)
#         end

#         # CPU version with zeros(T, ...) producing Array
#         function Base.$f(
#             ::Type{Grid},
#             nlat_half::Integer,
#             k::Integer...,
#         ) where {Grid<:AbstractGridArray{T}} where T
#             return Grid($f(T, get_npoints2D(Grid, nlat_half), k...), nlat_half)
#         end

#         # use Float64 if no type provided
#         function Base.$f(
#             ::Type{Grid},
#             nlat_half::Integer,
#             k::Integer...,
#         ) where {Grid<:AbstractGridArray}
#             return $f(Grid{Float64}, nlat_half, k...)
#         end
#     end
# end

# # zero element of an AbstractGridArray instance grid by creating new zero(grid.data)
# Base.zero(grid::Grid) where {Grid<:AbstractGridArray} =
#     nonparametric_type(Grid)(zero(grid.data), grid.nlat_half, grid.rings)

# # similar data but everything else identical
# function Base.similar(grid::Grid) where {Grid<:AbstractGridArray}
#     return nonparametric_type(Grid)(similar(grid.data), grid.nlat_half, grid.rings)
# end

# # data with new type T but everything else identical
# function Base.similar(grid::Grid, ::Type{T}) where {Grid<:AbstractGridArray, T}
#     return nonparametric_type(Grid)(similar(grid.data, T), grid.nlat_half, grid.rings)
# end

# # data with same type T but new size
# function Base.similar(
#     grid::Grid,
#     nlat_half::Integer,
#     k::Integer...
# ) where {Grid<:AbstractGridArray{T, N, ArrayType}} where {T, N, ArrayType}
#     similar_data = similar(grid.data, get_npoints2D(Grid, nlat_half), k...)
#     return nonparametric_type(Grid)(similar_data, nlat_half)
# end

# # data with new type T and new size
# function Base.similar(
#     grid::Grid,
#     ::Type{Tnew},
#     nlat_half::Integer,
#     k::Integer...
# ) where {Grid<:AbstractGridArray, Tnew}
#     similar_data = similar(grid.data, Tnew, get_npoints2D(Grid, nlat_half), k...)
#     return nonparametric_type(Grid)(similar_data, nlat_half)
# end

# # general version with ArrayType{T, N}(undef, ...) generator
# function (::Type{Grid})(
#     ::UndefInitializer,
#     nlat_half::Integer,
#     k::Integer...,
# ) where {Grid<:AbstractGridArray{T, N, ArrayType}} where {T, N, ArrayType}
#     ArrayType_ = nonparametric_type(ArrayType)
#     return Grid(ArrayType_{T, N}(undef, get_npoints2D(Grid, nlat_half), k...), nlat_half)
# end

# # CPU version with Array{T, N}(undef, ...) generator
# function (::Type{Grid})(
#     ::UndefInitializer,
#     nlat_half::Integer,
#     k::Integer...,
# ) where {Grid<:AbstractGridArray{T}} where T
#     return Grid(Array{T}(undef, get_npoints2D(Grid, nlat_half), k...), nlat_half)
# end

# # use Float64 if no type provided
# function (::Type{Grid})(
#     ::UndefInitializer,
#     nlat_half::Integer,
#     k::Integer...,
# ) where {Grid<:AbstractGridArray}
#     return Grid(Array{Float64}(undef, get_npoints2D(Grid, nlat_half), k...), nlat_half)
# end

# function Base.convert(
#     ::Type{Grid},
#     grid::AbstractGridArray,
# ) where {Grid<:AbstractGridArray{T, N, ArrayType}} where {T, N, ArrayType}
#     return Grid(ArrayType(grid.data))
# end

# # ITERATORS
# eachring(field::AbstractField) = eachring(field.grid)

# """
# $(TYPEDSIGNATURES)
# CartesianIndices for the 2nd to last dimension of an AbstractField,
# e.g. the vertical layer (or a time dimension, etc). This is
# to be used like

# for k in eachlayer(field)
#     for ring in eachring(field)
#         for ij in ring
#             field[ij, k]"""
# @inline eachlayer(field::AbstractField) = CartesianIndices(size(field)[2:end])

# # several arguments to check for matching grids
# function eachlayer(field1::AbstractGrid, fields::AbstractField...; kwargs...)
#     grids_match(field1, fields...; kwargs...) || throw(DimensionMismatch(field1, fields...))
#     return eachlayer(field1)
# end

# # equality and comparison, somehow needed as not covered by broadcasting
# Base.:(==)(F1::AbstractField, F2::AbstractField) = grids_match(F1, F2) && F1.data == F2.data
# Base.all(F::AbstractField) = all(F.data)
# Base.any(F::AbstractField) = any(F.data)

# """$(TYPEDSIGNATURES) True if both `A` and `B` are of the same nonparametric grid type
# (e.g. OctahedralGaussianArray, regardless type parameter `T` or underyling array type `ArrayType`)
# and of same resolution (`nlat_half`) and total grid points (`length`). Sizes of `(4,)` and `(4,1)`
# would match for example, but `(8,1)` and `(4,2)` would not (`nlat_half` not identical)."""
# function grids_match(
#     A::AbstractGridArray,
#     B::AbstractGridArray;
#     horizontal_only::Bool = false,
#     vertical_only::Bool = false,
# )
#     @assert ~(horizontal_only && vertical_only) "Conflicting options: horizontal_only = $horizontal_ony and vertical_only = $vertical_only"

#     horizontal_match = get_nlat_half(A) == get_nlat_half(B)
#     vertical_match = size(A)[2:end] == size(B)[2:end]
#     type_match = grids_match(typeof(A), typeof(B))

#     if horizontal_only
#         # type also has to match as two different grid types can have the same nlat_half
#         return horizontal_match && type_match
#     elseif vertical_only
#         return vertical_match
#     else
#         return horizontal_match && vertical_match && type_match
#     end
# end

# ## BROADCASTING
# # following https://docs.julialang.org/en/v1/manual/interfaces/#man-interfaces-broadcasting
# import Base.Broadcast: BroadcastStyle, Broadcasted, DefaultArrayStyle

# # {1} as grids are <:AbstractVector, Grid here is the non-parameteric Grid type!
# struct AbstractFieldStyle{N, Field} <: Broadcast.AbstractArrayStyle{N} end

# # important to remove Field{T} parameter T (eltype/number format) here to broadcast
# # automatically across the same grid type but with different T
# # e.g. FullGaussianGrid{Float32} and FullGaussianGrid{Float64}
# Base.BroadcastStyle(::Type{Grid}) where {Grid<:AbstractGridArray{T, N, ArrayType}} where {T, N, ArrayType} =
#     AbstractGridArrayStyle{N, nonparametric_type(Grid)}()

# # allocation for broadcasting, create a new Grid with undef of type/number format T
# function Base.similar(bc::Broadcasted{AbstractGridArrayStyle{N, Grid}}, ::Type{T}) where {N, Grid, T}
#     return Grid(Array{T}(undef, size(bc)))
# end

# # ::Val{0} for broadcasting with 0-dimensional, ::Val{1} for broadcasting with vectors, etc
# # when there's a dimension mismatch always choose the larger dimension
# AbstractGridArrayStyle{N, Grid}(::Val{N}) where {N, Grid} = AbstractGridArrayStyle{N, Grid}()
# AbstractGridArrayStyle{1, Grid}(::Val{2}) where {Grid} = AbstractGridArrayStyle{2, Grid}()
# AbstractGridArrayStyle{1, Grid}(::Val{0}) where {Grid} = AbstractGridArrayStyle{1, Grid}()
# AbstractGridArrayStyle{2, Grid}(::Val{3}) where {Grid} = AbstractGridArrayStyle{3, Grid}()
# AbstractGridArrayStyle{2, Grid}(::Val{1}) where {Grid} = AbstractGridArrayStyle{2, Grid}()
# AbstractGridArrayStyle{3, Grid}(::Val{4}) where {Grid} = AbstractGridArrayStyle{4, Grid}()
# AbstractGridArrayStyle{3, Grid}(::Val{2}) where {Grid} = AbstractGridArrayStyle{3, Grid}()

# ## GPU
# struct AbstractGPUFieldStyle{N, ArrayType, Field} <: GPUArrays.AbstractGPUArrayStyle{N} end

# function Base.BroadcastStyle(
#     ::Type{Field}
# ) where {Field<:AbstractField{T, N, ArrayType}} where {T, N, ArrayType <: GPUArrays.AbstractGPUArray}
#     return AbstractGPUFieldStyle{N, ArrayType, nonparametric_type(Field)}()
# end

# # ::Val{0} for broadcasting with 0-dimensional, ::Val{1} for broadcasting with vectors, etc
# # when there's a dimension mismatch always choose the larger dimension
# AbstractGPUGridArrayStyle{N, ArrayType, Grid}(::Val{N}) where {N, ArrayType, Grid} =
#     AbstractGPUGridArrayStyle{N, ArrayType, Grid}()

# AbstractGPUFieldStyle{1, ArrayType, Field}(::Val{2}) where {ArrayType, Field} = AbstractGPUGridArrayStyle{2, ArrayType, Grid}()
# AbstractGPUFieldStyle{1, ArrayType, Field}(::Val{0}) where {ArrayType, Field} = AbstractGPUGridArrayStyle{1, ArrayType, Grid}()
# AbstractGPUFieldStyle{2, ArrayType, Field}(::Val{3}) where {ArrayType, Field} = AbstractGPUGridArrayStyle{3, ArrayType, Grid}()
# AbstractGPUFieldStyle{2, ArrayType, Field}(::Val{1}) where {ArrayType, Field} = AbstractGPUGridArrayStyle{2, ArrayType, Grid}()
# AbstractGPUFieldStyle{3, ArrayType, Field}(::Val{4}) where {ArrayType, Field} = AbstractGPUGridArrayStyle{4, ArrayType, Grid}()
# AbstractGPUFieldStyle{3, ArrayType, Field}(::Val{2}) where {ArrayType, Field} = AbstractGPUGridArrayStyle{3, ArrayType, Grid}()

# function KernelAbstractions.get_backend(
#     g::Grid
# ) where {Grid <: AbstractGridArray{T, N, ArrayType}} where {T, N, ArrayType <: GPUArrays.AbstractGPUArray}
#     return KernelAbstractions.get_backend(g.data)
# end

# function Base.similar(
#     bc::Broadcasted{AbstractGPUGridArrayStyle{N, ArrayType, Grid}},
#     ::Type{T},
# ) where {N, ArrayType, Grid, T}
#     ArrayType_ = nonparametric_type(ArrayType)
#     return Grid(ArrayType_{T}(undef, size(bc)))
# end

# function Adapt.adapt_structure(to, grid::Grid) where {Grid <: AbstractGridArray}
#     Grid_ = nonparametric_type(Grid)
#     return Grid_(Adapt.adapt(to, grid.data), grid.nlat_half, grid.rings)
# end