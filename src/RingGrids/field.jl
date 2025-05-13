abstract type AbstractField{T, N, ArrayType, Grid} <: AbstractArray{T, N} end

struct Field{T, N, ArrayType <: AbstractArray{T, N}, Grid <: AbstractGrid} <: AbstractField{T, N, ArrayType, Grid}
    data::ArrayType
    grid::Grid
end

# for fields, also add info about the number of rings
function Base.array_summary(io::IO, field::AbstractField, inds::Tuple{Vararg{Base.OneTo}})
    print(io, Base.dims2string(length.(inds)), ", $(get_nlat(field))-ring ")
    Base.showarg(io, field, true)
end

isreduced(::Type{<:AbstractField}) = true
isreduced(::Type{<:AbstractField{T, N, A, <:AbstractFullGrid}}) where {T, N, A} = false
isreduced(field::Field) = isreduced(typeof(field))
isreduced(::Type{<:AbstractGrid}) = true
isreduced(::Type{<:AbstractFullGrid}) = false
isreduced(grid::AbstractGrid) = isreduced(typeof(grid))
isfull(F) = ~isreduced(F)

Base.size(f::Field, args...) = size(f.data, args...)

## CONVERSION
# convert an AbstractMatrix to the full grids, and vice versa
"""
($TYPEDSIGNATURES)
Initialize an instance of the grid from an Array. For keyword argument `input_as=Vector` (default)
the leading dimension is interpreted as a flat vector of all horizontal entries in one layer.
For `input_as==Matrx` the first two leading dimensions are interpreted as longitute and latitude.
This is only possible for full grids that are a subtype of `AbstractFullGridArray`.
"""
(Grid::Type{<:AbstractFullGridArray})(M::AbstractArray; input_as=Vector) = Grid(M, input_as)

function (Grid::Type{<:AbstractFullGridArray})(M::AbstractArray, input_as::Type{Matrix})
    # flatten the two horizontal dimensions into one, identical to vec(M) for M <: AbstractMatrix
    M_flat = reshape(M, :, size(M)[3:end]...)
    Grid(M_flat)
end 

Base.Array(grid::AbstractFullGridArray) = Array(reshape(grid.data, :, get_nlat(grid), size(grid.data)[2:end]...))
Base.Matrix(grid::AbstractFullGridArray) = Array(grid)

## INDEXING
# simply propagate all indices forward
Base.@propagate_inbounds Base.getindex(G::AbstractGridArray, ijk...) = getindex(G.data, ijk...)

@inline function Base.getindex(
    G::GridArray,
    col::Colon,
    k...,
) where {GridArray<:AbstractGridArray}
    GridArray_ = nonparametric_type(GridArray)  # obtain parameters from G.data
    return GridArray_(getindex(G.data, col, k...), G.nlat_half, G.rings)
end

# simply propagate all indices forward
Base.@propagate_inbounds Base.setindex!(G::AbstractGridArray, x, ijk...) = setindex!(G.data, x, ijk...)
Base.fill!(G::AbstractGridArray, x) = fill!(G.data, x)

Base.@propagate_inbounds Base.getindex(f::Field, args...) = getindex(f.data, args...)

@inline eachring(field::AbstractField) = eachring(field.grid)



## BROADCASTING
# following https://docs.julialang.org/en/v1/manual/interfaces/#man-interfaces-broadcasting
import Base.Broadcast: BroadcastStyle, Broadcasted, DefaultArrayStyle

# {1} as grids are <:AbstractVector, Grid here is the non-parameteric Grid type!
struct AbstractFieldStyle{N, Field} <: Broadcast.AbstractArrayStyle{N} end

# important to remove Field{T} parameter T (eltype/number format) here to broadcast
# automatically across the same grid type but with different T
# e.g. FullGaussianGrid{Float32} and FullGaussianGrid{Float64}
Base.BroadcastStyle(::Type{Grid}) where {Grid<:AbstractGridArray{T, N, ArrayType}} where {T, N, ArrayType} =
    AbstractGridArrayStyle{N, nonparametric_type(Grid)}()

# allocation for broadcasting, create a new Grid with undef of type/number format T
function Base.similar(bc::Broadcasted{AbstractGridArrayStyle{N, Grid}}, ::Type{T}) where {N, Grid, T}
    return Grid(Array{T}(undef, size(bc)))
end

# ::Val{0} for broadcasting with 0-dimensional, ::Val{1} for broadcasting with vectors, etc
# when there's a dimension mismatch always choose the larger dimension
AbstractGridArrayStyle{N, Grid}(::Val{N}) where {N, Grid} = AbstractGridArrayStyle{N, Grid}()
AbstractGridArrayStyle{1, Grid}(::Val{2}) where {Grid} = AbstractGridArrayStyle{2, Grid}()
AbstractGridArrayStyle{1, Grid}(::Val{0}) where {Grid} = AbstractGridArrayStyle{1, Grid}()
AbstractGridArrayStyle{2, Grid}(::Val{3}) where {Grid} = AbstractGridArrayStyle{3, Grid}()
AbstractGridArrayStyle{2, Grid}(::Val{1}) where {Grid} = AbstractGridArrayStyle{2, Grid}()
AbstractGridArrayStyle{3, Grid}(::Val{4}) where {Grid} = AbstractGridArrayStyle{4, Grid}()
AbstractGridArrayStyle{3, Grid}(::Val{2}) where {Grid} = AbstractGridArrayStyle{3, Grid}()

## GPU
struct AbstractGPUFieldStyle{N, ArrayType, Field} <: GPUArrays.AbstractGPUArrayStyle{N} end

function Base.BroadcastStyle(
    ::Type{Field}
) where {Field<:AbstractField{T, N, ArrayType}} where {T, N, ArrayType <: GPUArrays.AbstractGPUArray}
    return AbstractGPUFieldStyle{N, ArrayType, nonparametric_type(Field)}()
end

# ::Val{0} for broadcasting with 0-dimensional, ::Val{1} for broadcasting with vectors, etc
# when there's a dimension mismatch always choose the larger dimension
AbstractGPUGridArrayStyle{N, ArrayType, Grid}(::Val{N}) where {N, ArrayType, Grid} =
    AbstractGPUGridArrayStyle{N, ArrayType, Grid}()

AbstractGPUFieldStyle{1, ArrayType, Field}(::Val{2}) where {ArrayType, Field} = AbstractGPUGridArrayStyle{2, ArrayType, Grid}()
AbstractGPUFieldStyle{1, ArrayType, Field}(::Val{0}) where {ArrayType, Field} = AbstractGPUGridArrayStyle{1, ArrayType, Grid}()
AbstractGPUFieldStyle{2, ArrayType, Field}(::Val{3}) where {ArrayType, Field} = AbstractGPUGridArrayStyle{3, ArrayType, Grid}()
AbstractGPUFieldStyle{2, ArrayType, Field}(::Val{1}) where {ArrayType, Field} = AbstractGPUGridArrayStyle{2, ArrayType, Grid}()
AbstractGPUFieldStyle{3, ArrayType, Field}(::Val{4}) where {ArrayType, Field} = AbstractGPUGridArrayStyle{4, ArrayType, Grid}()
AbstractGPUFieldStyle{3, ArrayType, Field}(::Val{2}) where {ArrayType, Field} = AbstractGPUGridArrayStyle{3, ArrayType, Grid}()

function KernelAbstractions.get_backend(
    g::Grid
) where {Grid <: AbstractGridArray{T, N, ArrayType}} where {T, N, ArrayType <: GPUArrays.AbstractGPUArray}
    return KernelAbstractions.get_backend(g.data)
end

function Base.similar(
    bc::Broadcasted{AbstractGPUGridArrayStyle{N, ArrayType, Grid}},
    ::Type{T},
) where {N, ArrayType, Grid, T}
    ArrayType_ = nonparametric_type(ArrayType)
    return Grid(ArrayType_{T}(undef, size(bc)))
end

function Adapt.adapt_structure(to, grid::Grid) where {Grid <: AbstractGridArray}
    Grid_ = nonparametric_type(Grid)
    return Grid_(Adapt.adapt(to, grid.data), grid.nlat_half, grid.rings)
end