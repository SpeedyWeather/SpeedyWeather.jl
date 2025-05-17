const DEFAULT_NF = Float64

struct Field{T, N, ArrayType <: AbstractArray, Grid <: AbstractGrid} <: AbstractField{T, N, ArrayType, Grid}
    data::ArrayType
    grid::Grid

    # Inner constructor to check for matching grid and data
    function Field(data, grid)
        data_matches_grid(data, grid) || throw(DimensionMismatch(data, grid))
        return new{eltype(data), ndims(data), typeof(data), typeof(grid)}(data, grid)
    end
end

# horizontal dimension is always unravelled into a vector
const Field2D = Field{T, 1} where T
const Field3D = Field{T, 2} where T
const Field4D = Field{T, 3} where T

# default constructors
Field(grid::AbstractGrid, k...) = zeros(grid, k...)
Field(::Type{T}, grid::AbstractGrid, k...) where T = zeros(T, grid, k...)

grid_type(field::AbstractField) = grid_type(typeof(field))
array_type(::Type{Field{T, N, A, G}}) where {T, N, A, G} = A

# test number of horizontal grid points matches
data_matches_grid(data, grid) = size(data, 1) == get_npoints(grid)

function Base.DimensionMismatch(data::AbstractArray, grid::AbstractGrid)
    Grid_ = nonparametric_type(grid)
    nlat = get_nlat(grid)
    npoints = get_npoints(grid)
    return DimensionMismatch("$(summary(data)) cannot be used to create a Field with a $npoints-element, $nlat-ring $Grid_")
end

Base.DimensionMismatch(f1::AbstractField, f2::AbstractField) = DimensionMismatch("$(summary(f1)) does not match $(summary(f2))")
function Base.DimensionMismatch(f1::AbstractField, f2s::AbstractField...)
    length(f2s) == 0 && return DimensionMismatch("$(summary(f1)) does not match any other field")
    s = "$(summary(f1)) does not match either of "
    for f2 in f2s
        s *= "$(summary(f2)), "
    end
    DimensionMismatch(chop(s, tail=2))
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
# sizeof: don't count grid as possibly shared with other fields
Base.sizeof(field::AbstractField) = sizeof(field.data)  

for f in (:get_nlat, :get_nlat_half)
    @eval begin
        $f(field::AbstractField) = $f(field.grid)
    end
end

# grid is inherently 2D, for field distinguish between get_npoints (3D+) and get_npoints2D
get_npoints2D(field::AbstractField) = get_npoints(field.grid)
get_npoints(field::AbstractField) = get_npoints(field.grid, size(field)[2:end]...)

# needed for unalias
@inline Base.dataids(field::AbstractField) = Base.dataids(field.data)

## INDEXING
# simply propagate all indices forward
Base.@propagate_inbounds Base.getindex(field::AbstractField, ijk...) = getindex(field.data, ijk...)
Base.@propagate_inbounds Base.setindex!(field::AbstractField, x, ijk...) = setindex!(field.data, x, ijk...)
Base.fill!(field::AbstractField, x) = fill!(field.data, x)

# make [:, k...] not escape the Field
@inline Base.getindex(field::AbstractField, col::Colon, k...) = Field(field.data[col, k...], field.grid)
@inline eachring(field::AbstractField) = eachring(field.grid)

# ITERATORS

"""$(TYPEDSIGNATURES)
CartesianIndices for the 2nd to last dimension of an AbstractField,
e.g. the vertical layer (or a time dimension, etc). To be used like

    for k in eachlayer(field)
        for (j, ring) in enumerate(eachring(field))
            for ij in ring
                field[ij, k]"""
@inline eachlayer(field::AbstractField) = CartesianIndices(size(field)[2:end])

# several arguments to check for matching grids
function eachlayer(field1::AbstractField, fields::AbstractField...; kwargs...)
    fields_match(field1, fields...; kwargs...) || throw(DimensionMismatch(field1, fields...))
    return eachlayer(field1)
end

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

for f in (:zeros, :ones, :rand, :randn)
    @eval begin
        # zeros(grid, nlayers...)
        function Base.$f(grid::AbstractGrid, k...)
            data = array_type(grid.architecture)($f(DEFAULT_NF, get_npoints(grid), k...))
            return Field(data, grid)
        end
            
        # zeros(NF, grid, nlayers...)
        function Base.$f(::Type{T}, grid::AbstractGrid, k...) where T
            data = array_type(grid.architecture)($f(T, get_npoints(grid), k...))
            return Field(data, grid)
        end

        # zeros(Grid, nlat_half, nlayers...)
        # ::Integer added here to avoid ambiguity with array.jl in Base
        function Base.$f(Grid::Type{<:AbstractGrid}, nlat_half::Integer, k::Integer...; architecture=DEFAULT_ARCHITECTURE())
            grid = Grid(nlat_half, architecture)
            data = array_type(architecture)($f(DEFAULT_NF, get_npoints(grid), k...))
            return Field(data, grid)
        end

        # zeros(NF, Grid, nlat_half, nlayers...)
        function Base.$f(
            ::Type{T},
            Grid::Type{<:AbstractGrid},
            nlat_half,
            k...;
            architecture=DEFAULT_ARCHITECTURE(),
        ) where T
            grid = Grid(nlat_half, architecture)
            data = array_type(architecture)($f(T, get_npoints(grid), k...))
            return Field(data, grid)
        end

        function Base.$f(
            ::Type{F},
            nlat_half::Integer,
            k::Integer...;
            architecture=DEFAULT_ARCHITECTURE(),
        ) where {F<:AbstractField}
            Grid_ = grid_type(F)
            grid = nonparametric_type(Grid_)(nlat_half, architecture)
            data = array_type(architecture)($f(get_npoints(grid), k...))
            return Field(data, grid)
        end

        function Base.$f(
            ::Type{F},
            nlat_half::Integer,
            k::Integer...;
            architecture=DEFAULT_ARCHITECTURE(),
        ) where {F<:AbstractField{T}} where T
            Grid_ = grid_type(F)
            grid = nonparametric_type(Grid_)(nlat_half, architecture)
            data = array_type(architecture)($f(T, get_npoints(grid), k...))
            return Field(data, grid)
        end
    end
end

# zero element of a Field with new data but same grid
Base.zero(field::AbstractField) = Field(zero(field.data), field.grid)

# similar data but share grid
Base.similar(field::AbstractField) = Field(similar(field.data), field.grid)

# data with new type T but share grid
Base.similar(field::AbstractField, ::Type{T}) where {T} = Field(similar(field.data, T), field.grid)

# data with same type T but new size (=new grid)
function Base.similar(
    field::AbstractField,
    nlat_half::Integer,
    k::Integer...
)
    # use same architecture though
    new_grid = typeof(field.grid)(nlat_half, field.grid.architecture)
    similar_data = similar(field.data, get_npoints(new_grid), k...)
    return Field(similar_data, new_grid)
end

# data with new type T and new size
function Base.similar(
    field::AbstractField,
    ::Type{Tnew},
    nlat_half::Integer,
    k::Integer...
) where Tnew
    # use same architecture though
    new_grid = typeof(field.grid)(nlat_half, field.grid.architecture)
    similar_data = similar(field.data, Tnew, get_npoints(new_grid), k...)
    return Field(similar_data, new_grid)
end

# general version with ArrayType{T, N}(undef, ...) generator
function (::Type{F})(
    ::UndefInitializer,
    nlat_half::Integer,
    k::Integer...,
) where {F<:AbstractField{T, N, ArrayType, Grid}} where {T, N, ArrayType, Grid<:AbstractGrid{Architecture}} where Architecture
    # this (inevitably) creates a new grid and architecture instance
    grid = nonparametric_type(Grid)(nlat_half, Architecture())
    # TODO this should allocate on the device
    data = nonparametric_type(ArrayType){T}(undef, get_npoints(grid), k...)
    return Field(data, grid)
end

# in case Architecture is not provided use DEFAULT_ARCHITECTURE
function (::Type{F})(
    ::UndefInitializer,
    nlat_half::Integer,
    k::Integer...,
) where {F<:AbstractField{T, N, ArrayType, Grid}} where {T, N, ArrayType, Grid<:AbstractGrid}
    grid = nonparametric_type(Grid)(nlat_half, DEFAULT_ARCHITECTURE())
    data = nonparametric_type(ArrayType){T}(undef, get_npoints(grid), k...)
    return Field(data, grid)
end

# in case only grid is provided (e.g. FullGaussianField) use Float64, Array, DEFAULT_ARCHITECTURE
function (::Type{F})(
    ::UndefInitializer,
    nlat_half::Integer,
    k::Integer...,
) where {F<:AbstractField}
    Grid_ = grid_type(F)
    grid = nonparametric_type(Grid_)(nlat_half, DEFAULT_ARCHITECTURE())
    data = Array{DEFAULT_NF}(undef, get_npoints(grid), k...)
    return Field(data, grid)
end

# in case only number format is provided use Array and DEFAULT_ARCHITECTURE
function (::Type{F})(
    ::UndefInitializer,
    nlat_half::Integer,
    k::Integer...,
) where {F<:AbstractField{T}} where T
    Grid_ = grid_type(F)
    grid = nonparametric_type(Grid_)(nlat_half, DEFAULT_ARCHITECTURE())
    data = Array{T}(undef, get_npoints(grid), k...)
    return Field(data, grid)
end

#  same as above but with N (ignored though as obtained from integer arguments)
function (::Type{F})(
    ::UndefInitializer,
    nlat_half::Integer,
    k::Integer...,
) where {F<:AbstractField{T, N}} where {T, N}
    Grid_ = grid_type(F)
    grid = nonparametric_type(Grid_)(nlat_half, DEFAULT_ARCHITECTURE())
    data = Array{T}(undef, get_npoints(grid), k...)
    return Field(data, grid)
end

# in case ArrayType is provided use that!
function (::Type{F})(
    ::UndefInitializer,
    nlat_half::Integer,
    k::Integer...,
) where {F<:AbstractField{T, N, ArrayType}} where {T, N, ArrayType}
    Grid_ = grid_type(F)
    grid = nonparametric_type(Grid_)(nlat_half, DEFAULT_ARCHITECTURE())
    data = nonparametric_type(ArrayType){T}(undef, get_npoints(grid), k...)
    return Field(data, grid)
end

# TODO is that all that's needed?
function Base.convert(
    ::Type{F},
    field::AbstractField,
) where {F<:AbstractField{T, N, ArrayType, Grid}} where {T, N, ArrayType, Grid}
    return F(field.data, field.grid)
end

# equality and comparison, somehow needed as not covered by broadcasting
Base.:(==)(F1::AbstractField, F2::AbstractField) = fields_match(F1, F2) && F1.data == F2.data
Base.all(F::AbstractField) = all(F.data)
Base.any(F::AbstractField) = any(F.data)

grids_match(A::AbstractField, Bs::AbstractField...) = grids_match(A.grid, (B.grid for B in Bs)...)

"""$(TYPEDSIGNATURES) True if both `A` and `B` are of the same nonparametric grid type
(e.g. OctahedralGaussianArray, regardless type parameter `T` or underyling array type `ArrayType`)
and of same resolution (`nlat_half`) and total grid points (`length`). Sizes of `(4,)` and `(4,1)`
would match for example, but `(8,1)` and `(4,2)` would not (`nlat_half` not identical)."""
function fields_match(
    A::AbstractField,
    B::AbstractField;
    horizontal_only::Bool = false,
    vertical_only::Bool = false,
)
    @assert ~(horizontal_only && vertical_only) "Conflicting options: horizontal_only = $horizontal_ony and vertical_only = $vertical_only"

    horizontal_only && return grids_match(A, B)
    vertical_only && return size(A)[2:end] == size(B)[2:end]
    return grids_match(A, B) && size(A)[2:end] == size(B)[2:end]
end

function fields_match(A::AbstractField, Bs::AbstractField...; kwargs...)
    match = true    # single field A always matches itself
    for B in Bs     # check for all matching respectively with A
        match &= fields_match(A, B; kwargs...)
    end
    return match
end

## BROADCASTING
# following https://docs.julialang.org/en/v1/manual/interfaces/#man-interfaces-broadcasting
import Base.Broadcast: BroadcastStyle, Broadcasted, DefaultArrayStyle

# new broadcasting style add additional parameters of Field
struct FieldStyle{N, ArrayType, Grid} <: Broadcast.AbstractArrayStyle{N} end

# define broadcast style for Field from its parameters
Base.BroadcastStyle(::Type{F}) where {F<:AbstractField{T, N, ArrayType, Grid}} where {T, N, ArrayType, Grid} =
    FieldStyle{N, ArrayType, Grid}()

# allocation for broadcasting via similar, reusing grid from the first field of the broadcast arguments
# e.g. field1 + field2 creates a new field that share the grid of field1
# 2 .+ field1 creates a new field that share the grid of field1
function Base.similar(bc::Broadcasted{FieldStyle{N, ArrayType, Grid}}, ::Type{T}) where {N, ArrayType, Grid, T}
    for maybe_field in bc.args
        if maybe_field isa AbstractField
            return similar(maybe_field, T)
        end
    end
end

# ::Val{0} for broadcasting with 0-dimensional, ::Val{1} for broadcasting with vectors, etc
# when there's a dimension mismatch always choose the larger dimension
FieldStyle{N, ArrayType, Grid}(::Val{N}) where {N, ArrayType, Grid} = FieldStyle{N, ArrayType, Grid}()
FieldStyle{1, ArrayType, Grid}(::Val{2}) where {ArrayType, Grid} = FieldStyle{2, ArrayType, Grid}()
FieldStyle{1, ArrayType, Grid}(::Val{0}) where {ArrayType, Grid} = FieldStyle{1, ArrayType, Grid}()
FieldStyle{2, ArrayType, Grid}(::Val{3}) where {ArrayType, Grid} = FieldStyle{3, ArrayType, Grid}()
FieldStyle{2, ArrayType, Grid}(::Val{1}) where {ArrayType, Grid} = FieldStyle{2, ArrayType, Grid}()
FieldStyle{3, ArrayType, Grid}(::Val{4}) where {ArrayType, Grid} = FieldStyle{4, ArrayType, Grid}()
FieldStyle{3, ArrayType, Grid}(::Val{2}) where {ArrayType, Grid} = FieldStyle{3, ArrayType, Grid}()

## GPU (same but <: GPUArrays.AbstractGPUArrayStyle)
struct FieldGPUStyle{N, ArrayType, Grid} <: GPUArrays.AbstractGPUArrayStyle{N} end

# same as for FieldStyle but for constrain to ArrayType<:GPUArrays
function Base.BroadcastStyle(
    ::Type{F}
) where {F<:AbstractField{T, N, ArrayType, Grid}} where {T, N, ArrayType <: GPUArrays.AbstractGPUArray, Grid}
    return FieldGPUStyle{N, ArrayType, Grid}()
end

function Base.similar(bc::Broadcasted{FieldGPUStyle{N, ArrayType, Grid}}, ::Type{T}) where {N, ArrayType, Grid, T}
    for maybe_field in bc.args
        if maybe_field isa AbstractField
            return similar(maybe_field, T)
        end
    end
end


# ::Val{0} for broadcasting with 0-dimensional, ::Val{1} for broadcasting with vectors, etc
# when there's a dimension mismatch always choose the larger dimension
FieldGPUStyle{N, ArrayType, Grid}(::Val{N}) where {N, ArrayType, Grid} = FieldGPUStyle{N, ArrayType, Grid}()
FieldGPUStyle{1, ArrayType, Grid}(::Val{2}) where {ArrayType, Grid}    = FieldGPUStyle{2, ArrayType, Grid}()
FieldGPUStyle{1, ArrayType, Grid}(::Val{0}) where {ArrayType, Grid}    = FieldGPUStyle{1, ArrayType, Grid}()
FieldGPUStyle{2, ArrayType, Grid}(::Val{3}) where {ArrayType, Grid}    = FieldGPUStyle{3, ArrayType, Grid}()
FieldGPUStyle{2, ArrayType, Grid}(::Val{1}) where {ArrayType, Grid}    = FieldGPUStyle{2, ArrayType, Grid}()
FieldGPUStyle{3, ArrayType, Grid}(::Val{4}) where {ArrayType, Grid}    = FieldGPUStyle{4, ArrayType, Grid}()
FieldGPUStyle{3, ArrayType, Grid}(::Val{2}) where {ArrayType, Grid}    = FieldGPUStyle{3, ArrayType, Grid}()

function KernelAbstractions.get_backend(
    field::F
) where {F <: AbstractField{T, N, ArrayType}} where {T, N, ArrayType <: GPUArrays.AbstractGPUArray}
    return KernelAbstractions.get_backend(field.data)
end

function Adapt.adapt_structure(to, field::F) where {F <: AbstractField}
    # TODO this reuses the same grid but adapt can change the array type which should also change the architecture?
    return Field(Adapt.adapt(to, field.data), field.grid)
end