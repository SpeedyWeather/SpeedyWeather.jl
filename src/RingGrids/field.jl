"""Field contains data on a grid. There is only one `Field` type which can be used for all grids.
`data` is the data array, the first dimension is always the horizontal dimension (unravelled into a vector
for compatibility across full and reduced grids), the other dimensions can be used for the vertical and/or
time or other dimensions. The `grid` can be shared across multiple fields, e.g. a 2D and a 3D field
can share the same grid which just defines the discretization and the architecture (CPU/GPU) the grid is on.
$(TYPEDFIELDS)"""
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
(::Type{<:Field{T}})(data::AbstractArray, grid::AbstractGrid) where T = Field(T.(data), grid)

# TYPES
nonparametric_type(::Type{<:Field}) = Field
grid_type(field::AbstractField) = grid_type(typeof(field))
grid_type(::Type{Field{T, N, A, G}}) where {T, N, A, G} = G
field_type(field::AbstractField) = typeof(field)
field_type(::Type{F}) where {F<:AbstractField} = F
field_type(grid::AbstractGrid) = field_type(typeof(grid))
field_type(::Type{G}) where {G<:AbstractGrid} = Field{T, N, A, G} where {T, N, A}
full_grid_type(field::AbstractField) = full_grid_type(typeof(field.grid))
full_grid_type(::Type{F}) where {F<:AbstractField} = full_grid_type(grid_type(F))
Architectures.array_type(::Type{Field{T, N, A, G}}) where {T, N, A, G} = A
Architectures.array_type(field::AbstractField) = array_type(typeof(field))

# test number of horizontal grid points matches
data_matches_grid(data::AbstractArray, grid::AbstractGrid) = size(data, 1) == get_npoints(grid)

function Base.DimensionMismatch(data::AbstractArray, grid::AbstractGrid)
    Grid_ = nonparametric_type(grid)
    nlat = get_nlat(grid)
    npoints = get_npoints(grid)
    return DimensionMismatch("$(summary(data)) cannot be used to create a Field with a $npoints-element, $nlat-ring $Grid_")
end

Base.DimensionMismatch(f1::AbstractField, f2::AbstractField) = DimensionMismatch("$(summary(f1)) does not match $(summary(f2))")
Base.DimensionMismatch(f1::AbstractField, f2::AbstractArray) = DimensionMismatch("$(summary(f1)) does not match $(summary(f2))")
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

function Base.show(io::IO, ::MIME"text/plain", field::AbstractField)
    Base.array_summary(io, field, axes(field))

    if get(io, :limit, false)::Bool && displaysize(io)[1]-4 <= 0
        return print(io, " …")
    else
        println(io)
    end
    
    Base.print_array(io, field.data)
end

# TYPE: reduced (fewer points around poles) or full (constant number of longitudes)
isreduced(::Type{<:AbstractField}) = true
isreduced(::Type{<:AbstractFullField}) = false
isreduced(field::AbstractField) = isreduced(typeof(field))

# SIZE
Base.size(f::AbstractField, args...) = size(f.data, args...)
# sizeof: don't count grid as possibly shared with other fields
Base.sizeof(field::AbstractField) = sizeof(field.data)  

for f in (  :get_nlat, :get_nlat_half,
            :get_colat, :get_lat, :get_latd, :get_lon, :get_lond,
            :get_londlatds, :get_lonlats, :get_loncolats)
    @eval begin
        $f(field::AbstractField) = $f(field.grid)
    end
end

# grid is inherently 2D, for field distinguish between get_npoints (3D+) and get_npoints2D
get_npoints2D(field::AbstractField) = get_npoints(field.grid)
get_npoints(field::AbstractField) = get_npoints(field.grid, size(field)[2:end]...)
matrix_size(field::AbstractField) = matrix_size(field.grid)

# needed for unalias
@inline Base.dataids(field::AbstractField) = Base.dataids(field.data)

## INDEXING
# simply propagate all indices forward
Base.@propagate_inbounds Base.getindex(field::AbstractField, ijk...) = getindex(field.data, ijk...)
Base.@propagate_inbounds Base.setindex!(field::AbstractField, x, ijk...) = setindex!(field.data, x, ijk...)
Base.fill!(field::AbstractField, x) = fill!(field.data, x)

# make [:, k...] not escape the Field
@inline Base.getindex(field::AbstractField, col::Colon, k...) = Field(field.data[col, k...], field.grid)

# ITERATORS

"""$(TYPEDSIGNATURES)
Return an iterator over the rings of a field. These are precomputed ranges like [1:20, 21:44, ...] stored
in `field.grid.rings`. Every element are the unravellend indices `ij` of all 2D grid points that lie on that ring.
Intended use is like

    for (j, ring) in enumerate(eachring(field))
        for ij in ring
            field[ij, k]
            
with `j` being the ring index. j=1 around the north pole, j=nlat_half around the Equator or just north of it (for even nlat_half).
Several fields can be passed on to check for matching grids, e.g. for a 2D and a 3D field.
If `horizontal_only` is set to `true`, the function will only check the horizontal dimension. """
function eachring(field1::AbstractField, fields::AbstractField...; horizontal_only=true, kwargs...)
    fields_match(field1, fields...; horizontal_only, kwargs...) || throw(DimensionMismatch(field1, fields...))
    return eachring(field1)
end

eachring(field::AbstractField) = eachring(field.grid)

"""$(TYPEDSIGNATURES)
CartesianIndices for the 2nd to last dimension of a field (or fields),
e.g. the vertical layer (and possibly a time dimension, etc). E.g. for a Nx2x3 field,
with `N` horizontal grid points k iterates over (1, 1) to (2, 3), meaning both the 2nd 
and 3rd dimension. To be used like

    for k in eachlayer(field)
        for (j, ring) in enumerate(eachring(field))
            for ij in ring
                field[ij, k]
                
With `vertical_only=true` (default) only checks whether the non-horizontal dimensions of the fields match.
E.g. one can loop over two fields of each n layers on different grids."""
function eachlayer(field1::AbstractField, fields::AbstractField...; vertical_only=true, kwargs...)
    fields_match(field1, fields...; vertical_only, kwargs...) || throw(DimensionMismatch(field1, fields...))
    return eachlayer(field1)
end

eachlayer(field::AbstractField) = CartesianIndices(size(field)[2:end])

"""$(TYPEDSIGNATURES)
Iterator over all 2D grid points of a field (or fields), i.e. the horizontal dimension only.
Intended use is like    

    for ij in eachgridpoint(field)
        field[ij]
        
for a 2D field. For a 3D+ field, the iterator will only index and loop over the horizontal grid points.
If `horizontal_only` is set to `true` (default), the function will only check whether the horizontal dimensions match.
"""
function eachgridpoint(field1::AbstractField, fields::AbstractField...; horizontal_only=true, kwargs...)
    fields_match(field1, fields...; horizontal_only, kwargs...) || throw(DimensionMismatch(field1, fields...))
    return eachgridpoint(field1)
end

eachgridpoint(field::AbstractField) = eachgridpoint(field.grid)

## CONSTRUCTORS

# from data array
# pass keyword `input_as` on to positional argument for dispatch
(::Type{F})(data::AbstractArray; input_as=Vector) where F<:AbstractField = F(data, input_as)

# if no nlat_half provided calculate it
function (::Type{F})(data::AbstractArray, input_as::Type{Vector}) where F<:AbstractField
    npoints = size(data, 1)                     # from 1st dim of data
    Grid = grid_type(F)                         # from type of field
    nlat_half = get_nlat_half(Grid, npoints)    # get nlat_half of Grid
    grid = Grid(nlat_half)
    return Field(data, grid)
end

# reduced fields cannot be represented as matrix
function (::Type{F})(data::AbstractArray, input_as::Type{Matrix})  where F<:AbstractReducedField
    error("Only full fields can be created from a Matrix.")
end

# full fields just need a reshape
function (::Type{F})(data::AbstractArray, input_as::Type{Matrix}) where F<:AbstractFullField
    # flatten the two horizontal dimensions into one, identical to vec(data) for data <: AbstractMatrix
    data_flat = reshape(data, :, size(data)[3:end]...)
    return F(data_flat, input_as=Vector)
end 

# CONVERSION to Array

# pass `as` keyword on to positional argument for dispatch
Base.Array(field::AbstractField; as=Vector) = Array(field, as)

# all horizontal data is stored as unravelled into a Vector
Base.Array(field::AbstractField, as::Type{Vector}) = Array(field.data)

# full grids can be reshaped into a matrix, reduced grids cannot
Base.Array(field::AbstractFullField, as::Type{Matrix}) = Array(reshape(field.data, :, get_nlat(field), size(field)[2:end]...))
Base.Array(field::AbstractReducedField, as::Type{Matrix}) = error("$(typeof(field)), i.e. data on a reduced grid, cannot be converted to a Matrix.")

Base.Vector(field::AbstractField2D) = Array(field, as=Vector)
Base.Matrix(field::AbstractField2D) = Array(field, as=Matrix)

for f in (:zeros, :ones, :rand, :randn)
    @eval begin
        # zeros(grid, nlayers...)
        function Base.$f(grid::AbstractGrid, k::Integer...)
            data = array_type(grid.architecture)($f(DEFAULT_NF, get_npoints(grid), k...))
            return Field(data, grid)
        end
            
        # zeros(NF, grid, nlayers...)
        function Base.$f(::Type{T}, grid::AbstractGrid, k::Integer...) where T
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
            nlat_half::Integer,
            k::Integer...;
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
            Grid = grid_type(F)
            grid = nonparametric_type(Grid)(nlat_half, architecture)
            data = array_type(architecture)($f(get_npoints(grid), k...))
            return Field(data, grid)
        end

        function Base.$f(
            ::Type{F},
            nlat_half::Integer,
            k::Integer...;
            architecture=DEFAULT_ARCHITECTURE(),
        ) where {F<:AbstractField{T}} where T
            Grid = grid_type(F)
            grid = nonparametric_type(Grid)(nlat_half, architecture)
            data = array_type(architecture)($f(T, get_npoints(grid), k...))
            return Field(data, grid)
        end

        function Base.$f(
            ::Type{F},
            grid::AbstractGrid,
            k::Integer...,
        ) where {F<:AbstractField{T}} where T
            data = array_type(F)($f(T, get_npoints(grid), k...))
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
    Grid = grid_type(F)
    grid = nonparametric_type(Grid)(nlat_half, DEFAULT_ARCHITECTURE())
    data = Array{DEFAULT_NF}(undef, get_npoints(grid), k...)
    return Field(data, grid)
end

# in case only number format is provided use Array and DEFAULT_ARCHITECTURE
function (::Type{F})(
    ::UndefInitializer,
    nlat_half::Integer,
    k::Integer...,
) where {F<:AbstractField{T}} where T
    Grid = grid_type(F)
    grid = nonparametric_type(Grid)(nlat_half, DEFAULT_ARCHITECTURE())
    data = Array{T}(undef, get_npoints(grid), k...)
    return Field(data, grid)
end

#  same as above but with N (ignored though as obtained from integer arguments)
function (::Type{F})(
    ::UndefInitializer,
    nlat_half::Integer,
    k::Integer...,
) where {F<:AbstractField{T, N}} where {T, N}
    Grid = grid_type(F)
    grid = nonparametric_type(Grid)(nlat_half, DEFAULT_ARCHITECTURE())
    data = Array{T}(undef, get_npoints(grid), k...)
    return Field(data, grid)
end

# in case ArrayType is provided use that!
function (::Type{F})(
    ::UndefInitializer,
    nlat_half::Integer,
    k::Integer...,
) where {F<:AbstractField{T, N, ArrayType}} where {T, N, ArrayType}
    Grid = grid_type(F)
    grid = nonparametric_type(Grid)(nlat_half, DEFAULT_ARCHITECTURE())
    data = nonparametric_type(ArrayType){T}(undef, get_npoints(grid), k...)
    return Field(data, grid)
end

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

# Views that return a Field again (need to retain all horizontal grid points, hence `:, 1` for example)
# view(array, :) unravels like array[:] does, hence "::Colon, i, args..." used to enforce one argument after :
# exception is view(vector, :) which preserves the vector structure, equivalent here is the Field2D
# TODO extend Base.view?
field_view(field::AbstractField,  c::Colon, i, args...) = Field(view(field.data, c, i, args...), field.grid)
field_view(field::AbstractField2D, c::Colon) = Field(view(field.data, c), field.grid)
field_view(field::AbstractField, args...) = view(field, args...)   # fallback to normal view

## BROADCASTING
# following https://docs.julialang.org/en/v1/manual/interfaces/#man-interfaces-broadcasting
import Base.Broadcast: BroadcastStyle, Broadcasted, DefaultArrayStyle

# new broadcasting style, AbstractArrays use just {N}, add additional parameter Grid to
# not broadcast between different grids which have therefore different (incompatible) broadcast styles
# important to not have parameter T (eltype/number format) here to broadcast
# automatically across the same field type but with different T
# e.g. FullGaussianField{Float32} and FullGaussianField{Float64}
struct FieldStyle{N, Grid} <: Broadcast.AbstractArrayStyle{N} end

# define broadcast style for Field from its parameters
Base.BroadcastStyle(::Type{F}) where {F<:AbstractField{T, N, ArrayType, Grid}} where {T, N, ArrayType, Grid} =
    FieldStyle{N, nonparametric_type(Grid)}()

# find a field within Broadcasted to reuse its grid
find_field(bc::Base.Broadcast.Broadcasted) = find_field(bc.args) 
find_field(args::Tuple) = find_field(find_field(args[1]), Base.tail(args)) 
find_field(x) = x 
find_field(::Tuple{}) = nothing 
find_field(a::AbstractField, rest) = a 
find_field(::Any, rest) = find_field(rest) 

# allocation for broadcasting via similar, reusing grid from the first field of the broadcast arguments
# e.g. field1 + field2 creates a new field that share the grid of field1
# 2 .+ field1 creates a new field that share the grid of field1
function Base.similar(bc::Broadcasted{FieldStyle{N, Grid}}, ::Type{T}) where {N, Grid, T}
    field = find_field(bc)
    ArrayType_ = nonparametric_type(typeof(field.data))
    new_data = ArrayType_{T}(undef, size(bc))
    old_grid = field.grid
    return Field(new_data, old_grid)
end

# ::Val{0} for broadcasting with 0-dimensional, ::Val{1} for broadcasting with vectors, etc
# when there's a dimension mismatch always choose the larger dimension
FieldStyle{N, Grid}(::Val{N}) where {N, Grid} = FieldStyle{N, Grid}()
FieldStyle{1, Grid}(::Val{2}) where {Grid} = FieldStyle{2, Grid}()
FieldStyle{1, Grid}(::Val{0}) where {Grid} = FieldStyle{1, Grid}()
FieldStyle{2, Grid}(::Val{3}) where {Grid} = FieldStyle{3, Grid}()
FieldStyle{2, Grid}(::Val{1}) where {Grid} = FieldStyle{2, Grid}()
FieldStyle{3, Grid}(::Val{4}) where {Grid} = FieldStyle{4, Grid}()
FieldStyle{3, Grid}(::Val{2}) where {Grid} = FieldStyle{3, Grid}()

## GPU (same but <: GPUArrays.AbstractGPUArrayStyle)
struct FieldGPUStyle{N, Grid} <: GPUArrays.AbstractGPUArrayStyle{N} end

# same as for FieldStyle but for constrain to ArrayType<:GPUArrays
function Base.BroadcastStyle(
    ::Type{F}
) where {F<:AbstractField{T, N, ArrayType, Grid}} where {T, N, ArrayType <: GPUArrays.AbstractGPUArray, Grid}
    return FieldGPUStyle{N, Grid}()
end

function Base.similar(bc::Broadcasted{FieldGPUStyle{N, Grid}}, ::Type{T}) where {N, Grid, T}
    field = find_field(bc)
    ArrayType_ = nonparametric_type(typeof(field.data))
    new_data = ArrayType_{T}(undef, size(bc))
    old_grid = field.grid
    return Field(new_data, old_grid)
end

# ::Val{0} for broadcasting with 0-dimensional, ::Val{1} for broadcasting with vectors, etc
# when there's a dimension mismatch always choose the larger dimension
FieldGPUStyle{N, Grid}(::Val{N}) where {N, Grid} = FieldGPUStyle{N, Grid}()
FieldGPUStyle{1, Grid}(::Val{2}) where {Grid}    = FieldGPUStyle{2, Grid}()
FieldGPUStyle{1, Grid}(::Val{0}) where {Grid}    = FieldGPUStyle{1, Grid}()
FieldGPUStyle{2, Grid}(::Val{3}) where {Grid}    = FieldGPUStyle{3, Grid}()
FieldGPUStyle{2, Grid}(::Val{1}) where {Grid}    = FieldGPUStyle{2, Grid}()
FieldGPUStyle{3, Grid}(::Val{4}) where {Grid}    = FieldGPUStyle{4, Grid}()
FieldGPUStyle{3, Grid}(::Val{2}) where {Grid}    = FieldGPUStyle{3, Grid}()

function KernelAbstractions.get_backend(
    field::F
) where {F <: AbstractField{T, N, ArrayType}} where {T, N, ArrayType <: GPUArrays.AbstractGPUArray}
    return KernelAbstractions.get_backend(field.data)
end

function Adapt.adapt_structure(to, field::Field{T, N, ArrayType, Grid}) where {T, N, ArrayType, Grid}
    adapted_data = adapt(to, field.data)
    if ismatching(field.grid, typeof(adapted_data))
        return Field(adapted_data, adapt(to, field.grid))
    else # if not matching, create new grid with other architecture
        #@warn "Adapting field to new architecture with $(typeof(adapted_data))"
        return Field(adapted_data, adapt(to, Grid(field.grid, architecture(typeof(adapted_data)))))
    end
end

Architectures.architecture(field::AbstractField) = architecture(field.grid)
Architectures.on_architecture(arch, field::AbstractField) = Adapt.adapt(array_type(arch), field)