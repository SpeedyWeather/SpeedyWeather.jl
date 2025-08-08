import LinearAlgebra.transpose
import LinearAlgebra.transpose!

"""
ColumnField is a version of `Field` with data stored in a column-major format. 
It is used e.g. to represent data on a vertical column of a grid in column-based
parametrizations. 
$(TYPEDFIELDS)

A ColumnField can be easily created from a `Field` by transposing it:

```julia
field = Field(grid, k...)
column_field = transpose(field)
```
"""
struct ColumnField{T, N, ArrayType <: AbstractArray, Grid <: AbstractGrid} <: AbstractField{T, N, ArrayType, Grid}
    data::ArrayType
    grid::Grid

    function ColumnField(data, grid)
        data_matches_grid(data, grid; horizontal_dim=2) || throw(DimensionMismatch(data, grid))
        return new{eltype(data), ndims(data), typeof(data), typeof(grid)}(data, grid)
    end
end

const ColumnField2D = ColumnField{T, 1} where T
const ColumnField3D = ColumnField{T, 2} where T
const ColumnField4D = ColumnField{T, 3} where T

const FullColumnField = ColumnField{T, N, ArrayType, Grid} where {T, N, ArrayType, Grid<:AbstractFullGrid}
const ReducedColumnField = ColumnField{T, N, ArrayType, Grid} where {T, N, ArrayType, Grid<:AbstractReducedGrid}

# default constructors 
ColumnField(grid::AbstractGrid, k...) = transpose(zeros(grid, k...))
ColumnField(::Type{T}, grid::AbstractGrid, k...) where T = transpose(zeros(T, grid, k...))
(::Type{<:ColumnField{T}})(data::AbstractArray, grid::AbstractGrid) where T = ColumnField(T.(data), grid)

# TYPES
Architectures.nonparametric_type(::Type{<:ColumnField}) = ColumnField
grid_type(::Type{ColumnField{T, N, A, G}}) where {T, N, A, G} = G
Architectures.array_type(::Type{ColumnField{T, N, A, G}}) where {T, N, A, G} = A

# CONVERSION from Field 
LinearAlgebra.transpose(field::Field) = transpose_safe(field)
LinearAlgebra.transpose!(field::Field) = transpose_unsafe!(field, similar(field.data))

# and back to Field 
LinearAlgebra.transpose(field::ColumnField) = transpose_safe(field)
LinearAlgebra.transpose!(field::ColumnField) = transpose_unsafe!(field, similar(field.data))

# safe version which allocates a new field
function transpose_safe(field::Field)
    data_transposed = collect(permutedims(field.data, (2, 1, 3:ndims(field)...)))
    return ColumnField(data_transposed, field.grid)
end

# tranpose is defined here as swapping first two dimensions 
functison _size_of_transpose(data::AbstractArray)
    IJ, K, L... = size(data)
    return (K, IJ, L...)
end

# for vectors (i.e. Field2D) return original size
_size_of_transpose(data::AbstractVector) = size(data)

# unsafe version that leaves the origional data corrupted (=transposed even though `field` doesn't indicate this)
# but returns a ColumnField view onto the same memory
function transpose_unsafe!(field::Field, scratch::AbstractArray)
    @boundscheck length(scratch)) == length(field.data) || throw(BoundsError(scratch, field.data))
    scratch_reshaped = reshape(scratch, _size_of_transpose(field.data))     # so scratch input can be of size(field) or its tranpose
    permutedims!(scratch_reshaped, field.data, (2, 1, 3:ndims(field)...))   # transpose into scratch memory
    vec(field.data) .= vec(scratch_reshaped)                                # copy back (scratch memory is free again after this)
    return ColumnField(field.data, field.grid)                              # view on the same data of field
end

# safe version which allocates a new field
function transpose_safe(field::ColumnField)
    data_transposed = collect(permutedims(field.data, (2, 1, 3:ndims(field)...)))
    return Field(data_transposed, field.grid)
end

# unsafe version that leaves the origional data corrupted (=transposed even though `field` doesn't indicate this)
# but returns a ColumnField view onto the same memory
function transpose_unsafe!(field::ColumnField, scratch::AbstractArray)
    @boundscheck length(scratch) == length(field.data) || throw(BoundsError(scratch, field.data))
    scratch_reshaped = reshape(scratch, _size_transpose(field.data))        # so scratch input can be of size(field) or its tranpose
    permutedims!(scratch_reshaped, field.data, (2, 1, 3:ndims(field)...))   # transpose into scratch memory
    vec(field.data) .= vec(scratch_reshaped)                                # copy back (scratch memory is free again after this)
    return Field(field.data, field.grid)                                    # view on the same data of field
end

# SIMILAR AND UNDEF CONSTRUCTORS
# data with same type T but new size (=new grid)
function Base.similar(
    field::ColumnField,
    nlayers::Integer,
    nlat_half::Integer,
    k::Integer...
)
    # use same architecture though
    new_grid = typeof(field.grid)(nlat_half, field.grid.architecture)
    similar_data = similar(field.data, nlayers, get_npoints(new_grid), k...)
    return ColumnField(similar_data, new_grid)
end

# data with new type T and new size
function Base.similar(
    field::ColumnField,
    ::Type{Tnew},
    nlayers::Integer,
    nlat_half::Integer,
    k::Integer...
) where Tnew
    # use same architecture though
    new_grid = typeof(field.grid)(nlat_half, field.grid.architecture)
    similar_data = similar(field.data, Tnew, nlayers, get_npoints(new_grid), k...)
    return ColumnField(similar_data, new_grid)
end

# general version with ArrayType{T, N}(undef, ...) generator
function (::Type{F})(
    ::UndefInitializer,
    nlayers::Integer,
    nlat_half::Integer,
    k::Integer...,
) where {F<:ColumnField{T, N, ArrayType, Grid}} where {T, N, ArrayType, Grid<:AbstractGrid{Architecture}} where Architecture
    # this (inevitably) creates a new grid and architecture instance
    grid = nonparametric_type(Grid)(nlat_half, Architecture())
    # TODO this should allocate on the device
    data = nonparametric_type(ArrayType){T}(undef, nlayers, get_npoints(grid), k...)
    return ColumnField(data, grid)
end

# in case Architecture is not provided use DEFAULT_ARCHITECTURE
function (::Type{F})(
    ::UndefInitializer,
    nlayers::Integer,
    nlat_half::Integer,
    k::Integer...,
) where {F<:ColumnField{T, N, ArrayType, Grid}} where {T, N, ArrayType, Grid<:AbstractGrid}
    grid = nonparametric_type(Grid)(nlat_half, DEFAULT_ARCHITECTURE())
    data = nonparametric_type(ArrayType){T}(undef, nlayers, get_npoints(grid), k...)
    return ColumnField(data, grid)
end

# in case only grid is provided (e.g. FullGaussianField) use Float64, Array, DEFAULT_ARCHITECTURE
function (::Type{F})(
    ::UndefInitializer,
    nlayers::Integer,
    nlat_half::Integer,
    k::Integer...,
) where {F<:ColumnField}
    Grid = grid_type(F)
    grid = nonparametric_type(Grid)(nlat_half, DEFAULT_ARCHITECTURE())
    data = Array{DEFAULT_NF}(undef, nlayers, get_npoints(grid), k...)
    return ColumnField(data, grid)
end

# in case only number format is provided use Array and DEFAULT_ARCHITECTURE
function (::Type{F})(
    ::UndefInitializer,
    nlayers::Integer,
    nlat_half::Integer,
    k::Integer...,
) where {F<:ColumnField{T}} where T
    Grid = grid_type(F)
    grid = nonparametric_type(Grid)(nlat_half, DEFAULT_ARCHITECTURE())
    data = Array{T}(undef, nlayers, get_npoints(grid), k...)
    return ColumnField(data, grid)
end

#  same as above but with N (ignored though as obtained from integer arguments)
function (::Type{F})(
    ::UndefInitializer,
    nlayers::Integer,
    nlat_half::Integer,
    k::Integer...,
) where {F<:ColumnField{T, N}} where {T, N}
    Grid = grid_type(F)
    grid = nonparametric_type(Grid)(nlat_half, DEFAULT_ARCHITECTURE())
    data = Array{T}(undef, nlayers, get_npoints(grid), k...)
    return ColumnField(data, grid)
end

# in case ArrayType is provided use that!
function (::Type{F})(
    ::UndefInitializer,
    nlayers::Integer,
    nlat_half::Integer,
    k::Integer...,
) where {F<:ColumnField{T, N, ArrayType}} where {T, N, ArrayType}
    Grid = grid_type(F)
    grid = nonparametric_type(Grid)(nlat_half, DEFAULT_ARCHITECTURE())
    data = nonparametric_type(ArrayType){T}(undef, nlayers, get_npoints(grid), k...)
    return ColumnField(data, grid)
end

## Some ARITHMETICS with regular Fields 
add!(field::Field, column_field::ColumnField) = field .+= transpose(column_field)
add!(column_field::ColumnField, field::Field) = column_field .+= transpose(field)

## BROADCASTING (main functionality defined in field.jl)

# define broadcast style for ColumnField from its parameters
Base.BroadcastStyle(::Type{F}) where {F<:ColumnField{T, N, ArrayType, Grid}} where {T, N, ArrayType, Grid} =
    FieldStyle{N, nonparametric_type(Grid), true}()

# allocation for broadcasting via similar, reusing grid from the first field of the broadcast arguments
# e.g. field1 + field2 creates a new field that share the grid of field1
# 2 .+ field1 creates a new field that share the grid of field1
function Base.similar(bc::Broadcasted{FieldStyle{N, Grid, true}}, ::Type{T}) where {N, Grid, T}
    field = find_field(bc)
    ArrayType_ = nonparametric_type(typeof(field.data))
    new_data = ArrayType_{T}(undef, size(bc))
    old_grid = field.grid
    return ColumnField(new_data, old_grid)
end

# same as for FieldStyle but for constrain to ArrayType<:GPUArrays
function Base.BroadcastStyle(
    ::Type{F}
) where {F<:ColumnField{T, N, ArrayType, Grid}} where {T, N, ArrayType <: GPUArrays.AbstractGPUArray, Grid}
    return FieldGPUStyle{N, Grid, true}()
end

function Base.similar(bc::Broadcasted{FieldGPUStyle{N, Grid, true}}, ::Type{T}) where {N, Grid, T}
    field = find_field(bc)
    ArrayType_ = nonparametric_type(typeof(field.data))
    new_data = ArrayType_{T}(undef, size(bc))
    old_grid = field.grid
    return ColumnField(new_data, old_grid)
end

function Adapt.adapt_structure(to, field::ColumnField{T, N, ArrayType, Grid}) where {T, N, ArrayType, Grid}
    adapted_data = adapt(to, field.data)
    if ismatching(field.grid, typeof(adapted_data))
        return ColumnField(adapted_data, adapt(to, field.grid))
    else # if not matching, create new grid with other architecture
        #@warn "Adapting field to new architecture with $(typeof(adapted_data))"
        return ColumnField(adapted_data, adapt(to, Grid(field.grid, architecture(typeof(adapted_data)))))
    end
end
