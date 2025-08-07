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
ColumnField(grid::AbstractGrid, k...) = zeros(grid, k...)
ColumnField(::Type{T}, grid::AbstractGrid, k...) where T = zeros(T, grid, k...)
(::Type{<:ColumnField{T}})(data::AbstractArray, grid::AbstractGrid) where T = ColumnField(T.(data), grid)

# TYPES
Architectures.nonparametric_type(::Type{<:ColumnField}) = ColumnField
grid_type(::Type{ColumnField{T, N, A, G}}) where {T, N, A, G} = G
Architectures.array_type(::Type{ColumnField{T, N, A, G}}) where {T, N, A, G} = A

# conversion from Field 
transpose(field::Field) = transpose_safe(field)
transpose!(field::Field) = transpose_unsafe!(field, similar(transpose(field.data)))

# and back to Field 
transpose(field::ColumnField) = transpose_safe(field)
transpose!(field::ColumnField) = transpose_unsafe!(field, similar(transpose(field.data)))

# safe version which allocates a new field
function transpose_safe(field::Field)
    data_transposed = collect(permutedims(field.data, (2, 1, 3:ndims(field)...)))
    return ColumnField(data_transposed, field.grid)
end

# unsafe version that leaves the origional data corrupted (=transposed even though `field` doesn't indicate this)
# but returns a ColumnField view onto the same memory
function transpose_unsafe!(field::Field, scratch::AbstractArray)
    @boundscheck size(transpose(scratch)) == size(field.data) || throw(BoundsError(scratch, field.data))
    permutedims!(scratch, field.data, (2, 1, 3:ndims(field)...))    # transpose into scratch memory
    vec(field.data) .= vec(scratch)               # copy back (scratch memory is free again after this)
    return ColumnField(field.data, field.grid)    # view on the same data of field
end

# safe version which allocates a new field
function transpose_safe(field::ColumnField)
    data_transposed = collect(permutedims(field.data, (2, 1, 3:ndims(field)...)))
    return Field(data_transposed, field.grid)
end

# unsafe version that leaves the origional data corrupted (=transposed even though `field` doesn't indicate this)
# but returns a ColumnField view onto the same memory
function transpose_unsafe!(field::ColumnField, scratch::AbstractArray)
    @boundscheck size(transpose(scratch)) == size(field.data) || throw(BoundsError(scratch, field.data))
    permutedims!(scratch, field.data, (2, 1, 3:ndims(field)...))    # transpose into scratch memory
    vec(field.data) .= vec(scratch)               # copy back (scratch memory is free again after this)
    return Field(field.data, field.grid)    # view on the same data of field
end

# full grids can be reshaped into a matrix, reduced grids cannot
Base.Array(field::FullColumnField, as::Type{Matrix}) = Array(reshape(field.data, get_nlat(field), :, size(field)[2:end]...))

