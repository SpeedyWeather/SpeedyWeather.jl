# REORDERING: Ring to Nested or Matrix and vice versa
abstract type AbstractOrder end
struct RingOrder <: AbstractOrder end
struct NestedOrder <: AbstractOrder end
struct MatrixOrder <: AbstractOrder end

reorder(::RingOrder,   ij, grid::AbstractGrid) = nest2ring(ij, grid)
reorder(::NestedOrder, ij, grid::AbstractGrid) = ring2nest(ij, grid)
reorder(::MatrixOrder, ij, grid::AbstractGrid) = ring2xy(ij, grid)
reorder(order, field::AbstractField) = reorder!(similar(field), order, field)

function reorder!(
    out::AbstractField,
    order::AbstractOrder,
    field::AbstractField,
)
    @boundscheck out.grid == field.grid || throw(BoundsError("Reordering requires identical grids, got $(out.grid) and $(field.grid)."))
    @assert ispow2(get_nlat_half(field.grid)) "Reordering only supported for nlat_half power of 2, got $(get_nlat_half(field.grid))."

    arch = architecture(field)
    # Ensure worksize is always 2D by adding layer dimension dim=1
    worksize = ndims(field) == 1 ? (size(field, 1), 1) : size(field)
    launch!(arch, RingGridWorkOrder, worksize, reorder_kernel!, out, order, field)
    return out
end

@kernel inbounds=true function reorder_kernel!(out, @Const(order), @Const(field))
    ij, k = @index(Global, NTuple)
    # TODO the recomputes the reordering for every layer k, maybe distribute only over ij?
    out_indices = reorder(order, ij, field.grid)    
    out[out_indices, k] = field[ij, k]
end

ring_order(  field::AbstractField) = reorder(RingOrder(),   field)
nested_order(field::AbstractField) = reorder(NestedOrder(), field)
matrix_order(field::AbstractField) = reorder(MatrixOrder(), field)

ring_order!(  out, field::AbstractField) = reorder!(out, RingOrder(),   field)
nested_order!(out, field::AbstractField) = reorder!(out, NestedOrder(), field)
matrix_order!(out, field::AbstractField) = reorder!(out, MatrixOrder(), field)

Matrix(field::AbstractField) = reshape(field.data, matrix_size(field)...)

