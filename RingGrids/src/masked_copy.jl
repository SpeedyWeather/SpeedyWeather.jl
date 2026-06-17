"""$(TYPEDSIGNATURES)
Return a vector of linear indices into the ring grid at positions where `mask` is `false`.
The returned index vector maps from a contiguous masked-subset index (1…n) to the
corresponding flat grid-point index in the ring grid, so that

    field[ringgrid2masked[i], k] == masked_array[i, k]

after a `masked_copy!` call. The result lives on the same device as `mask`.
"""
function mask_to_indices(mask::AbstractField2D{<:Bool})
    # from ring grid to masked subset
    ij = on_architecture(architecture(mask), collect(1:length(mask)))
    ringgrid2masked = ij[~mask]
    return ringgrid2masked
end

"""$(TYPEDSIGNATURES)
Copy a subset of grid points from the ring-grid `Field` `src` into the plain array `dst`,
where `indices` maps each row of `dst` to the corresponding flat grid-point index in `src`.
`indices` is typically produced by [`mask_to_indices`](@ref).

The first dimension of `dst` and `indices` must match (the number of masked points);
remaining dimensions are treated as vertical layers and are copied verbatim.
"""
function masked_copy!(dst::AbstractArray, src::AbstractField, indices)
    arch = architecture(dst)
    @boundscheck size(indices, 1) == size(dst, 1) || throw(BoundsError(indices, dst))

    launch!(arch, ArrayWorkOrder, size(dst), masked_copy_kernel!, dst, src, indices)
    return dst
end

@kernel inbounds = true function masked_copy_kernel!(dst::AbstractArray, src::AbstractField, indices::AbstractVector{<:Integer})
    I = @index(Global, Cartesian)   # destination array index
    ij = indices[I[1]]              # corresponding source array index
    dst[I] = src[ij, I[2:end]...]                
end

"""$(TYPEDSIGNATURES)
Reverse of the field→array method: scatter the rows of the plain array `src` back into
the ring-grid `Field` `dst` at the positions given by `indices`. Grid points not referenced
by `indices` are left unchanged. `indices` is typically produced by [`mask_to_indices`](@ref)."""
function masked_copy!(dst::AbstractField, src::AbstractArray, indices)
    arch = architecture(src)
    @boundscheck size(indices, 1) == size(src, 1) || throw(BoundsError(indices, dst))

    launch!(arch, ArrayWorkOrder, size(src), masked_copy_kernel!, dst, src, indices)
    return dst
end

@kernel inbounds = true function masked_copy_kernel!(dst::AbstractField, src::AbstractArray, indices::AbstractVector{<:Integer})
    I = @index(Global, Cartesian)   # destination array index
    ij = indices[I[1]]              # corresponding source array index
    dst[ij, I[2:end]...] = src[I]   
end