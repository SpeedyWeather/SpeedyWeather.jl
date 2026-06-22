"""$(TYPEDSIGNATURES)
Return a vector of linear indices into the ring grid at positions where `mask` is `false`.
The returned index vector maps from a contiguous unmasked-subset index (1…n) to the
corresponding flat grid-point index in the ring grid, so that

    field[indices[i], k] == unmasked_array[i, k]

after a `copy_unmasked!` call. The result lives on the same device as `mask`."""
function unmasked_indices(mask::AbstractField2D{<:Bool})
    ij = on_architecture(architecture(mask), collect(1:length(mask)))
    return ij[.!mask.data]
end

"""$(TYPEDSIGNATURES)
Copy a subset of grid points from the ring-grid `Field` `src` into the plain array `dst`,
where `indices` maps each row of `dst` to the corresponding flat grid-point index in `src`.
`indices` is typically produced by [`unmasked_indices`](@ref).

The first dimension of `dst` and `indices` must match (the number of unmasked points);
remaining dimensions are treated as vertical layers and are copied verbatim.
"""
function copy_unmasked!(dst::AbstractArray, src::AbstractField, indices)
    arch = architecture(dst)
    @boundscheck size(indices, 1) == size(dst, 1) || throw(BoundsError(indices, dst))
    @boundscheck size(dst, 1) <= size(src, 1) || throw(BoundsError(indices, dst))

    # technically without this check the kernel is unsafe but `maximum` requires reading in the array ...
    # @boundscheck 1 <= maximum(indices) <= size(src, 1) || throw(BoundsError(indices, src))

    # call the gather kernel as dst is expected to be smaller as elements in src are masked
    launch!(arch, ArrayWorkOrder, size(dst), copy_unmasked_gather_kernel!, dst, src, indices)
    return dst
end

@kernel inbounds = true function copy_unmasked_gather_kernel!(dst, src, indices)
    I = @index(Global, NTuple)
    ij = indices[I[1]]
    dst[I...] = src[ij, I[2:end]...]
end

"""$(TYPEDSIGNATURES)
Reverse of the field→array method: scatter the rows of the plain array `src` back into
the ring-grid `Field` `dst` at the positions given by `indices`. Grid points not referenced
by `indices` are left unchanged. `indices` is typically produced by [`unmasked_indices`](@ref)."""
function copy_unmasked!(dst::AbstractField, src::AbstractArray, indices)
    arch = architecture(src)
    
    @boundscheck size(indices, 1) == size(src, 1) || throw(BoundsError(indices, dst))
    @boundscheck size(src, 1) <= size(dst, 1) || throw(BoundsError(indices, dst))

    # technically without this check the kernel is unsafe but `maximum` requires reading in the array ...
    # @boundscheck maximum(indices) <= size(dst, 1) || throw(BoundsError(indices, src))

    # call scatter kernel as dst is larger than src which only contains the unmasked elements in dst
    launch!(arch, ArrayWorkOrder, size(src), copy_unmasked_scatter_kernel!, dst, src, indices)
    return dst
end

@kernel inbounds = true function copy_unmasked_scatter_kernel!(dst, src, indices)
    I = @index(Global, NTuple)
    ij = indices[I[1]]
    dst[ij, I[2:end]...] = src[I...]
end
