function mask_to_indices(mask::AbstractField{<:Bool})
    indices = on_architecture(architecture(mask), collect(1:length(mask)))
    return indices[mask]
end

"""To be used like

    grid = HEALPixGrid(24)
    src = rand(grid)
    mask = rand(Bool, grid)
    indices = mask_to_indices(mask)
    dst = similar(indices, eltype(src))
    masked_copy!(dst, src, indices)
"""
function masked_copy!(dst, src::AbstractField, indices)
    arch = architecture(dst)
    launch!(arch, LinearWorkOrder, (length(dst),), masked_copy_kernel!, dst, src, indices)
    return dst
end

@kernel inbounds = true function masked_copy_kernel!(dst, src, indices)
    i = @index(Global, Linear)     # destination array index
    k = indices[i]                 # corresponding source array index
    dst[i] = src[k]                
end