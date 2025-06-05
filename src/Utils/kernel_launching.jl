import ..device

# Support for 1D, To-DO: do we need higher dim workgroups like Oceananigans? Not in one horizontal level (because LTA/Rings), but maybe when we include vertical?
heuristic_workgroup(Wx) = min(Wx, 256)

# horizontal + vertical, 3D 
function heuristic_workgroup(Wxy, Wz) 
    minWy = min(Wz, 16)
    return (min(Wxy,cld(256,minWy)),minWy)
end 

# TO-DO: in the future we might want to distinguish between LTA/RingGrids/different grids in some way
"""
($(TYPEDSIGNATURES))
Returns the `workgroup` and `worksize` for launching a kernel over `dims`
on `grid` that excludes peripheral nodes.
The `workgroup` is a tuple specifying the threads per block in each
dimension. The `worksize` specifies the range of the loop in each dimension.
"""
function work_layout(dims_type::Symbol, worksize::NTuple{N, Int}) where N

    # To-Do: introduce `dims_type` e.g: `:mk` for kernel over m, k or `:lmk` for kernel over lm, k, ....

    if dims_type == :lmk # kernel over sph number 'lm' and layers 
        # compare how jack has done that in 
        # number of (M, k) items?
        workgroup = heuristic_workgroup(worksize...)
    elseif dims_type == :lmk_inner_points    # kernel over sph number 'lm' and layers, but leaving out the lm=1 element
        # number of (M, k) items?
        return heuristic_workgroup(worksize[1]-2, worksize[2:end]...), (worksize[1]-2, worksize[2:end]...)
    elseif dims_type == :diagonal # kernel just over the diagonal of a (m+1)x m or m x m LTA
        return heuristic_workgroup(worksize[2], worksize[3:end]...), (worksize[2], worksize[3:end]...)
    elseif dims_type == :lastrow # kernel just over the last row of a (m+1)x m or m x m LTA
        return heuristic_workgroup(worksize[2], worksize[3:end]...), (worksize[2], worksize[3:end]...)
    elseif dims_type == :lmk_main # kernel over everything in a LTA that's not the last row or diagonal
        N_elements = worksize[2]*(worksize[2]+1)÷2 - worksize[2] 
        return heuristic_workgroup(N_elements, worksize[3:end]...), (N_elements, worksize[3:end]...)
    else 
        error("Not yet implemented")
    end 
    return workgroup, worksize
end



"""
$(TYPEDSIGNATURES)
Configure `kernel!` to launch over the `dims` of `grid` on
the architecture `arch`.

# Arguments
============

- `arch`: The architecture on which the kernel will be launched.
- `dims_type`: The dimensions on which the kernel will be executed (TO-DO: in the future we might distinguish LTA/RingGrids here)
- `worksize`: The size that defines the work distribution.
- `kernel!`: The kernel function to be executed.
"""
@inline function configure_kernel(arch, dims_type, worksize, kernel!)

    workgroup, worksize = work_layout(dims_type, worksize)

    dev  = device(arch)
    loop = kernel!(dev, workgroup, worksize)

    return loop, worksize
end
  
"""
launch!(arch, dims_type, worksize, kernel!, kernel_args...)

Launches `kernel!` with arguments `kernel_args`
over the `dims_type` with `worksize` on the architecture `arch`.
Kernels run on the default stream.

See [configure_kernel](@ref) for more information.
"""
@inline launch!(args...; kwargs...) = _launch!(args...; kwargs...)

# Inner interface
@inline function _launch!(arch, dims_type, worksize, kernel!, kernel_args...)

    loop!, worksize = configure_kernel(arch, dims_type, worksize, kernel!)

    # Don't launch kernels with no size
    haswork = if worksize isa Number
        worksize > 0       
    else
        true
    end

    if haswork
        loop!(kernel_args...)
    end

    return nothing
end