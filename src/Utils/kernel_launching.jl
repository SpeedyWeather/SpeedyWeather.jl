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
function work_layout(data_type, worksize::NTuple{N, Int}) where N
    workgroup = heuristic_workgroup(worksize...)
    return workgroup, worksize
end

"""
$(TYPEDSIGNATURES)
Configure `kernel!` to launch over the `dims` of `grid` on
the architecture `arch`.

# Arguments
============

- `arch`: The architecture on which the kernel will be launched.
- `data_type`: The data_type on which the kernel will be executed (TO-DO: in the future we might distinguish LTA/RingGrids here)
- `workspec`: The workspec that defines the work distribution.
- `kernel!`: The kernel function to be executed.
"""
@inline function configure_kernel(arch, data_type, workspec, kernel!)

    workgroup, worksize = work_layout(data_type, workspec)

    dev  = device(arch)
    loop = kernel!(dev, workgroup, worksize)

    return loop, worksize
end
  
"""
launch!(arch, data_type, workspec, kernel!, kernel_args...)

Launches `kernel!` with arguments `kernel_args`
over the `dims` of `data_type` on the architecture `arch`.
Kernels run on the default stream.

See [configure_kernel](@ref) for more information.
"""
@inline launch!(args...; kwargs...) = _launch!(args...; kwargs...)

# Inner interface
@inline function _launch!(arch, data_type, workspec, kernel!, kernel_args...)

    loop!, worksize = configure_kernel(arch, data_type, workspec, kernel!)

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
