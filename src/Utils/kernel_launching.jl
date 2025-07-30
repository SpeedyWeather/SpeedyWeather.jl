import ..device

# Support for 1D, To-DO: do we need higher dim workgroups like Oceananigans? Not in one horizontal level (because LTA/Rings), but maybe when we include vertical?
heuristic_workgroup(Wx) = min(Wx, 256)

# horizontal + vertical, 3D 
function heuristic_workgroup(Wxy, Wz) 
    minWy = min(Wz, 16)
    return (min(Wxy,cld(256,minWy)),minWy)
end 

function heuristic_workgroup(Wx, Wy, Wz)
    minWz = min(Wz, 8)
    minWy = min(Wy, 8)
    return (min(Wx,cld(cld(256,minWy),minWz)),minWy,minWz)
end

# WORK ORDER TYPES FOR TYPE-BASED DISPATCH
"""
Abstract type for different work order patterns in kernel launching.
"""
abstract type AbstractWorkOrder end

"""
Work order for kernels over spectral number 'lm' and layers.
"""
struct SpectralWorkOrder <: AbstractWorkOrder end

"""
Work order for kernels over RingGrids 'ij' indices and layers.
"""
struct RingGridWorkOrder <: AbstractWorkOrder end

"""
Work order for kernels over spectral number 'lm' and layers, excluding the lm=1 element.
"""
struct SpectralInnerWorkOrder <: AbstractWorkOrder end

"""
Work order for kernels over the diagonal of a (m+1)×m or m×m LowerTriangularArray.
"""
struct DiagonalWorkOrder <: AbstractWorkOrder end

"""
Work order for kernels over a regular 3D array.
"""
struct Array3DWorkOrder <: AbstractWorkOrder end

"""
Work order for kernels over linear/eachindex iteration.
"""
struct LinearWorkOrder <: AbstractWorkOrder end

# TYPE-BASED DISPATCH METHODS

"""
$(TYPEDSIGNATURES)
Returns the `workgroup` and `worksize` for launching a kernel with spectral work order.
"""
function work_layout(::Type{SpectralWorkOrder}, worksize::NTuple{N, Int}) where N
    workgroup = heuristic_workgroup(worksize...)
    return workgroup, worksize
end

"""
$(TYPEDSIGNATURES)
Returns the `workgroup` and `worksize` for launching a kernel with ring grid work order.
"""
function work_layout(::Type{RingGridWorkOrder}, worksize::NTuple{N, Int}) where N
    workgroup = heuristic_workgroup(worksize...)
    return workgroup, worksize
end

"""
$(TYPEDSIGNATURES)
Returns the `workgroup` and `worksize` for launching a kernel with spectral inner work order,
excluding the lm=1 element.
"""
function work_layout(::Type{SpectralInnerWorkOrder}, worksize::NTuple{N, Int}) where N
    return heuristic_workgroup(worksize[1]-2, worksize[2:end]...), (worksize[1]-2, worksize[2:end]...)
end

"""
$(TYPEDSIGNATURES)
Returns the `workgroup` and `worksize` for launching a kernel over the diagonal of a LowerTriangularArray.
"""
function work_layout(::Type{DiagonalWorkOrder}, worksize::NTuple{N, Int}) where N
    return heuristic_workgroup(worksize[2], worksize[3:end]...), (worksize[2], worksize[3:end]...)
end

"""
$(TYPEDSIGNATURES)
Returns the `workgroup` and `worksize` for launching a kernel over a regular 3D array.
"""
function work_layout(::Type{Array3DWorkOrder}, worksize::NTuple{N, Int}) where N
    return heuristic_workgroup(worksize...), worksize
end

"""
$(TYPEDSIGNATURES)
Returns the `workgroup` and `worksize` for launching a kernel with linear iteration.
"""
function work_layout(::Type{LinearWorkOrder}, worksize::NTuple{N, Int}) where N
    return heuristic_workgroup(prod(worksize)), (prod(worksize),)
end

"""
$(TYPEDSIGNATURES)
Configure `kernel!` to launch over the `dims` of `grid` on
the architecture `arch`.

# Arguments
============

- `arch`: The architecture on which the kernel will be launched.
- `dims_type`: The dimensions on which the kernel will be executed, a subtype of `AbstractWorkOrder`.
- `worksize`: The size that defines the work distribution.
- `kernel!`: The kernel function to be executed.
"""
@inline function configure_kernel(arch, work_order::Type{<:AbstractWorkOrder}, worksize, kernel!)
    workgroup, worksize = work_layout(work_order, worksize)
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

# Inner interface for type-based dispatch (preferred)
@inline function _launch!(arch, work_order::Type{<:AbstractWorkOrder}, worksize, kernel!, kernel_args...)

    loop!, worksize = configure_kernel(arch, work_order, worksize, kernel!)
    
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