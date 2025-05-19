import SpeedyWeather.SpectralTransforms: ismatching

# ∇!, function barrier for GPU arrays 
∇!(
    dpdx::LowerTriangularArray{T,N,<:CuArray},     # Output: zonal gradient
    dpdy::LowerTriangularArray{T,N,<:CuArray},     # Output: meridional gradient
    p::LowerTriangularArray{T,N,<:CuArray},        # Input: spectral coefficients
    S::SpectralTransform;           # includes precomputed arrays
    kwargs...
    ) where {T,N} = ∇_KA!(dpdx, dpdy, p, S; kwargs...)

# this is available directly for all array types for testing, later the function barrier above might be removed 
function ∇_KA!(
    dpdx::LowerTriangularArray,     # Output: zonal gradient
    dpdy::LowerTriangularArray,     # Output: meridional gradient
    p::LowerTriangularArray,        # Input: spectral coefficients
    S::SpectralTransform;           # includes precomputed arrays
    radius = DEFAULT_RADIUS,        # scale with radius if provided, otherwise unit sphere
)
    (; grad_y1, grad_y2, lm2m_indices) = S
    @boundscheck ismatching(S, p) || throw(DimensionMismatch(S, p))

    dpdx .= @. complex(0, lm2m_indices - 1)*p

    # first and last element aren't covered by the kernel because they would access p[0], p[end+1]
    dpdy[1,:] .= grad_y2[1] .* p[2, :]
    launch!(S.architecture, :lmk_inner_points, size(dpdy), dpdy_kernel!, dpdy, p, grad_y1, grad_y2)
    dpdy[end,:] .= grad_y1[end] .* p[end-1,:]

    # 1/radius factor if not unit sphere
    if radius != 1
        R⁻¹ = inv(radius)
        dpdx .*= R⁻¹
        dpdy .*= R⁻¹
    end

    return dpdx, dpdy
end

@kernel inbounds=true function dpdy_kernel!(dpdy, p, @Const(grad_y1), @Const(grad_y2))
    I = @index(Global, Cartesian)
    lm = I[1] + 1   # +1 because we leave out the lm=1 element
    k = ndims(p) == 1 ? CartesianIndex() : I[2]

    gy1 = grad_y1[lm]
    gy2 = grad_y2[lm]
    
    # compared to the old CPU only version, some of the gy1 and gy2 are zero 
    # that's why we don't need to check for l==m (diagonal) or l==p.m (last row)
    dpdy[lm, k] = gy1*p[lm-1, k] + gy2*p[lm+1, k]
end

# _divergence! function barrier for GPU arrays 
_divergence!(
    kernel,
    div::LowerTriangularArray{T,N,<:CuArray},
    u::LowerTriangularArray{T,N,<:CuArray},
    v::LowerTriangularArray{T,N,<:CuArray},
    S::SpectralTransform;
    radius = DEFAULT_RADIUS,
) = _divergence_KA!(kernel, div, u, v, S; radius)
    

function divergence_KA!(
    div::LowerTriangularArray,
    u::LowerTriangularArray,
    v::LowerTriangularArray,
    S::SpectralTransform;
    flipsign::Bool=false,
    add::Bool=false,
    kwargs...,
)
    # Dispatch to the appropriate specialized version based on flipsign and add
    if flipsign
        if add
            _divergence_KA_flipsign_add!(div, u, v, S; kwargs...)
        else
            _divergence_KA_flipsign_noadd!(div, u, v, S; kwargs...)
        end
    else
        if add
            _divergence_KA_noflipsign_add!(div, u, v, S; kwargs...)
        else
            _divergence_KA_noflipsign_noadd!(div, u, v, S; kwargs...)
        end
    end
    return div
end

# Common implementation for all four cases
function _divergence_KA_common!(
    div::LowerTriangularArray,
    u::LowerTriangularArray,
    v::LowerTriangularArray,
    S::SpectralTransform,
    kernel_func;
    radius = DEFAULT_RADIUS,
)
    (; grad_y_vordiv1, grad_y_vordiv2, lm2m_indices) = S

    @boundscheck ismatching(S, div) || throw(DimensionMismatch(S, div))

    # /radius scaling if not unit sphere
    r_factor = radius != 1 ? inv(radius) : 1.0

    return div, grad_y_vordiv1, grad_y_vordiv2, lm2m_indices, r_factor
end

# Case 1: flipsign=false, add=false (default) - div = ∇⋅(u, v)
function _divergence_KA_noflipsign_noadd!(
    div::LowerTriangularArray,
    u::LowerTriangularArray,
    v::LowerTriangularArray,
    S::SpectralTransform;
    radius = DEFAULT_RADIUS,
)
    div, grad_y_vordiv1, grad_y_vordiv2, lm2m_indices, r_factor = 
        _divergence_KA_common!(div, u, v, S, nothing; radius=radius)

    # First element - special case
    for k in axes(div, 2)
        div[1, k] = grad_y_vordiv2[1] * v[2, k]
    end
    
    # Main kernel
    launch!(S.architecture, :lmk_inner_points, size(div), _divergence_kernel_noflipsign_noadd!, 
            div, u, v, grad_y_vordiv1, grad_y_vordiv2, lm2m_indices)

    # Last element
    div[end,:] .= 0

    # Apply radius scaling if needed
    if radius != 1
        div .*= r_factor
    end

    return div
end

# Case 2: flipsign=false, add=true - div += ∇⋅(u, v)
function _divergence_KA_noflipsign_add!(
    div::LowerTriangularArray,
    u::LowerTriangularArray,
    v::LowerTriangularArray,
    S::SpectralTransform;
    radius = DEFAULT_RADIUS,
)
    div, grad_y_vordiv1, grad_y_vordiv2, lm2m_indices, r_factor = 
        _divergence_KA_common!(div, u, v, S, nothing; radius=radius)

    # First element - special case
    for k in axes(div, 2)
        div[1, k] += grad_y_vordiv2[1] * v[2, k]
    end
    
    # Main kernel
    launch!(S.architecture, :lmk_inner_points, size(div), _divergence_kernel_noflipsign_add!, 
            div, u, v, grad_y_vordiv1, grad_y_vordiv2, lm2m_indices)

    # Last element - no change needed for add=true

    # Apply radius scaling if needed
    if radius != 1
        div .*= r_factor
    end

    return div
end

# Case 3: flipsign=true, add=false - div = -∇⋅(u, v)
function _divergence_KA_flipsign_noadd!(
    div::LowerTriangularArray,
    u::LowerTriangularArray,
    v::LowerTriangularArray,
    S::SpectralTransform;
    radius = DEFAULT_RADIUS,
)
    div, grad_y_vordiv1, grad_y_vordiv2, lm2m_indices, r_factor = 
        _divergence_KA_common!(div, u, v, S, nothing; radius=radius)

    # First element - special case
    for k in axes(div, 2)
        div[1, k] = -grad_y_vordiv2[1] * v[2, k]
    end
    
    # Main kernel
    launch!(S.architecture, :lmk_inner_points, size(div), _divergence_kernel_flipsign_noadd!, 
            div, u, v, grad_y_vordiv1, grad_y_vordiv2, lm2m_indices)

    # Last element
    div[end,:] .= 0

    # Apply radius scaling if needed
    if radius != 1
        div .*= r_factor
    end

    return div
end

# Case 4: flipsign=true, add=true - div -= ∇⋅(u, v)
function _divergence_KA_flipsign_add!(
    div::LowerTriangularArray,
    u::LowerTriangularArray,
    v::LowerTriangularArray,
    S::SpectralTransform;
    radius = DEFAULT_RADIUS,
)
    div, grad_y_vordiv1, grad_y_vordiv2, lm2m_indices, r_factor = 
        _divergence_KA_common!(div, u, v, S, nothing; radius=radius)

    # First element - special case
    for k in axes(div, 2)
        div[1, k] -= grad_y_vordiv2[1] * v[2, k]
    end
    
    # Main kernel
    launch!(S.architecture, :lmk_inner_points, size(div), _divergence_kernel_flipsign_add!, 
            div, u, v, grad_y_vordiv1, grad_y_vordiv2, lm2m_indices)

    # Last element - no change needed for add=true

    # Apply radius scaling if needed
    if radius != 1
        div .*= r_factor
    end

    return div
end

# Specialized kernels for each case
@kernel inbounds=true function _divergence_kernel_noflipsign_noadd!(div, @Const(u), @Const(v), @Const(grad_y_vordiv1), @Const(grad_y_vordiv2), @Const(lm2m_indices))
    I = @index(Global, Cartesian)
    lm = I[1] + 1   # +1 because we leave out the lm=1 element
    
    k = ndims(div) == 1 ? CartesianIndex() : I[2]
    m = lm2m_indices[lm] 
    
    # div = a-b+c where a = im*(m-1)*u, b = grad_y_vordiv1*v[lm-1], c = grad_y_vordiv2*v[lm+1]
    ∂u∂λ = im*(m-1)*u[lm, k]
    ∂v∂θ1 = grad_y_vordiv1[lm]*v[lm-1, k]
    ∂v∂θ2 = grad_y_vordiv2[lm]*v[lm+1, k]
    div[lm, k] = ∂u∂λ - ∂v∂θ1 + ∂v∂θ2
end

@kernel inbounds=true function _divergence_kernel_noflipsign_add!(div, @Const(u), @Const(v), @Const(grad_y_vordiv1), @Const(grad_y_vordiv2), @Const(lm2m_indices))
    I = @index(Global, Cartesian)
    lm = I[1] + 1   # +1 because we leave out the lm=1 element
    
    k = ndims(div) == 1 ? CartesianIndex() : I[2]
    m = lm2m_indices[lm] 
    
    # div += a-b+c where a = im*(m-1)*u, b = grad_y_vordiv1*v[lm-1], c = grad_y_vordiv2*v[lm+1]
    ∂u∂λ = im*(m-1)*u[lm, k]
    ∂v∂θ1 = grad_y_vordiv1[lm]*v[lm-1, k]
    ∂v∂θ2 = grad_y_vordiv2[lm]*v[lm+1, k]
    div[lm, k] += ∂u∂λ - ∂v∂θ1 + ∂v∂θ2
end

@kernel inbounds=true function _divergence_kernel_flipsign_noadd!(div, @Const(u), @Const(v), @Const(grad_y_vordiv1), @Const(grad_y_vordiv2), @Const(lm2m_indices))
    I = @index(Global, Cartesian)
    lm = I[1] + 1   # +1 because we leave out the lm=1 element
    
    k = ndims(div) == 1 ? CartesianIndex() : I[2]
    m = lm2m_indices[lm] 
    
    # div = -(a-b+c) where a = im*(m-1)*u, b = grad_y_vordiv1*v[lm-1], c = grad_y_vordiv2*v[lm+1]
    ∂u∂λ = im*(m-1)*u[lm, k]
    ∂v∂θ1 = grad_y_vordiv1[lm]*v[lm-1, k]
    ∂v∂θ2 = grad_y_vordiv2[lm]*v[lm+1, k]
    div[lm, k] = -(∂u∂λ - ∂v∂θ1 + ∂v∂θ2)
end

@kernel inbounds=true function _divergence_kernel_flipsign_add!(div, @Const(u), @Const(v), @Const(grad_y_vordiv1), @Const(grad_y_vordiv2), @Const(lm2m_indices))
    I = @index(Global, Cartesian)
    lm = I[1] + 1   # +1 because we leave out the lm=1 element
    
    k = ndims(div) == 1 ? CartesianIndex() : I[2]
    m = lm2m_indices[lm] 
    
    # div -= a-b+c where a = im*(m-1)*u, b = grad_y_vordiv1*v[lm-1], c = grad_y_vordiv2*v[lm+1]
    ∂u∂λ = im*(m-1)*u[lm, k]
    ∂v∂θ1 = grad_y_vordiv1[lm]*v[lm-1, k]
    ∂v∂θ2 = grad_y_vordiv2[lm]*v[lm+1, k]
    div[lm, k] -= ∂u∂λ - ∂v∂θ1 + ∂v∂θ2
end

# Legacy function for backward compatibility
function _divergence_KA!(  
    kernel,
    div::LowerTriangularArray,
    u::LowerTriangularArray,
    v::LowerTriangularArray,
    S::SpectralTransform;
    radius = DEFAULT_RADIUS,
)
    @warn "_divergence_KA! with a kernel function is deprecated, use the specialized versions instead"
    
    (; grad_y_vordiv1, grad_y_vordiv2, lm2m_indices) = S

    @boundscheck ismatching(S, div) || throw(DimensionMismatch(S, div))

    # first and last points are handled separately to avoid access to out-of-bounds array elements
    div[1,:] .= kernel.(div[1,:], 0, 0, grad_y_vordiv2[1].*v[2,:])
    
    launch!(S.architecture, :lmk_inner_points, size(div), _divergence_kernel!, kernel, div, u, v, grad_y_vordiv1, grad_y_vordiv2, lm2m_indices)

    # Last element
    div[end,:] .= 0

    # /radius scaling if not unit sphere
    if radius != 1
        div .*= inv(radius)
    end

    return div
end

@kernel inbounds=true function _divergence_kernel!(@Const(kernel_func), div, @Const(u), @Const(v), @Const(grad_y_vordiv1), @Const(grad_y_vordiv2), @Const(lm2m_indices))
    I = @index(Global, Cartesian)
    lm = I[1] + 1   # +1 because we leave out the lm=1 element
    
    k = ndims(div) == 1 ? CartesianIndex() : I[2]
    m = lm2m_indices[lm] 
    
    # This is the legacy kernel that uses a function closure
    div[lm, k] = kernel_func(div[lm, k], Complex(0,m-1)*u[lm, k], grad_y_vordiv1[lm]*v[lm-1, k], grad_y_vordiv2[lm]*v[lm+1, k])
end