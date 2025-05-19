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
    (; grad_y1, grad_y2, lm2l_indices, lm2m_indices) = S
    @boundscheck ismatching(S, p) || throw(DimensionMismatch(S, p))

    dpdx .= @. complex(0,(lm2m_indices - 1))*p
    launch!(S.architecture, :lmk, size(dpdy), dpdy_kernel!, dpdy, p, grad_y1, grad_y2, lm2l_indices, lm2m_indices)

    # 1/radius factor if not unit sphere
    if radius != 1
        R⁻¹ = inv(radius)
        dpdx .*= R⁻¹
        dpdy .*= R⁻¹
    end

    return dpdx, dpdy
end

@kernel inbounds=true function dpdy_kernel!(dpdy, p, @Const(grad_y1), @Const(grad_y2), @Const(lm2l_indices), @Const(lm2m_indices))
    I = @index(Global, Cartesian)
    lm = I[1]
    k = ndims(p) == 1 ? CartesianIndex() : I[2]

    l = lm2l_indices[lm]
    m = lm2m_indices[lm]
    gy1 = grad_y1[lm]
    gy2 = grad_y2[lm]
    
    if l == m
        # DIAGONAL (separated to avoid access to l-1, m which is above the diagonal)
        # meridional gradient: p[lm-1]=0 on diagonal
        dpdy[I] = gy2*p[lm+1, k]
    elseif l == p.m
        # LAST ROW (separated to avoid out-of-bounds access to lmax+1)
        dpdy[I] = gy1*p[lm-1, k]
    else
        # all other cases
        dpdy[I] = gy1*p[lm-1, k] + gy2*p[lm+1, k]
    end
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
    
function _divergence_KA!(  
    kernel,
    div::LowerTriangularArray,
    u::LowerTriangularArray,
    v::LowerTriangularArray,
    S::SpectralTransform;
    radius = DEFAULT_RADIUS,
)
    (; grad_y_vordiv1, grad_y_vordiv2, lm2ij_indices) = S
  
    @boundscheck ismatching(S, div) || throw(DimensionMismatch(S, div))

    launch!(S.architecture, :lmk, size(div), _divergence_kernel!, kernel, div, u, v, grad_y_vordiv1, grad_y_vordiv2, lm2ij_indices)

    # /radius scaling if not unit sphere
    if radius != 1
        div .*= inv(radius)
    end

    return div
end

@kernel function _divergence_kernel!(kernel_func, div, u, v, grad_y_vordiv1, grad_y_vordiv2, @Const(lm2ij_indices))

    I = @index(Global, Cartesian)
    lm = I[1]

    # To-Do: not really ideal, but I don't know how to do it better right now (except for two completely different kernels)
    k = ndims(div) == 1 ? CartesianIndex() : I[2]

    l = lm2ij_indices[lm, 1]
    m = lm2ij_indices[lm, 2]

    if l==size(div, 1, as=Matrix) # Last row, only vectors make use of the lmax+1 row, set to zero for scalars div, curl
        div[I] = 0
    else 
        ∂u∂λ  = ((m-1)*im)*u[I]

        # distinguish DIAGONAL (to avoid access to v[l-1, m])
        ∂v∂θ1 = l==m ? 0 : grad_y_vordiv1[lm] * v[lm-1, k] 

        ∂v∂θ2 = grad_y_vordiv2[lm] * v[lm+1, k]  
        div[I] = kernel_func(div[I], ∂u∂λ, ∂v∂θ1, ∂v∂θ2)
    end 
end 