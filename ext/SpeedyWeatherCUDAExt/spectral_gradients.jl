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

    dpdx .= @. complex(0,(lm2m_indices - 1))*p

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
    
function _divergence_KA!(  
    kernel,
    div::LowerTriangularArray,
    u::LowerTriangularArray,
    v::LowerTriangularArray,
    S::SpectralTransform;
    radius = DEFAULT_RADIUS,
)
    (; grad_y_vordiv1, grad_y_vordiv2, lm2ij_indices, lm2m_indices) = S
  
    @boundscheck ismatching(S, div) || throw(DimensionMismatch(S, div))

    ∂u∂λ = @. complex(0,(lm2m_indices - 1))*u

    ∂v∂θ1 = similar(div)
    ∂v∂θ2 = similar(div)

    launch!(S.architecture, :lmk, size(∂v∂θ1), ∂v∂θ1_kernel!, ∂v∂θ1, u, v, grad_y_vordiv1, grad_y_vordiv2, lm2ij_indices)


    kernel.(div, ∂u∂λ, ∂v∂θ1, ∂v∂θ2)




    # /radius scaling if not unit sphere
    if radius != 1
        div .*= inv(radius)
    end

    return div
end

@kernel inline=true function ∂v∂θ1_kernel!(∂v∂θ1, @Const(u), @Const(v), @Const(grad_y_vordiv1), @Const(grad_y_vordiv2), @Const(lm2ij_indices))

    I = @index(Global, Cartesian)
    lm = I[1]

    # To-Do: not really ideal, but I don't know how to do it better right now (except for two completely different kernels)
    k = ndims(div) == 1 ? CartesianIndex() : I[2]

    l = lm2l_indices[lm]
    m = lm2m_indices[lm]
    grad_y_vordiv = = grad_y_vordiv1

    if l==∂v∂θ1.m # Last row, only vectors make use of the lmax+1 row, set to zero for scalars div, curl
        ∂v∂θ1[I] = 0
    else 
        # distinguish DIAGONAL (to avoid access to v[l-1, m])
        ∂v∂θ1 = l==m ? 0 : grad_y_vordiv1[lm] * v[lm-1, k] 

        ∂v∂θ2 = grad_y_vordiv2[lm] * v[lm+1, k]  
        div[I] = kernel_func(div[I], ∂u∂λ[I], ∂v∂θ1, ∂v∂θ2)
    end 
end 


function _divergence_KA_old!(  
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

@kernel inline=true function _divergence_kernel!(kernel_func, div, u, v, grad_y_vordiv1, grad_y_vordiv2, @Const(lm2ij_indices))

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