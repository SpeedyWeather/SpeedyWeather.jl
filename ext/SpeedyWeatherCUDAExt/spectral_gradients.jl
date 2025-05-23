import SpeedyWeather.SpectralTransforms: ismatching, _divergence!

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
    (; grad_y1, grad_y2) = S
    (; m_indices) = p.spectrum 

    @boundscheck ismatching(S, p) || throw(DimensionMismatch(S, p))

    dpdx .= @. complex(0, m_indices - 1)*p

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
SpeedyWeather.SpectralTransforms._divergence!(
    kernel,
    div::LowerTriangularArray{T,N,<:CuArray},
    u::LowerTriangularArray{T,N,<:CuArray},
    v::LowerTriangularArray{T,N,<:CuArray},
    S::SpectralTransform;
    radius = DEFAULT_RADIUS,
) = _divergence_KA!(kernel, div, u, v, S; radius)
    
