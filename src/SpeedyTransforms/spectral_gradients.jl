const DEFAULT_RADIUS = 1

"""
$(TYPEDSIGNATURES)
Curl of a vector `u, v` written into `curl`, `curl = ∇×(u, v)`.
`u, v` are expected to have a 1/coslat-scaling included, otherwise `curl` is scaled.
Acts on the unit sphere, i.e. it omits 1/radius scaling as all gradient operators
unless the `radius` keyword argument is provided. `flipsign` option calculates -∇×(u, v) instead.
`add` option calculates `curl += ∇×(u, v)` instead. `flipsign` and `add` can be combined.
This functions only creates the kernel and calls the generic divergence function _divergence!
subsequently with flipped u, v -> v, u for the curl."""
function curl!(
    curl::LowerTriangularArray,
    u::LowerTriangularArray,
    v::LowerTriangularArray,
    S::SpectralTransform;
    flipsign::Bool=false,
    add::Bool=false,
    kwargs...,
)
    # = -(∂λ - ∂θ) or (∂λ - ∂θ), adding or overwriting the output curl 
    kernel = flipsign ? (add ? (o, a, b, c) -> o-(a+b-c) : (o, a, b, c) -> -(a+b-c)) :
                        (add ? (o, a, b, c) -> o+(a+b-c) : (o, a, b, c) -> a+b-c)  
    _divergence!(kernel, curl, v, u, S; kwargs...)      # flip u, v -> v, u
end

"""
$(TYPEDSIGNATURES)
Divergence of a vector `u, v` written into `div`, `div = ∇⋅(u, v)`. 
`u, v` are expected to have a 1/coslat-scaling included, otherwise `div` is scaled.
Acts on the unit sphere, i.e. it omits 1/radius scaling as all gradient operators,
unless the `radius` keyword argument is provided. `flipsign` option calculates -∇⋅(u, v) instead.
`add` option calculates `div += ∇⋅(u, v)` instead. `flipsign` and `add` can be combined.
This functions only creates the kernel and calls the generic divergence function _divergence! subsequently."""
function divergence!(
    div::LowerTriangularArray,
    u::LowerTriangularArray,
    v::LowerTriangularArray,
    S::SpectralTransform;
    flipsign::Bool=false,
    add::Bool=false,
    kwargs...,
)
    # = -(∂λ + ∂θ) or (∂λ + ∂θ), adding or overwriting the output div
    kernel = flipsign ? (add ? (o, a, b, c) -> o-(a-b+c) : (o, a, b, c) -> -(a-b+c)) :
                        (add ? (o, a, b, c) ->  o+(a-b+c) : (o, a, b, c) -> a-b+c)    
    _divergence!(kernel, div, u, v, S; kwargs...)
end

"""
$(TYPEDSIGNATURES)
Generic divergence function of vector `u`, `v` that writes into the output into `div`.
Generic as it uses the kernel `kernel` such that curl, div, add or flipsign
options are provided through `kernel`, but otherwise a single function is used.
Acts on the unit sphere, i.e. it omits 1/radius scaling as all gradient operators,
unless the `radius` keyword argument is provided."""
function _divergence!(  
    kernel,
    div::LowerTriangularArray,
    u::LowerTriangularArray,
    v::LowerTriangularArray,
    S::SpectralTransform;
    radius = DEFAULT_RADIUS,
)
    (; grad_y_vordiv1, grad_y_vordiv2 ) = S
  
    @boundscheck ismatching(S, div) || throw(DimensionMismatch(S, div))
    lmax, mmax = size(div, OneBased, as=Matrix)

    for k in eachmatrix(div, u, v)          # also checks size compatibility
        lm = 0
        @inbounds for m in 1:mmax           # 1-based l, m
            
            # DIAGONAL (separate to avoid access to v[l-1, m])
            lm += 1                                 
            ∂u∂λ  = ((m-1)*im)*u[lm, k]
            ∂v∂θ1 = 0                # always above the diagonal
            ∂v∂θ2 = grad_y_vordiv2[lm] * v[lm+1, k]
            div[lm, k] = kernel(div[lm, k], ∂u∂λ, ∂v∂θ1, ∂v∂θ2)
            
            # BELOW DIAGONAL (but skip last row)
            for l in m+1:lmax-1
                lm += 1
                ∂u∂λ  = ((m-1)*im)*u[lm, k]
                ∂v∂θ1 = grad_y_vordiv1[lm] * v[lm-1, k]
                ∂v∂θ2 = grad_y_vordiv2[lm] * v[lm+1, k]  # this pulls in data from the last row though
                div[lm, k] = kernel(div[lm, k], ∂u∂λ, ∂v∂θ1, ∂v∂θ2)
            end

            # Last row, only vectors make use of the lmax+1 row, set to zero for scalars div, curl
            lm += 1
            div[lm, k] = 0
        end
    end

    # /radius scaling if not unit sphere
    if radius != 1
        div .*= inv(radius)
    end

    return div
end

# KA / GPU version, called by the SpeedyWeatherCUDAExt when input with CuArrays
function _divergence_KA!(  
    kernel,
    div::LowerTriangularArray,
    u::LowerTriangularArray,
    v::LowerTriangularArray,
    S::SpectralTransform;
    radius = DEFAULT_RADIUS,
)
    (; grad_y_vordiv1, grad_y_vordiv2, lm2m_indices) = S
  
    @boundscheck ismatching(S, div) || throw(DimensionMismatch(S, div))

    # first element is special because it doesn't can't access v[0]
    div[1, :] .= kernel.(div[1, :], 0, 0, grad_y_vordiv2[1]*v[2, :])

    launch!(S.architecture, :lmk_inner_points, size(div), _divergence_kernel!, kernel, div, u, v, grad_y_vordiv1, grad_y_vordiv2, lm2m_indices)

    # last element is always zero (and we can access v[lmax+1])
    div[end,:] .= 0

    # /radius scaling if not unit sphere
    if radius != 1
        div .*= inv(radius)
    end

    return div
end

@kernel function _divergence_kernel!(kernel_func, div, u, v, grad_y_vordiv1, grad_y_vordiv2, @Const(lm2m_indices))

    I = @index(Global, Cartesian)
    lm = I[1] + 1

    k = ndims(div) == 1 ? CartesianIndex() : I[2]

    m = lm2m_indices[lm]

    ∂u∂λ  = ((m-1)*im)*u[I]
    ∂v∂θ1 = grad_y_vordiv1[lm] * v[lm-1, k] 
    ∂v∂θ2 = grad_y_vordiv2[lm] * v[lm+1, k]  
    div[I] = kernel_func(div[I], ∂u∂λ, ∂v∂θ1, ∂v∂θ2)
end 

"""
$(TYPEDSIGNATURES)
Divergence (∇⋅) of two vector components `u, v` which need to have size (n+1)xn,
the last row will be set to zero in the returned `LowerTriangularMatrix`.
This function requires both `u, v` to be transforms of fields that are scaled with
`1/cos(lat)`. Acts on the unit sphere, i.e. it omits 1/radius scaling unless
`radius` keyword argument is provided.
An example usage is therefore

    RingGrids.scale_coslat⁻¹!(u_grid)
    RingGrids.scale_coslat⁻¹!(v_grid)
    u = transform(u_grid, one_more_degree=true)
    v = transform(v_grid, one_more_degree=true)
    div = divergence(u, v, radius = 6.371e6)
    div_grid = transform(div)
"""
function divergence(u::LowerTriangularArray,
                    v::LowerTriangularArray;
                    kwargs...)

    S = SpectralTransform(u)
    return divergence(u, v, S; kwargs...)
end

# use SpectralTransform if provided
function divergence(u::LowerTriangularArray,
                    v::LowerTriangularArray,
                    S::SpectralTransform;
                    kwargs...)
    div = similar(u)
    return divergence!(div, u, v, S; add=false, flipsign=false, kwargs...)
end

# called by divergence or curl
function _div_or_curl(  
    kernel!,
    u::Grid,
    v::Grid;
    kwargs...,
) where {Grid<:AbstractGridArray}

    u_grid = copy(u)
    v_grid = copy(v)

    RingGrids.scale_coslat⁻¹!(u_grid)
    RingGrids.scale_coslat⁻¹!(v_grid)

    S = SpectralTransform(u_grid, one_more_degree=true)
    us = transform(u_grid, S)
    vs = transform(v_grid, S)

    div_or_vor = similar(us)
    kernel!(div_or_vor, us, vs, S; add=false, flipsign=false, kwargs...)
    return div_or_vor
end

"""
$(TYPEDSIGNATURES)
Divergence (∇⋅) of two vector components `u, v` on a grid.
Applies 1/coslat scaling, transforms to spectral space and returns
the spectral divergence. Acts on the unit sphere, i.e. it omits 1/radius scaling unless
`radius` keyword argument is provided.
"""
divergence(u::Grid, v::Grid; kwargs...) where {Grid<:AbstractGridArray} =
    _div_or_curl(divergence!, u, v; kwargs...)

"""
$(TYPEDSIGNATURES)
Curl (∇×) of two vector components `u, v` on a grid.
Applies 1/coslat scaling, transforms to spectral space and returns
the spectral curl. Acts on the unit sphere, i.e. it omits 1/radius scaling unless
`radius` keyword argument is provided."""
curl(u::Grid, v::Grid; kwargs...) where {Grid<:AbstractGridArray} =
    _div_or_curl(curl!, u, v; kwargs...)

"""
$(TYPEDSIGNATURES)
Curl (∇×) of two vector components `u, v` of size (n+1)xn, the last row
will be set to zero in the returned `LowerTriangularMatrix`. This function
requires both `u, v` to be transforms of fields that are scaled with
`1/cos(lat)`. Acts on the unit sphere, i.e. it omits 1/radius scaling unless
`radius` keyword argument is provided. An example usage is therefore

    RingGrids.scale_coslat⁻¹!(u_grid)
    RingGrids.scale_coslat⁻¹!(v_grid)
    u = transform(u_grid)
    v = transform(v_grid)
    vor = curl(u, v, radius=6.371e6)
    vor_grid = transform(div)
"""
function curl(  u::LowerTriangularArray,
                v::LowerTriangularArray;
                kwargs...)

    S = SpectralTransform(u)
    return curl(u, v, S; kwargs...)
end

# use SpectralTransform if provided
function curl(  u::LowerTriangularArray,
                v::LowerTriangularArray,
                S::SpectralTransform;
                kwargs...)
    vor = similar(u)
    curl!(vor, u, v, S; add=false, flipsign=false, kwargs...)
    return vor
end

"""
$(TYPEDSIGNATURES)
Get U, V (=(u, v)*coslat) from vorticity ζ spectral space (divergence D=0)
Two operations are combined into a single linear operation. First, invert the
spherical Laplace ∇² operator to get stream function from vorticity. Then
compute zonal and meridional gradients to get U, V.
Acts on the unit sphere, i.e. it omits any radius scaling as all inplace gradient operators,
unless the `radius` keyword argument is provided."""
function UV_from_vor!(  
    U::LowerTriangularArray,
    V::LowerTriangularArray,
    vor::LowerTriangularArray,
    S::SpectralTransform;
    radius = DEFAULT_RADIUS,
)
    (; vordiv_to_uv_x, vordiv_to_uv1, vordiv_to_uv2 ) = S
    @boundscheck ismatching(S, U) || throw(DimensionMismatch(S, U))
    
    # maximum degree l, order m of spherical harmonics (1-based)
    lmax, mmax = size(U, OneBased, as=Matrix)

    for k in eachmatrix(U, V, vor)                          # also checks size compatibility
        lm = 0
        @inbounds for m in 1:mmax-1                         # 1-based l, m, exclude last column

            # DIAGONAL (separated to avoid access to l-1, m which is above the diagonal)
            lm += 1

            # U = -∂/∂lat(Ψ) and V = V = ∂/∂λ(Ψ) combined with Laplace inversion ∇⁻², omit radius R scaling
            U[lm, k] = vordiv_to_uv2[lm] * vor[lm+1, k]     # - vordiv_to_uv1[lm]*vor[l-1, m] <- is zero
            V[lm, k] = im*vordiv_to_uv_x[lm] * vor[lm, k]

            # BELOW DIAGONAL
            for l in m+1:lmax-2                             # skip last two rows
                lm += 1

                # U = -∂/∂lat(Ψ) and V = V = ∂/∂λ(Ψ) combined with Laplace inversion ∇⁻², omit radius R scaling
                # U[lm] = vordiv_to_uv2[lm]*vor[lm+1] - vordiv_to_uv1[lm]*vor[lm-1]
                U[lm, k] = muladd(vordiv_to_uv2[lm], vor[lm+1, k], -vordiv_to_uv1[lm]*vor[lm-1, k])
                V[lm, k] = im*vordiv_to_uv_x[lm] * vor[lm, k]
            end

            # SECOND LAST ROW
            lm += 1
            U[lm, k] = -vordiv_to_uv1[lm] * vor[lm-1, k]    # meridional gradient again (but only 2nd term from above)
            V[lm, k] = im*vordiv_to_uv_x[lm] * vor[lm, k]   # zonal gradient again (as above)

            # LAST ROW (separated to avoid out-of-bounds access to l+2, m)
            lm += 1
            U[lm, k] = -vordiv_to_uv1[lm] * vor[lm-1, k]    # meridional gradient again (but only 2nd term from above)
            V[lm, k] = 0                                    # set explicitly to 0 as Ψ does not contribute to last row of V
        end

        # LAST COLUMN
        @inbounds begin
            lm += 1                     # second last row
            U[lm, k] = 0
            V[lm, k] = im*vordiv_to_uv_x[lm] * vor[lm, k]

            lm += 1                     # last row
            U[lm, k] = -vordiv_to_uv1[lm] * vor[lm-1, k]
            V[lm, k] = 0
        end
    end

    # *radius scaling if not unit sphere (*radius² for ∇⁻² then /radius to get from stream function to velocity)
    if radius != 1
        U .*= radius
        V .*= radius
    end

    return U, V
end

"""
$(TYPEDSIGNATURES)
Get U, V (=(u, v)*coslat) from vorticity ζ and divergence D in spectral space.
Two operations are combined into a single linear operation. First, invert the
spherical Laplace ∇² operator to get stream function from vorticity and
velocity potential from divergence. Then compute zonal and meridional gradients
to get U, V.
Acts on the unit sphere, i.e. it omits any radius scaling as all inplace gradient operators.
"""
function UV_from_vordiv!(   
    U::LowerTriangularArray,
    V::LowerTriangularArray,
    vor::LowerTriangularArray,
    div::LowerTriangularArray,
    S::SpectralTransform;
    radius = DEFAULT_RADIUS,
)
    (; vordiv_to_uv_x, vordiv_to_uv1, vordiv_to_uv2 ) = S
    @boundscheck ismatching(S, U) || throw(DimensionMismatch(S, U))

    # maximum degree l, order m of spherical harmonics (1-based)
    lmax, mmax = size(U, OneBased, as=Matrix)

    for k in eachmatrix(U, V, vor, div)                 # also checks size compatibility
        lm = 0
        @inbounds for m in 1:mmax-1                     # 1-based l, m, skip last column

            # DIAGONAL (separated to avoid access to l-1, m which is above the diagonal)
            lm += 1
            
            # div, vor contribution to meridional gradient
            ∂ζθ =  vordiv_to_uv2[lm]*vor[lm+1, k]       # lm-1 term is zero
            ∂Dθ = -vordiv_to_uv2[lm]*div[lm+1, k]       # lm-1 term is zero
            
            # the following is moved into the muladd        
            # ∂Dλ = im*vordiv_to_uv_x[lm]*div[lm]       # divergence contribution to zonal gradient
            # ∂ζλ = im*vordiv_to_uv_x[lm]*vor[lm]       # vorticity contribution to zonal gradient

            z = im*vordiv_to_uv_x[lm]
            U[lm, k] = muladd(z, div[lm, k], ∂ζθ)       # = ∂Dλ + ∂ζθ
            V[lm, k] = muladd(z, vor[lm, k], ∂Dθ)       # = ∂ζλ + ∂Dθ

            # BELOW DIAGONAL (all terms)
            for l in m+1:lmax-2                         # skip last two rows (lmax-1, lmax)
                lm += 1
                
                # div, vor contribution to meridional gradient
                # ∂ζθ = vordiv_to_uv2[lm]*vor[lm+1] - vordiv_to_uv1[lm]*vor[lm-1]
                # ∂Dθ = vordiv_to_uv1[lm]*div[lm-1] - vordiv_to_uv2[lm]*div[lm+1]
                ∂ζθ = muladd(vordiv_to_uv2[lm], vor[lm+1, k], -vordiv_to_uv1[lm]*vor[lm-1, k])
                ∂Dθ = muladd(vordiv_to_uv1[lm], div[lm-1, k], -vordiv_to_uv2[lm]*div[lm+1, k])

                # The following is moved into the muladd
                # ∂Dλ = im*vordiv_to_uv_x[lm]*div[lm]   # divergence contribution to zonal gradient
                # ∂ζλ = im*vordiv_to_uv_x[lm]*vor[lm]   # vorticity contribution to zonal gradient

                z = im*vordiv_to_uv_x[lm]
                U[lm, k] = muladd(z, div[lm, k], ∂ζθ)   # = ∂Dλ + ∂ζθ
                V[lm, k] = muladd(z, vor[lm, k], ∂Dθ)   # = ∂ζλ + ∂Dθ            
            end

            # SECOND LAST ROW (separated to imply that vor, div are zero in last row)
            lm += 1
            U[lm, k] = im*vordiv_to_uv_x[lm]*div[lm, k] - vordiv_to_uv1[lm]*vor[lm-1, k]
            V[lm, k] = im*vordiv_to_uv_x[lm]*vor[lm, k] + vordiv_to_uv1[lm]*div[lm-1, k]

            # LAST ROW (separated to avoid out-of-bounds access to lmax+1)
            lm += 1
            U[lm, k] = -vordiv_to_uv1[lm]*vor[lm-1, k]  # only last term from 2nd last row
            V[lm, k] =  vordiv_to_uv1[lm]*div[lm-1, k]  # only last term from 2nd last row
        end

        # LAST COLUMN
        @inbounds begin
            lm += 1                                         # second last row
            U[lm, k] = im*vordiv_to_uv_x[lm]*div[lm, k]     # other terms are zero
            V[lm, k] = im*vordiv_to_uv_x[lm]*vor[lm, k]     # other terms are zero

            lm += 1                                         # last row
            U[lm, k] = -vordiv_to_uv1[lm]*vor[lm-1, k]      # other terms are zero
            V[lm, k] =  vordiv_to_uv1[lm]*div[lm-1, k]      # other terms are zero
        end
    end

    # *radius scaling if not unit sphere (*radius² for ∇⁻², then /radius to get from stream function to velocity)
    if radius != 1
        U .*= radius
        V .*= radius
    end

    return U, V
end

"""
$(TYPEDSIGNATURES)
Laplace operator ∇² applied to the spectral coefficients `alms` in spherical
coordinates. The eigenvalues which are precomputed in `S`.
∇²! is the in-place version which directly stores the output in the first argument `∇²alms`.
Acts on the unit sphere, i.e. it omits any radius scaling as all inplace gradient operators,
unless the `radius` keyword argument is provided.

Keyword arguments
=================

  - `add=true` adds the ∇²(alms) to the output
  - `flipsign=true` computes -∇²(alms) instead
  - `inverse=true` computes ∇⁻²(alms) instead

Default is `add=false`, `flipsign=false`, `inverse=false`. These options can be combined."""
function ∇²!(
    ∇²alms::LowerTriangularArray,   # Output: (inverse) Laplacian of alms
    alms::LowerTriangularArray,     # Input: spectral coefficients
    S::SpectralTransform;           # precomputed eigenvalues
    add::Bool=false,                # add to output array or overwrite
    flipsign::Bool=false,           # -∇² or ∇²
    inverse::Bool=false,            # ∇⁻² or ∇²
    radius = DEFAULT_RADIUS,        # scale with radius if provided, otherwise unit sphere
)
    @boundscheck ismatching(S, ∇²alms) || throw(DimensionMismatch(S, ∇²alms))

    # use eigenvalues⁻¹/eigenvalues for ∇⁻²/∇² based but name both eigenvalues
    eigenvalues = inverse ? S.eigenvalues⁻¹ : S.eigenvalues

    kernel = flipsign ? (add ? (o,a) -> (o-a) : (o, a) -> -a) : 
                        (add ? (o,a) -> (o+a) : (o, a) -> a)
    
    launch!(S.architecture, :lmk, size(∇²alms), ∇²_kernel!, ∇²alms, alms, eigenvalues, kernel, S.lm2l_indices)

    # /radius² or *radius² scaling if not unit sphere
    if radius != 1
        R_plusminus_squared = inverse ? radius^2 : inv(radius^2)
        ∇²alms .*= R_plusminus_squared
    end

    return ∇²alms
end

@kernel function ∇²_kernel!(∇²alms, alms, @Const(eigenvalues), kernel_func, @Const(lm2l_indices))

    I = @index(Global, Cartesian) # I[1] == lm, I[2] == k
                                  # we use cartesian index instead of NTuple here
                                  # because this works for 2D and 3D matrices
    l = lm2l_indices[I[1]]

    ∇²alms[I] = kernel_func(∇²alms[I], alms[I]*eigenvalues[l])
end 

"""
$(TYPEDSIGNATURES)
Laplace operator ∇² applied to input `alms`, using precomputed eigenvalues from `S`.
Acts on the unit sphere, i.e. it omits 1/radius^2 scaling unless
`radius` keyword argument is provided."""
function ∇²(
    alms::LowerTriangularArray,     # Input: spectral coefficients
    S::SpectralTransform;           # precomputed eigenvalues
    kwargs...,
)
    ∇²alms = similar(alms)
    ∇²!(∇²alms, alms, S; add=false, flipsign=false, inverse=false, kwargs...)
    return ∇²alms
end

"""
$(TYPEDSIGNATURES)
Returns the Laplace operator ∇² applied to input `alms`.
Acts on the unit sphere, i.e. it omits 1/radius^2 scaling unless
`radius` keyword argument is provided."""
∇²(alms::LowerTriangularArray; kwargs...) = ∇²(alms, SpectralTransform(alms); kwargs...)

"""
$(TYPEDSIGNATURES)
InverseLaplace operator ∇⁻² applied to input `alms`, using precomputed
eigenvalues from `S`. Acts on the unit sphere, i.e. it omits radius^2 scaling unless
`radius` keyword argument is provided."""
function ∇⁻²(
    ∇²alms::LowerTriangularArray,   # Input: spectral coefficients
    S::SpectralTransform;           # precomputed eigenvalues
    kwargs...,
)
    alms = similar(∇²alms)
    ∇⁻²!(alms, ∇²alms, S; add=false, flipsign=false, kwargs...)
    return alms
end

"""
$(TYPEDSIGNATURES)
Returns the inverse Laplace operator ∇⁻² applied to input `alms`.
Acts on the unit sphere, i.e. it omits radius^2 scaling unless
`radius` keyword argument is provided."""
∇⁻²(∇²alms::LowerTriangularArray; kwargs...) = ∇⁻²(∇²alms, SpectralTransform(∇²alms); kwargs...)

"""$(TYPEDSIGNATURES) Calls `∇²!(∇⁻²alms, alms, S; add, flipsign, inverse=true)`."""
function ∇⁻²!(
    ∇⁻²alms::LowerTriangularArray,  # Output: inverse Laplacian of alms
    alms::LowerTriangularArray,     # Input: spectral coefficients
    S::SpectralTransform;           # precomputed eigenvalues
    add::Bool = false,              # add to output array or overwrite
    flipsign::Bool = false,         # -∇⁻² or ∇⁻²
    kwargs...,
)
    inverse = true
    return ∇²!(∇⁻²alms, alms, S; add, flipsign, inverse, kwargs...)
end

"""$(TYPEDSIGNATURES) Applies the gradient operator ∇ applied to input `p` and stores the result
in `dpdx` (zonal derivative) and `dpdy` (meridional derivative). The gradient operator acts
on the unit sphere and therefore omits the 1/radius scaling unless `radius` keyword argument is provided."""
function ∇!(
    dpdx::LowerTriangularArray,     # Output: zonal gradient
    dpdy::LowerTriangularArray,     # Output: meridional gradient
    p::LowerTriangularArray,        # Input: spectral coefficients
    S::SpectralTransform;           # includes precomputed arrays
    radius = DEFAULT_RADIUS,        # scale with radius if provided, otherwise unit sphere
)
    (; grad_y1, grad_y2, lm2m_indices) = S
    @boundscheck ismatching(S, p) || throw(DimensionMismatch(S, p))

    @. dpdx = complex(0, lm2m_indices - 1)*p

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

@kernel inbounds=true function dpdy_kernel!(dpdy, @Const(p), @Const(grad_y1), @Const(grad_y2))
    I = @index(Global, Cartesian)
    lm = I[1] + 1   # +1 because we leave out the lm=1 element
    k = ndims(p) == 1 ? CartesianIndex() : I[2]

    gy1 = grad_y1[lm]
    gy2 = grad_y2[lm]
    
    # compared to the old CPU only version, some of the gy1 and gy2 are zero 
    # that's why we don't need to check for l==m (diagonal) or l==p.m (last row)
    dpdy[lm, k] = gy1*p[lm-1, k] + gy2*p[lm+1, k]
end

"""$(TYPEDSIGNATURES) The zonal and meridional gradient of `p`
using an existing `SpectralTransform` `S`. Acts on the unit sphere,
i.e. it omits 1/radius scaling unless `radius` keyword argument is provided."""
function ∇(p::LowerTriangularArray, S::SpectralTransform; kwargs...)
    dpdx = similar(p)
    dpdy = similar(p)
    ∇!(dpdx, dpdy, p, S; kwargs...)
    return dpdx, dpdy
end

"""$(TYPEDSIGNATURES) The zonal and meridional gradient of `p`.
Precomputes a `SpectralTransform` `S`. Acts on the unit-sphere,
i.e. it omits 1/radius scaling unless `radius` keyword argument is provided."""
function ∇(p::LowerTriangularArray; kwargs...)
    S = SpectralTransform(p, one_more_degree=true)
    return ∇(p, S; kwargs...)
end

"""$(TYPEDSIGNATURES) The zonal and meridional gradient of `grid`.
Transform to spectral space, takes the gradient and unscales the 1/coslat
scaling in the gradient. Acts on the unit-sphere, i.e. it omits 1/radius scaling unless
`radius` keyword argument is provided. Makes use of an existing spectral transform `S`."""
function ∇(grid::AbstractGridArray, S::SpectralTransform; kwargs...)
    p = transform(grid, S)
    dpdx, dpdy = ∇(p, S; kwargs...)
    dpdx_grid = transform(dpdx, S, unscale_coslat=true)
    dpdy_grid = transform(dpdy, S, unscale_coslat=true)
    return dpdx_grid, dpdy_grid
end

"""$(TYPEDSIGNATURES) The zonal and meridional gradient of `grid`.
Transform to spectral space, takes the gradient and unscales the 1/coslat
scaling in the gradient. Acts on the unit-sphere, i.e. it omits 1/radius scaling unless
`radius` keyword argument is provided."""
function ∇(grid::AbstractGridArray; kwargs...)
    S = SpectralTransform(grid, one_more_degree=true)
    return ∇(grid, S; kwargs...)
end
