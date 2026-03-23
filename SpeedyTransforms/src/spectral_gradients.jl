const DEFAULT_RADIUS = 1

"""
    KernelOP{mode, flipsign, add}

Type for dispatching on kernel operations in spectral gradient calculations.
- `mode`: `true` for curl, `false` for divergence
- `flipsign`: `true` or `false` to negate the result
- `add`: `true` or `false` to add to the output instead of overwriting
"""
struct KernelOP{mode, flipsign, add} end

# Curl operations (mode == true): a+b-c
# Standard curl (no flipsign, no add)
@inline (::KernelOP{true, false, false})(o, a, b, c) = a + b - c

# Curl with flipsign (no add)
@inline (::KernelOP{true, true, false})(o, a, b, c) = -(a + b - c)

# Curl with add (no flipsign)
@inline (::KernelOP{true, false, true})(o, a, b, c) = o + (a + b - c)

# Curl with flipsign and add
@inline (::KernelOP{true, true, true})(o, a, b, c) = o - (a + b - c)

# Divergence operations (mode == false): a-b+c
# Standard divergence (no flipsign, no add)
@inline (::KernelOP{false, false, false})(o, a, b, c) = a - b + c

# Divergence with flipsign (no add)
@inline (::KernelOP{false, true, false})(o, a, b, c) = -(a - b + c)

# Divergence with add (no flipsign)
@inline (::KernelOP{false, false, true})(o, a, b, c) = o + (a - b + c)

# Divergence with flipsign and add
@inline (::KernelOP{false, true, true})(o, a, b, c) = o - (a - b + c)


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
        S::AbstractSpectralTransform;
        flipsign::Bool = false,
        add::Bool = false,
        kwargs...,
    )
    # = -(∂λ - ∂θ) or (∂λ - ∂θ), adding or overwriting the output curl
    kernel = KernelOP{true, flipsign, add}()
    return _divergence!(kernel, curl, v, u, S; kwargs...)      # flip u, v -> v, u
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
        S::AbstractSpectralTransform;
        flipsign::Bool = false,
        add::Bool = false,
        kwargs...,
    )
    # = -(∂λ + ∂θ) or (∂λ + ∂θ), adding or overwriting the output div
    kernel = KernelOP{false, flipsign, add}()
    return _divergence!(kernel, div, u, v, S; kwargs...)
end

function _divergence!(
        kernel,
        div::LowerTriangularArray,
        u::LowerTriangularArray,
        v::LowerTriangularArray,
        S::AbstractSpectralTransform;
        radius = DEFAULT_RADIUS,
    )
    (; grad_x_vordiv, grad_y_vordiv1, grad_y_vordiv2) = S.gradients

    @boundscheck ismatching(S, div) || throw(DimensionMismatch(S, div))

    launch!(architecture(div), SpectralWorkOrder, size(div), _divergence_kernel!, kernel, div, u, v, grad_x_vordiv, grad_y_vordiv1, grad_y_vordiv2)

    # radius scaling if not unit sphere
    if radius != 1
        div .*= inv(radius)
    end

    return div
end

@kernel inbounds = true function _divergence_kernel!(kernel_func::KernelOP{mode, flipsign, add}, div, u, v, grad_x_vordiv, grad_y_vordiv1, grad_y_vordiv2) where {mode, flipsign, add}

    I = @index(Global, Cartesian)
    lm = I[1]
    lmmax = size(div, 1)
    k = ndims(div) == 1 ? CartesianIndex() : I[2]

    if lm == 1
        div[I] = kernel_func(div[I], 0, 0, grad_y_vordiv2[1] * v[2, k])
    elseif lm == lmmax
        div[I] = 0
    else
        ∂u∂λ = im * grad_x_vordiv[lm] * u[I]
        ∂v∂θ1 = grad_y_vordiv1[lm] * v[lm - 1, k]
        ∂v∂θ2 = grad_y_vordiv2[lm] * v[lm + 1, k]
        div[I] = kernel_func(div[I], ∂u∂λ, ∂v∂θ1, ∂v∂θ2)
    end
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
function divergence(
        u::LowerTriangularArray,
        v::LowerTriangularArray;
        kwargs...
    )

    S = SpectralTransform(u)
    return divergence(u, v, S; kwargs...)
end

# use SpectralTransform if provided
function divergence(
        u::LowerTriangularArray,
        v::LowerTriangularArray,
        S::AbstractSpectralTransform;
        kwargs...
    )
    div = similar(u)
    return divergence!(div, u, v, S; add = false, flipsign = false, kwargs...)
end

# called by divergence or curl
function _div_or_curl(
        kernel!,
        u::AbstractField,
        v::AbstractField;
        kwargs...,
    )
    u_grid = copy(u)
    v_grid = copy(v)

    RingGrids.scale_coslat⁻¹!(u_grid)
    RingGrids.scale_coslat⁻¹!(v_grid)

    S = SpectralTransform(u_grid, one_more_degree = true)
    us = transform(u_grid, S)
    vs = transform(v_grid, S)

    div_or_vor = similar(us)
    kernel!(div_or_vor, us, vs, S; add = false, flipsign = false, kwargs...)
    return div_or_vor
end

"""
$(TYPEDSIGNATURES)
Divergence (∇⋅) of two vector components `u, v` on a grid.
Applies 1/coslat scaling, transforms to spectral space and returns
the spectral divergence. Acts on the unit sphere, i.e. it omits 1/radius scaling unless
`radius` keyword argument is provided.
"""
divergence(u::AbstractField, v::AbstractField; kwargs...) = _div_or_curl(divergence!, u, v; kwargs...)

"""
$(TYPEDSIGNATURES)
Curl (∇×) of two vector components `u, v` on a grid.
Applies 1/coslat scaling, transforms to spectral space and returns
the spectral curl. Acts on the unit sphere, i.e. it omits 1/radius scaling unless
`radius` keyword argument is provided."""
curl(u::AbstractField, v::AbstractField; kwargs...) = _div_or_curl(curl!, u, v; kwargs...)

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
function curl(
        u::LowerTriangularArray,
        v::LowerTriangularArray;
        kwargs...
    )

    S = SpectralTransform(u)
    return curl(u, v, S; kwargs...)
end

# use SpectralTransform if provided
function curl(
        u::LowerTriangularArray,
        v::LowerTriangularArray,
        S::AbstractSpectralTransform;
        kwargs...
    )
    vor = similar(u)
    curl!(vor, u, v, S; add = false, flipsign = false, kwargs...)
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
        S::AbstractSpectralTransform;
        radius = DEFAULT_RADIUS,
    )
    (; vordiv_to_uv_x, vordiv_to_uv1, vordiv_to_uv2) = S.gradients
    @boundscheck ismatching(S, U) || throw(DimensionMismatch(S, U))

    launch!(architecture(U), SpectralWorkOrder, size(U), _UV_from_vor_kernel!, U, V, vor, vor.spectrum.l_indices, vor.spectrum.lmax, vordiv_to_uv_x, vordiv_to_uv1, vordiv_to_uv2)

    # *radius scaling if not unit sphere (*radius² for ∇⁻² then /radius to get from stream function to velocity)
    if radius != 1
        U .*= radius
        V .*= radius
    end

    return U, V
end

@kernel inbounds = true function _UV_from_vor_kernel!(U, V, vor, l_indices, lmax, vordiv_to_uv_x, vordiv_to_uv1, vordiv_to_uv2)
    I = @index(Global, Cartesian)
    lm = I[1]
    k = ndims(vor) == 1 ? CartesianIndex() : I[2]
    l = l_indices[lm]

    # Get the coefficients for the current lm index
    z = im * vordiv_to_uv_x[lm]
    vordiv_uv1 = vordiv_to_uv1[lm]
    vordiv_uv2 = vordiv_to_uv2[lm]

    # Handle different cases based on position in the triangular matrix
    if lm == 1  # First element (diagonal)
        # U = -∂/∂lat(Ψ) - no lm-1 term for first element
        U[I] = vordiv_uv2 * vor[lm + 1, k]
        # V = ∂/∂λ(Ψ)
        V[I] = z * vor[I]
    elseif l == (lmax - 1)     # extra in case vor[lmax,:] != 0, see comment in UV_from_vordiv!
        U[I] = -vordiv_uv1 * vor[lm - 1, k]    # meridional gradient again (but only 2nd term from above)
        V[I] = z * vor[I]          # zonal gradient again (as above)
    elseif l == lmax         # extra in case vor[lmax,:] != 0, see comment in UV_from_vordiv!
        U[I] = -vordiv_uv1 * vor[lm - 1, k]
        V[I] = 0
    else  # General case (below diagonal)
        # U = -∂/∂lat(Ψ) combined with Laplace inversion ∇⁻²
        U[I] = muladd(vordiv_uv2, vor[lm + 1, k], -vordiv_uv1 * vor[lm - 1, k])
        # V = ∂/∂λ(Ψ)
        V[I] = z * vor[I]
    end
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
    (; vordiv_to_uv_x, vordiv_to_uv1, vordiv_to_uv2) = S.gradients
    @boundscheck ismatching(S, U) || throw(DimensionMismatch(S, U))

    launch!(architecture(U), SpectralWorkOrder, size(U), _UV_from_vordiv_kernel!, U, V, vor, div, vor.spectrum.l_indices, vor.spectrum.lmax, vordiv_to_uv_x, vordiv_to_uv1, vordiv_to_uv2)

    # *radius scaling if not unit sphere (*radius² for ∇⁻², then /radius to get from stream function to velocity)
    if radius != 1
        U .*= radius
        V .*= radius
    end

    return U, V
end

@kernel inbounds = true function _UV_from_vordiv_kernel!(U, V, vor, div, l_indices, lmax, vordiv_to_uv_x, vordiv_to_uv1, vordiv_to_uv2)
    I = @index(Global, Cartesian)
    lm = I[1]
    k = ndims(vor) == 1 ? CartesianIndex() : I[2]
    l = l_indices[lm]

    # Get the coefficients for the current lm index
    z = im * vordiv_to_uv_x[lm]
    vordiv_uv1 = vordiv_to_uv1[lm]
    vordiv_uv2 = vordiv_to_uv2[lm]

    # Handle different cases based on position in the triangular matrix
    if lm == 1  # First element (diagonal)
        # Meridional gradient contributions
        ∂ζθ = vordiv_uv2 * vor[lm + 1, k]       # lm-1 term is zero
        ∂Dθ = -vordiv_uv2 * div[lm + 1, k]       # lm-1 term is zero

        # Zonal gradient contributions
        U[I] = muladd(z, div[I], ∂ζθ)          # = ∂Dλ + ∂ζθ
        V[I] = muladd(z, vor[I], ∂Dθ)          # = ∂ζλ + ∂Dθ
    elseif l == (lmax - 1)  # Second last row
        U[I] = muladd(z, div[I], -vordiv_uv1 * vor[lm - 1, k])
        V[I] = muladd(z, vor[I], vordiv_uv1 * div[lm - 1, k])
    elseif l == lmax  # Last row
        U[I] = -vordiv_uv1 * vor[lm - 1, k]      # only last term from 2nd last row
        V[I] = vordiv_uv1 * div[lm - 1, k]      # only last term from 2nd last row
    else  # General case (below diagonal)
        # Meridional gradient contributions
        ∂ζθ = muladd(vordiv_uv2, vor[lm + 1, k], -vordiv_uv1 * vor[lm - 1, k])
        ∂Dθ = muladd(vordiv_uv1, div[lm - 1, k], -vordiv_uv2 * div[lm + 1, k])

        # Zonal gradient contributions
        U[I] = muladd(z, div[I], ∂ζθ)          # = ∂Dλ + ∂ζθ
        V[I] = muladd(z, vor[I], ∂Dθ)          # = ∂ζλ + ∂Dθ
    end
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

Default is `add=false`, `flipsign=false`, `inverse=false`. These options can be combined.
"""
function ∇²!(
        ∇²alms::LowerTriangularArray,       # Output: (inverse) Laplacian of alms
        alms::LowerTriangularArray,         # Input: spectral coefficients
        S::AbstractSpectralTransform;       # precomputed eigenvalues
        add::Bool = false,                  # add to output array or overwrite
        flipsign::Bool = false,             # -∇² or ∇²
        inverse::Bool = false,              # ∇⁻² or ∇²
        radius = DEFAULT_RADIUS,            # scale with radius if provided, otherwise unit sphere
    )
    @boundscheck ismatching(S, ∇²alms) || throw(DimensionMismatch(S, ∇²alms))

    # use eigenvalues⁻¹/eigenvalues for ∇⁻²/∇² based but name both eigenvalues
    eigenvalues = inverse ? S.gradients.eigenvalues⁻¹ : S.gradients.eigenvalues

    kernel = flipsign ? (add ? (o, a) -> (o - a) : (o, a) -> -a) :
        (add ? (o, a) -> (o + a) : (o, a) -> a)

    launch!(architecture(∇²alms), SpectralWorkOrder, size(∇²alms), ∇²_kernel!, ∇²alms, alms, eigenvalues, kernel, alms.spectrum.l_indices)

    # /radius² or *radius² scaling if not unit sphere
    if radius != 1
        R_plusminus_squared = inverse ? radius^2 : inv(radius^2)
        ∇²alms .*= R_plusminus_squared
    end

    return ∇²alms
end

@kernel function ∇²_kernel!(∇²alms, alms, eigenvalues, kernel_func, l_indices)

    I = @index(Global, Cartesian) # I[1] == lm, I[2] == k
    # we use cartesian index instead of NTuple here
    # because this works for 2D and 3D matrices
    l = l_indices[I[1]]

    ∇²alms[I] = kernel_func(∇²alms[I], alms[I] * eigenvalues[l])
end

"""
$(TYPEDSIGNATURES)
Laplace operator ∇² applied to input `alms`, using precomputed eigenvalues from `S`.
Acts on the unit sphere, i.e. it omits 1/radius^2 scaling unless
`radius` keyword argument is provided."""
function ∇²(
        alms::LowerTriangularArray,             # Input: spectral coefficients
        S::AbstractSpectralTransform;           # precomputed eigenvalues
        kwargs...,
    )
    ∇²alms = similar(alms)
    ∇²!(∇²alms, alms, S; add = false, flipsign = false, inverse = false, kwargs...)
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
        ∇²alms::LowerTriangularArray,           # Input: spectral coefficients
        S::AbstractSpectralTransform;           # precomputed eigenvalues
        kwargs...,
    )
    alms = similar(∇²alms)
    ∇⁻²!(alms, ∇²alms, S; add = false, flipsign = false, kwargs...)
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
        ∇⁻²alms::LowerTriangularArray,          # Output: inverse Laplacian of alms
        alms::LowerTriangularArray,             # Input: spectral coefficients
        S::AbstractSpectralTransform;           # precomputed eigenvalues
        add::Bool = false,                      # add to output array or overwrite
        flipsign::Bool = false,                 # -∇⁻² or ∇⁻²
        kwargs...,
    )
    inverse = true
    return ∇²!(∇⁻²alms, alms, S; add, flipsign, inverse, kwargs...)
end

"""$(TYPEDSIGNATURES) Applies the gradient operator ∇ applied to input `p` and stores the result
in `dpdx` (zonal derivative) and `dpdy` (meridional derivative). The gradient operator acts
on the unit sphere and therefore omits the 1/radius scaling unless `radius` keyword argument is provided."""
function ∇!(
        dpdx::LowerTriangularArray,             # Output: zonal gradient
        dpdy::LowerTriangularArray,             # Output: meridional gradient
        p::LowerTriangularArray,                # Input: spectral coefficients
        S::AbstractSpectralTransform;           # includes precomputed arrays
        radius = DEFAULT_RADIUS,                # scale with radius if provided, otherwise unit sphere
    )
    (; grad_y1, grad_y2) = S.gradients
    (; m_indices) = p.spectrum
    @boundscheck ismatching(S, p) || throw(DimensionMismatch(S, p))

    # TODO: there's currently a scalar indexing error when using p direclty instead of p.data, this should be fixed
    @. dpdx = complex(0, m_indices - 1) * p.data

    launch!(architecture(dpdy), SpectralWorkOrder, size(dpdy), dpdy_kernel!, dpdy, p.data, grad_y1, grad_y2)

    # 1/radius factor if not unit sphere
    if radius != 1
        R⁻¹ = inv(radius)
        dpdx .*= R⁻¹
        dpdy .*= R⁻¹
    end

    return dpdx, dpdy
end

@kernel inbounds = true function dpdy_kernel!(dpdy, p, grad_y1, grad_y2)
    I = @index(Global, Cartesian)
    lm = I[1]
    k = ndims(p) == 1 ? CartesianIndex() : I[2]
    lmmax = size(dpdy, 1)

    gy1 = grad_y1[lm]
    gy2 = grad_y2[lm]

    # compared to the old CPU only version, some of the gy1 and gy2 are zero
    # that's why we don't need to check for l==m (diagonal) or l==p.m (last row)
    if lm == 1
        dpdy[lm, k] = gy2 * p[lm + 1, k]
    elseif lm == lmmax
        dpdy[lm, k] = gy1 * p[lm - 1, k]
    else
        dpdy[lm, k] = gy1 * p[lm - 1, k] + gy2 * p[lm + 1, k]
    end
end

"""$(TYPEDSIGNATURES)
Applies the gradient operator ∇ to the 2D input `p` and stores the zonal derivative
in layer 1 and the meridional derivative in layer 2 of the 3D output array `dpdxy`.
The gradient operator acts on the unit sphere and therefore omits the 1/radius scaling
unless `radius` keyword argument is provided.
This batched form avoids two separate spectral→grid transforms by computing both
gradient components into a single 3D array that can be transformed in one call."""
function ∇!(
        dpdxy::LowerTriangularArray,            # Output: 3D, layer 1 = zonal, layer 2 = meridional gradient
        p::LowerTriangularArray,                # Input: 2D spectral coefficients
        S::AbstractSpectralTransform;           # includes precomputed arrays
        radius = DEFAULT_RADIUS,                # scale with radius if provided, otherwise unit sphere
    )
    @boundscheck ismatching(S, p) || throw(DimensionMismatch(S, p))
    @boundscheck size(dpdxy, 2) >= 2 || throw(DimensionMismatch("dpdxy must have at least 2 layers"))

    (; grad_y1, grad_y2) = S.gradients
    (; m_indices) = p.spectrum

    launch!(architecture(dpdxy), SpectralWorkOrder, size(dpdxy), ∇_kernel!, dpdxy, p.data, m_indices, grad_y1, grad_y2)

    # 1/radius factor if not unit sphere
    if radius != 1
        dpdxy .*= inv(radius)
    end

    return dpdxy
end

@kernel inbounds = true function ∇_kernel!(dpdxy, p, m_indices, grad_y1, grad_y2)
    I = @index(Global, Cartesian)
    lm = I[1]
    k = I[2]
    lmmax = size(dpdxy, 1)

    # only compute for layer 1 (x) and layer 2 (y); ignore higher layers if present
    if k == 1
        # zonal gradient: dpdx = im*(m-1)*p
        dpdxy[lm, 1] = complex(0, m_indices[lm] - 1) * p[lm]
    elseif k == 2
        # meridional gradient using grad_y precomputed coefficients
        gy1 = grad_y1[lm]
        gy2 = grad_y2[lm]
        if lm == 1
            dpdxy[lm, 2] = gy2 * p[lm + 1]
        elseif lm == lmmax
            dpdxy[lm, 2] = gy1 * p[lm - 1]
        else
            dpdxy[lm, 2] = gy1 * p[lm - 1] + gy2 * p[lm + 1]
        end
    end
end

"""$(TYPEDSIGNATURES) The zonal and meridional gradient of `p`
using an existing `SpectralTransform` `S`. Acts on the unit sphere,
i.e. it omits 1/radius scaling unless `radius` keyword argument is provided."""
function ∇(p::LowerTriangularArray, S::AbstractSpectralTransform; kwargs...)
    dpdx = similar(p)
    dpdy = similar(p)
    ∇!(dpdx, dpdy, p, S; kwargs...)
    return dpdx, dpdy
end

"""$(TYPEDSIGNATURES) The zonal and meridional gradient of `p`.
Precomputes a `SpectralTransform` `S`. Acts on the unit-sphere,
i.e. it omits 1/radius scaling unless `radius` keyword argument is provided."""
function ∇(p::LowerTriangularArray; kwargs...)
    S = SpectralTransform(p, one_more_degree = true)
    return ∇(p, S; kwargs...)
end

"""$(TYPEDSIGNATURES) The zonal and meridional gradient of `grid`.
Transform to spectral space, takes the gradient and unscales the 1/coslat
scaling in the gradient. Acts on the unit-sphere, i.e. it omits 1/radius scaling unless
`radius` keyword argument is provided. Makes use of an existing spectral transform `S`."""
function ∇(field::AbstractField, S::AbstractSpectralTransform; kwargs...)
    p = transform(field, S)
    dpdx, dpdy = ∇(p, S; kwargs...)
    dpdx_grid = transform(dpdx, S, unscale_coslat = true)
    dpdy_grid = transform(dpdy, S, unscale_coslat = true)
    return dpdx_grid, dpdy_grid
end

"""$(TYPEDSIGNATURES) The zonal and meridional gradient of `grid`.
Transform to spectral space, takes the gradient and unscales the 1/coslat
scaling in the gradient. Acts on the unit-sphere, i.e. it omits 1/radius scaling unless
`radius` keyword argument is provided."""
function ∇(field::AbstractField; kwargs...)
    S = SpectralTransform(field, one_more_degree = true)
    return ∇(field, S; kwargs...)
end
