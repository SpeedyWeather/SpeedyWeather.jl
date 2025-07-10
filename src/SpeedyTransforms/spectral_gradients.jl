const DEFAULT_RADIUS = 1

"""
$(TYPEDSIGNATURES)
Recursion factors `ϵ` as a function of degree `l` and order `m` (0-based) of the spherical harmonics.
ϵ(l, m) = sqrt((l^2-m^2)/(4*l^2-1))."""
recursion_factor(l, m) = sqrt((l^2-m^2)/(4*l^2-1))

"""
$(TYPEDSIGNATURES)     
Returns a matrix of recursion factors `ϵ` up to degree `lmax` and order `mmax` (1-based) of the `spectrum` in number format `NF`."""
function recursion_factors(
    ::Type{NF},             # number format NF
    spectrum::Spectrum,
) where NF

    ϵ_lm = zeros(NF, spectrum)                          # store in lower triangular matrix
    for (m, degrees) in enumerate(orders(spectrum))     # loop over 1-based l, m
        for l in degrees
            ϵ_lm[l, m] = recursion_factor(l-1, m-1)     # convert to 0-based l, m for function arguments
        end
    end
    return ϵ_lm
end

# if number format not provided use DEFAULT_NF
recursion_factors(spectrum::Spectrum) = recursion_factors(DEFAULT_NF, spectrum)

"""$(TYPEDSIGNATURES)
Precomputed factors used for meridional gradients.
On unit sphere, hence 1/radius-scaling is omitted."""
function meridional_gradient_factors(
    ::Type{NF},             # number format NF
    spectrum::Spectrum
) where NF

    # get recursion factors for spherical harmonics
    (;lmax, mmax) = spectrum
    ϵ_lm = recursion_factors(NF, Spectrum(lmax+1, mmax))

    # meridional gradient for scalars (coslat scaling included)
    grad_y1 = zeros(NF, spectrum)                           # term 1, mul with harmonic l-1, m
    grad_y2 = zeros(NF, spectrum)                           # term 2, mul with harmonic l+1, m

    for (m, degrees) in enumerate(orders(spectrum))         # 1-based degree l, order m
        for l in degrees        
            grad_y1[l, m] = -(l-2)*ϵ_lm[l, m]
            grad_y2[l, m] = (l+1)*ϵ_lm[l+1, m]
        end
        # explicitly set the last row to zero, so that kernels yield correct gradient in last row 
        grad_y2[degrees[end], m] = 0
    end

    # meridional gradient used to get from u, v/coslat to vorticity and divergence
    grad_y_vordiv1 = zeros(NF, spectrum)                    # term 1, mul with harmonic l-1, m
    grad_y_vordiv2 = zeros(NF, spectrum)                    # term 2, mul with harmonic l+1, m
 
    for (m, degrees) in enumerate(orders(spectrum))         # 1-based degree l, order m
        for l in degrees           
            grad_y_vordiv1[l, m] = l*ϵ_lm[l, m]
            grad_y_vordiv2[l, m] = (l-1)*ϵ_lm[l+1, m]
        end
        # explicitly set the last row to zero, so that kernels yield correct gradient in last row 
        grad_y_vordiv1[degrees[end], m] = 0
        grad_y_vordiv2[degrees[end], m] = 0
    end

    return grad_y1, grad_y2, grad_y_vordiv1, grad_y_vordiv2
end

# if number format not provided use default
meridional_gradient_factors(spectrum::Spectrum) = meridional_gradient_factors(DEFAULT_NF, spectrum)

# zonal gradient used to get from u, v/coslat to vorticity and divergence
function zonal_gradient_factors(
    spectrum::Spectrum,
)
    grad_x_vordiv = zeros(Int, spectrum)                    # allocate as int, can be converted later
    for (m, degrees) in enumerate(orders(spectrum))         # 1-based degree l, order m
        for l in degrees                   
            grad_x_vordiv[l, m] = m-1                       # just the 0-based order (=zonal wavenumber), omit imaginary unit
        end # last row zero to get vor and div correct
        grad_x_vordiv[degrees[end], m] = 0
    end

    return grad_x_vordiv
end

# zonal integration (sort of) to get from vorticity and divergence to u, v*coslat
function zonal_integration_factors(
    ::Type{NF},             # number format NF
    spectrum::Spectrum,
) where NF
    vordiv_to_uv_x = zeros(NF, spectrum)
    for (m, degrees) in enumerate(orders(spectrum))         # 1-based degree l, order m
        for l in degrees
            # -m/(l*(l+1)) is the 0-based integration factor, 1-based with m->m-1, l->l-1
            vordiv_to_uv_x[l, m] = -(m-1)/(l*(l-1))              
        end
    end
    vordiv_to_uv_x[1] = 0                                   # remove NaN from 0/0
    return vordiv_to_uv_x
end

# if number format not provided use default
zonal_integration_factors(spectrum::Spectrum) = zonal_integration_factors(DEFAULT_NF, spectrum)

# meridional integration (sort of) to get from vorticity and divergence to u, v*coslat
function meridional_integration_factors(
    ::Type{NF},             # number format NF
    spectrum::Spectrum,
) where NF

    # get recursion factors for spherical harmonics
    (;lmax, mmax) = spectrum
    ϵ_lm = recursion_factors(NF, Spectrum(lmax+1, mmax))

    vordiv_to_uv1 = zeros(NF, spectrum)                     # term 1, to be mul with harmonic l-1, m
    vordiv_to_uv2 = zeros(NF, spectrum)                     # term 2, to be mul with harmonic l+1, m

    for (m, degrees) in enumerate(orders(spectrum))         # 1-based degree l, order m
        for l in degrees             
            vordiv_to_uv1[l, m] = ϵ_lm[l, m]/(l-1)
            vordiv_to_uv2[l, m] = ϵ_lm[l+1, m]/l
        end
        # explicitly set the last row of vordiv_to_uv2 to zero, so that kernels yield correct gradient in last row 
        vordiv_to_uv2[degrees[end], m] = 0
    end

    vordiv_to_uv1[1] = 0                 # remove NaN from 0/0
    return vordiv_to_uv1, vordiv_to_uv2
end

# if number format not provided use default
meridional_integration_factors(spectrum::Spectrum) = meridional_integration_factors(DEFAULT_NF, spectrum)

"""$(TYPEDSIGNATURES)
Get the eigenvalues of the spherical harmonics up to degree `lmax` and order `mmax`
as determined by the `spectrum` in number format `NF`."""
function get_eigenvalues(::Type{NF}, spectrum::Spectrum) where NF
    (; lmax, mmax) = spectrum
    eigenvalues = [NF(-l*(l+1)) for l in 0:lmax-1]      # 0:lmax-1 for 0-based indexing
end

get_eigenvalues(spectrum::Spectrum) = get_eigenvalues(DEFAULT_NF, spectrum)

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
@inline (::KernelOP{true, false, false})(o, a, b, c) = a+b-c

# Curl with flipsign (no add)
@inline (::KernelOP{true, true, false})(o, a, b, c) = -(a+b-c)

# Curl with add (no flipsign)
@inline (::KernelOP{true, false, true})(o, a, b, c) = o+(a+b-c)

# Curl with flipsign and add
@inline (::KernelOP{true, true, true})(o, a, b, c) = o-(a+b-c)

# Divergence operations (mode == false): a-b+c
# Standard divergence (no flipsign, no add)
@inline (::KernelOP{false, false, false})(o, a, b, c) = a-b+c

# Divergence with flipsign (no add)
@inline (::KernelOP{false, true, false})(o, a, b, c) = -(a-b+c)

# Divergence with add (no flipsign)
@inline (::KernelOP{false, false, true})(o, a, b, c) = o+(a-b+c)

# Divergence with flipsign and add
@inline (::KernelOP{false, true, true})(o, a, b, c) = o-(a-b+c)


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
    kernel = KernelOP{true, flipsign, add}()
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
    kernel = KernelOP{false, flipsign, add}()
    _divergence!(kernel, div, u, v, S; kwargs...)
end

function _divergence!(  
    kernel,
    div::LowerTriangularArray,
    u::LowerTriangularArray,
    v::LowerTriangularArray,
    S::SpectralTransform;
    radius = DEFAULT_RADIUS,
)
    (; grad_x_vordiv, grad_y_vordiv1, grad_y_vordiv2) = S
  
    @boundscheck ismatching(S, div) || throw(DimensionMismatch(S, div))

    launch!(architecture(div), :lmk, size(div), _divergence_kernel!, kernel, div, u, v, grad_x_vordiv, grad_y_vordiv1, grad_y_vordiv2)
    synchronize(architecture(div))
    
    # radius scaling if not unit sphere
    if radius != 1
        div .*= inv(radius)
    end

    return div
end

@kernel inbounds=true function _divergence_kernel!(kernel_func::KernelOP{mode, flipsign, add}, div, u, v, grad_x_vordiv, grad_y_vordiv1, grad_y_vordiv2) where {mode, flipsign, add}

    I = @index(Global, Cartesian)
    lm = I[1]
    lmmax = size(div, 1)
    k = ndims(div) == 1 ? CartesianIndex() : I[2]

    if lm == 1
        div[I] = kernel_func(div[I], 0, 0, grad_y_vordiv2[1]*v[2, k])
    elseif lm == lmmax    
        div[I] = 0
    else
        ∂u∂λ  = im*grad_x_vordiv[lm]*u[I]
        ∂v∂θ1 = grad_y_vordiv1[lm] * v[lm-1, k] 
        ∂v∂θ2 = grad_y_vordiv2[lm] * v[lm+1, k]  
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
    u::AbstractField,
    v::AbstractField;
    kwargs...,
)
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
    
    launch!(architecture(U), :lmk, size(U), _UV_from_vor_kernel!, U, V, vor, vor.spectrum.l_indices, vor.spectrum.lmax, vordiv_to_uv_x, vordiv_to_uv1, vordiv_to_uv2)
    synchronize(architecture(U))
    
    # *radius scaling if not unit sphere (*radius² for ∇⁻² then /radius to get from stream function to velocity)
    if radius != 1
        U .*= radius
        V .*= radius
    end

    return U, V
end

@kernel inbounds=true function _UV_from_vor_kernel!(U, V, vor, @Const(l_indices), lmax, vordiv_to_uv_x, vordiv_to_uv1, vordiv_to_uv2)    
    I = @index(Global, Cartesian)
    lm = I[1]
    k = ndims(vor) == 1 ? CartesianIndex() : I[2]
    l = l_indices[lm]

    # Get the coefficients for the current lm index
    z = im*vordiv_to_uv_x[lm]
    vordiv_uv1 = vordiv_to_uv1[lm]
    vordiv_uv2 = vordiv_to_uv2[lm]
    
    # Handle different cases based on position in the triangular matrix
    if lm == 1  # First element (diagonal)
        # U = -∂/∂lat(Ψ) - no lm-1 term for first element
        U[I] = vordiv_uv2 * vor[lm+1, k]
        # V = ∂/∂λ(Ψ)
        V[I] = z * vor[I]
    elseif l==(lmax-1)     # extra in case vor[lmax,:] != 0, see comment in UV_from_vordiv!
        U[I] = -vordiv_uv1 * vor[lm-1, k]    # meridional gradient again (but only 2nd term from above)
        V[I] = z * vor[I]          # zonal gradient again (as above)
    elseif l==lmax         # extra in case vor[lmax,:] != 0, see comment in UV_from_vordiv!
        U[I] = -vordiv_uv1 * vor[lm-1, k]
        V[I] = 0
    else  # General case (below diagonal)
        # U = -∂/∂lat(Ψ) combined with Laplace inversion ∇⁻²
        U[I] = muladd(vordiv_uv2, vor[lm+1, k], -vordiv_uv1*vor[lm-1, k])
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

# New kernel-based implementation, currently unused because of problems with Enzyme
function UV_from_vordiv_kernel!(   
    U::LowerTriangularArray,
    V::LowerTriangularArray,
    vor::LowerTriangularArray,
    div::LowerTriangularArray,
    S::SpectralTransform;
    radius = DEFAULT_RADIUS,
)
    (; vordiv_to_uv_x, vordiv_to_uv1, vordiv_to_uv2 ) = S
    @boundscheck ismatching(S, U) || throw(DimensionMismatch(S, U))

    launch!(architecture(U), :lmk, size(U), _UV_from_vordiv_kernel!, U, V, vor, div, vor.spectrum.l_indices, vor.spectrum.lmax, vordiv_to_uv_x, vordiv_to_uv1, vordiv_to_uv2)
    synchronize(architecture(U))
    
    # *radius scaling if not unit sphere (*radius² for ∇⁻², then /radius to get from stream function to velocity)
    if radius != 1
        U .*= radius
        V .*= radius
    end

    return U, V
end

@kernel inbounds=true function _UV_from_vordiv_kernel!(U, V, vor, div, @Const(l_indices), lmax, vordiv_to_uv_x, vordiv_to_uv1, vordiv_to_uv2)    
    I = @index(Global, Cartesian)
    lm = I[1]
    k = ndims(vor) == 1 ? CartesianIndex() : I[2]
    l = l_indices[lm]
    
    # Get the coefficients for the current lm index
    z = vordiv_to_uv_x[lm]
    vordiv_uv1 = vordiv_to_uv1[lm]
    vordiv_uv2 = vordiv_to_uv2[lm]
    
    # Handle different cases based on position in the triangular matrix
    if lm == 1  # First element (diagonal)
        # Meridional gradient contributions
        ∂ζθ =  vordiv_uv2 * vor[lm+1, k]       # lm-1 term is zero
        ∂Dθ = -vordiv_uv2 * div[lm+1, k]       # lm-1 term is zero
        
        # Zonal gradient contributions
        U[I] = muladd(z, div[I], ∂ζθ)          # = ∂Dλ + ∂ζθ
        V[I] = muladd(z, vor[I], ∂Dθ)          # = ∂ζλ + ∂Dθ
    elseif l == (lmax-1)  # Second last row
        U[I] = muladd(z, div[I], -vordiv_uv1 * vor[lm-1, k])
        V[I] = muladd(z, vor[I], vordiv_uv1 * div[lm-1, k])
    elseif l == lmax  # Last row
        U[I] = -vordiv_uv1 * vor[lm-1, k]      # only last term from 2nd last row
        V[I] =  vordiv_uv1 * div[lm-1, k]      # only last term from 2nd last row
    else  # General case (below diagonal)
        # Meridional gradient contributions
        ∂ζθ = muladd(vordiv_uv2, vor[lm+1, k], -vordiv_uv1 * vor[lm-1, k])
        ∂Dθ = muladd(vordiv_uv1, div[lm-1, k], -vordiv_uv2 * div[lm+1, k])
        
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
    
    launch!(architecture(∇²alms), :lmk, size(∇²alms), ∇²_kernel!, ∇²alms, alms, eigenvalues, kernel, alms.spectrum.l_indices)
    synchronize(architecture(∇²alms))

    # /radius² or *radius² scaling if not unit sphere
    if radius != 1
        R_plusminus_squared = inverse ? radius^2 : inv(radius^2)
        ∇²alms .*= R_plusminus_squared
    end

    return ∇²alms
end

@kernel function ∇²_kernel!(∇²alms, alms, @Const(eigenvalues), kernel_func, @Const(l_indices))

    I = @index(Global, Cartesian) # I[1] == lm, I[2] == k
                                  # we use cartesian index instead of NTuple here
                                  # because this works for 2D and 3D matrices
    l = l_indices[I[1]]

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
    (; grad_y1, grad_y2) = S
    (; m_indices) = p.spectrum
    @boundscheck ismatching(S, p) || throw(DimensionMismatch(S, p))

    # TODO: there's currently a scalar indexing error when using p direclty instead of p.data, this should be fixed
    @. dpdx = complex(0, m_indices - 1)*p.data

    launch!(architecture(dpdy), :lmk, size(dpdy), dpdy_kernel!, dpdy, p.data, grad_y1, grad_y2)
    synchronize(architecture(dpdy))
    
    # 1/radius factor if not unit sphere
    if radius != 1
        R⁻¹ = inv(radius)
        dpdx .*= R⁻¹
        dpdy .*= R⁻¹
    end

    return dpdx, dpdy
end

@kernel inbounds=true function dpdy_kernel!(dpdy, p, grad_y1, grad_y2)
    I = @index(Global, Cartesian)
    lm = I[1]
    k = ndims(p) == 1 ? CartesianIndex() : I[2]
    lmmax = size(dpdy, 1)

    gy1 = grad_y1[lm]
    gy2 = grad_y2[lm]

    # compared to the old CPU only version, some of the gy1 and gy2 are zero 
    # that's why we don't need to check for l==m (diagonal) or l==p.m (last row)
    if lm == 1
        dpdy[lm, k] = gy2*p[lm+1, k]
    elseif lm == lmmax
        dpdy[lm, k] = gy1*p[lm-1, k]
    else
        dpdy[lm, k] = gy1*p[lm-1, k] + gy2*p[lm+1, k]
    end
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
function ∇(field::AbstractField, S::SpectralTransform; kwargs...)
    p = transform(field, S)
    dpdx, dpdy = ∇(p, S; kwargs...)
    dpdx_grid = transform(dpdx, S, unscale_coslat=true)
    dpdy_grid = transform(dpdy, S, unscale_coslat=true)
    return dpdx_grid, dpdy_grid
end

"""$(TYPEDSIGNATURES) The zonal and meridional gradient of `grid`.
Transform to spectral space, takes the gradient and unscales the 1/coslat
scaling in the gradient. Acts on the unit-sphere, i.e. it omits 1/radius scaling unless
`radius` keyword argument is provided."""
function ∇(field::AbstractField; kwargs...)
    S = SpectralTransform(field, one_more_degree=true)
    return ∇(field, S; kwargs...)
end
