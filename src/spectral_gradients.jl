"""
    gradient_latitude!( coslat_u::AbstractArray{Complex{NF}},   # output: cos(lat)*zonal velocity u
                        Ψ::AbstractArray{Complex{NF}},          # input: streamfunction Ψ
                        ϵlms::AbstractArray{NF},                # recursion factors
                        R::Real=1                               # radius of the sphere/Earth
                        ) where {NF<:AbstractFloat}             # number format NF

Meridional gradient in spectral space of spherical harmonic coefficients `Ψ` on a sphere with
radius R. Returns `coslat_u`, i.e. the gradient ∂Ψ/∂lat with an additional cosine of latitude scaling.
This function uses the recursion relation (0-based degree l, order m)

    (coslat u)_lm = -1/R*(-(l-1)*ϵ_lm*Ψ_l-1,m + (l+2)*ϵ_l+1,m*Ψ_l+1,m ).
    
As u = -1/R*∂Ψ/∂lat, this function can be generally used to compute the gradient in latitude."""
function gradient_latitude!(coslat_u::AbstractMatrix{Complex{NF}},  # output: cos(lat)*zonal velocity u
                            Ψ::AbstractMatrix{Complex{NF}},         # input: streamfunction Ψ
                            S::SpectralTransform{NF};               # use precomputed recursion factors
                            flipsign::Bool=false,                   # flip the sign to obtain u from Ψ
                            add::Bool=false                         # coslat_u += (add) or = (overwrite)
                            ) where {NF<:AbstractFloat}             # number format NF

    # u needs one more degree/meridional mode l for each m than Ψ due to the recursion
    # Ψ can have size n+1 x n but then the last row is not used in the loop
    lmax_out, mmax_out = size(coslat_u)
    lmax_in,  mmax_in  = size(Ψ)

    @boundscheck mmax_out == mmax_in || throw(BoundsError)
    mmax = mmax_out - 1                 # 0-based max order m of harmonics

    @boundscheck abs(lmax_out-lmax_in) <= 1 || throw(BoundsError)
    lmax = lmax_in - 1                  # 0-based max degree l of harmonics
    output_larger = lmax_out - lmax_in

    if flipsign                     # used to get u from streamfunction Ψ (u ~ -∂Ψ/∂lat)
        grad_y1 = S.minus_grad_y1
        grad_y2 = S.minus_grad_y2
    else                            # no sign flip otherwise
        grad_y1 = S.grad_y1
        grad_y2 = S.grad_y2
    end

    g = grad_y2[1,1]*Ψ[2,1]
    coslat_u[1,1] = add ? coslat_u[1,1] + g : g         # l=m=0 mode only with term 2

    for m in 1:mmax+1
        for l in max(2,m):lmax
            g = grad_y1[l,m]*Ψ[l-1,m] + grad_y2[l,m]*Ψ[l+1,m]
            coslat_u[l,m] = add ? coslat_u[l,m] + g : g
        end
        for l in lmax+1:lmax+1+output_larger            # execute 2x for out > in, 1x for out == in and
            g = grad_y1[l,m]*Ψ[l-1,m]                   # 0x for out < in
            coslat_u[l,m] = add ? coslat_u[l,m] + g : g
        end
    end

    return coslat_u
end

function gradient_latitude( Ψ::AbstractMatrix{Complex{NF}}, # input: streamfunction Ψ
                            S::SpectralTransform{NF};       # precomputed gradient arrays
                            one_more_l::Bool=true,          # allocate output with one more degree l?
                            flipsign::Bool=false            # flip the sign to obtain u from Ψ?
                            ) where {NF<:AbstractFloat}     # number format NF
    _,mmax = size(Ψ) .- 1                                   # max degree l, order m of spherical harmonics
    coslat_u = zeros(Complex{NF},mmax+one_more_l+1,mmax+1)  # preallocate output, one more l for recursion
    return gradient_latitude!(coslat_u,Ψ,S;flipsign)        # call in-place version
end

"""
    coslat_v = gradient_longitude!( coslat_v::AbstractMatrix{Complex{NF}},
                                    Ψ::AbstractMatrix{Complex{NF}};
                                    radius::Real=1
                                    ) where {NF<:AbstractFloat}

Zonal gradient in spectral space of spherical harmonic coefficients `Ψ` on a sphere with radius `radius`.
While the zonal gradient has a 1/cos(lat) scaling in spherical coordinates, this functions omits the scaling
such that the returned array is scaled with coslat.
"""
function gradient_longitude!(   coslat_v::AbstractMatrix{Complex{NF}},  # output: cos(latitude)*meridional velocity
                                Ψ::AbstractMatrix{Complex{NF}},         # input: spectral coefficients of stream function
                                radius::Real=1;                         # radius of the sphere/Earth
                                add::Bool=false,                        # coslat_u += (add) or = (overwrite)
                                flipsign::Bool=false                    # flip sign of output?
                                ) where {NF<:AbstractFloat}             # number format NF

    lmax_out, mmax_out = size(coslat_v)
    lmax_in,  mmax_in  = size(Ψ)

    @boundscheck mmax_out == mmax_in || throw(BoundsError)
    mmax = mmax_out - 1                 # 0-based max order m of harmonics

    @boundscheck abs(lmax_out-lmax_in) <= 1 || throw(BoundsError)
    lmax = min(lmax_out,lmax_in) - 1    # 0-based max degree l of harmonics
    output_larger = lmax_out > lmax_in

    iradius⁻¹ = convert(Complex{NF},(-1)^flipsign*im/radius)            # = ±imaginary/radius converted to NF

    @inbounds for m in 1:mmax+1                     # loop over all coefficients, order m
        for l in m:lmax+1                           # degree l
            g = (m-1)*iradius⁻¹*Ψ[l,m]              # gradient in lon = *i*m/radius but 1-based order
            coslat_v[l,m] = add ? coslat_v[l,m] + g : g 
        end
    end

    # if grad in lon does not project onto the last degree l, set explicitly to zero in that case
    if output_larger & ~add
        for m in 1:mmax+1        
            coslat_v[end,m] = zero(Complex{NF}) 
        end
    end

    return coslat_v
end

"""Gradient in longitude in spectral space. Input: coefficients `alms` of the spherical harmonics."""
function gradient_longitude(Ψ::AbstractMatrix{NF},  # input array: spectral coefficients
                            R::Real=1;              # radius of the sphere/Earth
                            one_more_l::Bool=true   # allocate output with one more degree l
                            ) where NF              # number format NF
    _,mmax = size(Ψ)
    coslat_v = zeros(NF,mmax+one_more_l,mmax)       # preallocate output array (gradient in longitude)
    return gradient_longitude!(coslat_v,Ψ,R)        # call in-place version
end

"""
    curl!(  curl::AbstractMatrix{Complex{NF}},
            u::AbstractMatrix{Complex{NF}},
            v::AbstractMatrix{Complex{NF}},
            S::SpectralTransform{NF};
            flipsign::Bool=false,
            add::Bool=false
            ) where {NF<:AbstractFloat}

Curl of a vector `u,v` written into `curl`, `curl = ∇×(u,v)`. 
`u,v` are expected to have a 1/coslat-scaling included, then
`curl` is not scaled. `flipsign` option calculates -∇×(u,v) instead.
`add` option calculates `curl += ∇×(u,v)` instead. `flipsign` and `add`
can be combined. This functions only creates the kernel and calls
the generic divergence function _divergence! subsequently with flipped
u,v -> v,u for the curl."""
function curl!( curl::AbstractMatrix{Complex{NF}},
                u::AbstractMatrix{Complex{NF}},
                v::AbstractMatrix{Complex{NF}},
                S::SpectralTransform{NF};
                flipsign::Bool=false,
                add::Bool=false
                ) where {NF<:AbstractFloat}

    # = -(∂λ - ∂θ) or (∂λ - ∂θ), adding or overwriting the output curl
    kernel(o,a,b,c) = flipsign ? (add ? o-(a+b-c) : -(a+b-c)) :
                                 (add ? o+(a+b-c) :   a+b-c )    
    _divergence!(kernel,curl,v,u,S)             # flip u,v -> v,u
end

"""
    divergence!(div::AbstractMatrix{Complex{NF}},
                u::AbstractMatrix{Complex{NF}},
                v::AbstractMatrix{Complex{NF}},
                S::SpectralTransform{NF};
                flipsign::Bool=false,
                add::Bool=false
                ) where {NF<:AbstractFloat}

Divergence of a vector `u,v` written into `div`, `div = ∇⋅(u,v)`. 
`u,v` are expected to have a 1/coslat-scaling included, then
`div` is not scaled. `flipsign` option calculates -∇⋅(u,v) instead.
`add` option calculates `div += ∇⋅(u,v)` instead. `flipsign` and `add`
can be combined. This functions only creates the kernel and calls
the generic divergence function _divergence! subsequently."""
function divergence!(   div::AbstractMatrix{Complex{NF}},
                        u::AbstractMatrix{Complex{NF}},
                        v::AbstractMatrix{Complex{NF}},
                        S::SpectralTransform{NF};
                        flipsign::Bool=false,
                        add::Bool=false
                        ) where {NF<:AbstractFloat}

    # = -(∂λ + ∂θ) or (∂λ + ∂θ), adding or overwriting the output div
    kernel(o,a,b,c) = flipsign ? (add ? o-(a-b+c) : -(a+b-c)) :
                                 (add ? o+(a-b+c) :   a-b+c )                
    _divergence!(kernel,div,u,v,S)
end

"""
    _divergence!(   kernel,
                    div::AbstractMatrix{Complex{NF}},
                    u::AbstractMatrix{Complex{NF}},
                    v::AbstractMatrix{Complex{NF}},
                    S::SpectralTransform{NF}
                    ) where {NF<:AbstractFloat}

Generic divergence function of vector `u`,`v` that writes into the output into `div`.
Generic as it uses the kernel `kernel` such that curl, div, add or flipsign
options are provided through `kernel`, but otherwise a single function is used."""
function _divergence!(  kernel,
                        div::AbstractMatrix{Complex{NF}},
                        u::AbstractMatrix{Complex{NF}},
                        v::AbstractMatrix{Complex{NF}},
                        S::SpectralTransform{NF}
                        ) where {NF<:AbstractFloat}

    lmax,mmax = size(div) .- (1,1)                  # 0-based lmax,mmax 
    @boundscheck size(u) == (lmax+2,mmax+1) || throw(BoundsError)
    @boundscheck size(v) == (lmax+2,mmax+1) || throw(BoundsError)

    @unpack grad_y_vordiv1,grad_y_vordiv2 = S

    z = zero(Complex{NF})
    div[1,1] = kernel(div[1,1],z,z,z)               # l=m=0 harmonic is zero

    @inbounds for m in 1:mmax+1                     # 1-based l,m
        for l in max(2,m):lmax+1                    # skip l=m=0 harmonic (mean) to avoid access to v[0,1]
            ∂u∂λ  = ((m-1)*im)*u[l,m]
            ∂v∂θ1 = grad_y_vordiv1[l,m]*v[l-1,m]
            ∂v∂θ2 = grad_y_vordiv2[l,m]*v[l+1,m]
            div[l,m] = kernel(div[l,m],∂u∂λ,∂v∂θ1,∂v∂θ2)
        end
    end

    return nothing
end

"""Compute curl and divergence of vectors u,v in one go. Currently not used."""
function curl_div!( curl::AbstractMatrix{Complex{NF}},
                    div::AbstractMatrix{Complex{NF}},
                    u::AbstractMatrix{Complex{NF}},
                    v::AbstractMatrix{Complex{NF}},
                    S::SpectralTransform{NF}
                    ) where {NF<:AbstractFloat}

    lmax,mmax = size(div) .- (1,1)                  # 0-based lmax,mmax
    @boundscheck size(curl) == size(div) || throw(BoundsError)
    @boundscheck size(u) == (lmax+2,mmax+1) || throw(BoundsError)
    @boundscheck size(v) == (lmax+2,mmax+1) || throw(BoundsError)

    @unpack grad_y_vordiv1,grad_y_vordiv2 = S

    div[1,1]  = zero(Complex{NF})                   # l=m=0 harmonic is zero
    curl[1,1] = zero(Complex{NF})                   # l=m=0 harmonic is zero
    @inbounds for m in 1:mmax+1                     # 1-based l,m
        for l in max(2,m):lmax+1                    # skip l=m=0 harmonic (mean) to avoid access to v[0,1]
            ∂u∂λ  = ((m-1)*im)*u[l,m]
            ∂v∂λ  = ((m-1)*im)*v[l,m]
            ∂v∂θ1 = grad_y_vordiv1[l,m]*v[l-1,m]
            ∂u∂θ1 = grad_y_vordiv1[l,m]*u[l-1,m]
            ∂v∂θ2 = grad_y_vordiv2[l,m]*v[l+1,m]
            ∂u∂θ2 = grad_y_vordiv2[l,m]*u[l+1,m]
            div[l,m]  = ∂u∂λ - ∂v∂θ1 + ∂v∂θ2
            curl[l,m] = ∂v∂λ + ∂u∂θ1 - ∂u∂θ2
        end
    end

    return nothing
end

"""
    UV_from_vor!(   U::AbstractMatrix{Complex{NF}},
                    V::AbstractMatrix{Complex{NF}},
                    vor::AbstractMatrix{Complex{NF}},
                    S::SpectralTransform{NF}
                    ) where {NF<:AbstractFloat}

Get U,V (=(u,v)*coslat) from vorticity ζ spectral space (divergence D=0)
Two operations are combined into a single linear operation. First, invert the
spherical Laplace ∇² operator to get stream function from vorticity. Then
compute zonal and meridional gradients to get U,V."""
function UV_from_vor!(  U::AbstractMatrix{Complex{NF}},
                        V::AbstractMatrix{Complex{NF}},
                        vor::AbstractMatrix{Complex{NF}},
                        S::SpectralTransform{NF}
                        ) where {NF<:AbstractFloat}

    lmax,mmax = size(vor) .- (1,1)                  # 0-based lmax,mmax
    @boundscheck size(U) == (lmax+2,mmax+1) || throw(BoundsError)
    @boundscheck size(V) == (lmax+2,mmax+1) || throw(BoundsError)

    @unpack vordiv_to_uv_x,vordiv_to_uv1,vordiv_to_uv2 = S

    U[1,1] = vordiv_to_uv2[1,1]*vor[2,1]            # l=m=0 harmonic has only one contribution
    V[1,1] = zero(Complex{NF})                      # obtained via zonal derivative (*i*m) = 0 

    @inbounds for m in 1:mmax+1                     # 1-based l,m
        for l in max(2,m):lmax                      # skip l=m=0 harmonic (mean) to avoid access to v[0,1]
            # meridional gradient, u = -∂/∂lat(Ψ), omit radius R scaling
            U[l,m] = vordiv_to_uv2[l,m]*vor[l+1,m] - vordiv_to_uv1[l,m]*vor[l-1,m]

            # zonal gradient, V = ∂/∂λ(Ψ), omit radius R scaling
            V[l,m] = im*vordiv_to_uv_x[l,m]*vor[l,m]
        end
    end

    @inbounds for m in 1:mmax+1
        l = lmax+1                                  # second last row, l+1 index is out of bounds (=0 entry)
        U[l,m] = -vordiv_to_uv1[l,m]*vor[l-1,m]     # meridional gradient again (but only 2nd term from above)
        V[l,m] = im*vordiv_to_uv_x[l,m]*vor[l,m]    # zonal gradient again (as above)

        l = lmax+2                                  # last row, l,l+1 indices are out of bounds (=0 entries)
        U[l,m] = -vordiv_to_uv1[l,m]*vor[l-1,m]     # meridional gradient again (but only 2nd term from above)
        V[l,m] = zero(Complex{NF})                  # set explicitly to 0 as Ψ does not contribute to last row of V
    end
end

"""
    UV_from_vordiv!(U::AbstractMatrix{Complex{NF}},
                    V::AbstractMatrix{Complex{NF}},
                    vor::AbstractMatrix{Complex{NF}},
                    div::AbstractMatrix{Complex{NF}},
                    S::SpectralTransform{NF}
                    ) where {NF<:AbstractFloat}

Get U,V (=(u,v)*coslat) from vorticity ζ and divergence D in spectral space.
Two operations are combined into a single linear operation. First, invert the
spherical Laplace ∇² operator to get stream function from vorticity and
velocity potential from divergence. Then compute zonal and meridional gradients
to get U,V."""
function UV_from_vordiv!(   U::AbstractMatrix{Complex{NF}},
                            V::AbstractMatrix{Complex{NF}},
                            vor::AbstractMatrix{Complex{NF}},
                            div::AbstractMatrix{Complex{NF}},
                            S::SpectralTransform{NF}
                            ) where {NF<:AbstractFloat}

    lmax,mmax = size(vor) .- (1,1)                  # 0-based lmax,mmax
    @boundscheck size(vor) == size(div) || throw(BoundsError)
    @boundscheck size(U) == (lmax+2,mmax+1) || throw(BoundsError)
    @boundscheck size(V) == (lmax+2,mmax+1) || throw(BoundsError)

    @unpack vordiv_to_uv_x,vordiv_to_uv1,vordiv_to_uv2 = S

    U[1,1] =  vordiv_to_uv2[1,1]*vor[2,1]           # l=m=0 harmonic has only one contribution
    V[1,1] = -vordiv_to_uv2[1,1]*div[2,1]

    @inbounds for m in 1:mmax+1                     # 1-based l,m
        for l in max(2,m):lmax                      # skip l=m=0 harmonic (mean) to avoid access to v[0,1]
            ∂Dλ = im*vordiv_to_uv_x[l,m]*div[l,m]   # divergence contribution to zonal gradient
            ∂ζλ = im*vordiv_to_uv_x[l,m]*vor[l,m]   # vorticity contribution to zonal gradient

            # div,vor contribution to meridional gradient
            ∂ζθ = vordiv_to_uv2[l,m]*vor[l+1,m] - vordiv_to_uv1[l,m]*vor[l-1,m]
            ∂Dθ = vordiv_to_uv1[l,m]*div[l-1,m] - vordiv_to_uv2[l,m]*div[l+1,m]
            U[l,m] = ∂Dλ + ∂ζθ
            V[l,m] = ∂ζλ + ∂Dθ
        end
    end

    @inbounds for m in 1:mmax+1
        l = lmax+1                                  # second last row, l+1 index is out of bounds (=0 entry)
        U[l,m] = im*vordiv_to_uv_x[l,m]*div[l,m] - vordiv_to_uv1[l,m]*vor[l-1,m]
        V[l,m] = im*vordiv_to_uv_x[l,m]*vor[l,m] + vordiv_to_uv1[l,m]*div[l-1,m]

        l = lmax+2                                  # last row, l,l+1 indices are out of bounds (=0 entries)
        U[l,m] = -vordiv_to_uv1[l,m]*vor[l-1,m]
        V[l,m] =  vordiv_to_uv1[l,m]*div[l-1,m]
    end
end

# alias functions to scale the latitude of lat-lon map `A`
scale_coslat!(  A::AbstractMatrix,G::Geometry) = _scale_lat!(A,G.coslat)
scale_coslat²!( A::AbstractMatrix,G::Geometry) = _scale_lat!(A,G.coslat²)
scale_coslat⁻¹!(A::AbstractMatrix,G::Geometry) = _scale_lat!(A,G.coslat⁻¹)
scale_coslat⁻²!(A::AbstractMatrix,G::Geometry) = _scale_lat!(A,G.coslat⁻²)

"""
    _scale_lat!(A::AbstractMatrix{NF},v::AbstractVector) where {NF<:AbstractFloat}

Generic latitude scaling applied to `A` in-place with latitude-like vector `v`."""
function _scale_lat!(A::AbstractMatrix{NF},v::AbstractVector) where {NF<:AbstractFloat}
    nlon,nlat = size(A)
    @boundscheck nlat == length(v) || throw(BoundsError)

    @inbounds for j in 1:nlat
        vj = convert(NF,v[j])
        for i in 1:nlon
            A[i,j] *= vj
        end
    end
end 

# In-place update of spectral coeffs alms with their (inverse) Laplace operator-ed version 
∇⁻²!(alms::AbstractMatrix{Complex{NF}},S::SpectralTransform{NF}) where NF = alms .* S.eigenvalues⁻¹
∇²!( alms::AbstractMatrix{Complex{NF}},S::SpectralTransform{NF}) where NF = alms .* S.eigenvalues

"""
    ∇⁻²!(   ∇⁻²alms::AbstractMatrix{Complex},
            alms::AbstractMatrix{Complex},
            S::SpectralTransform)

Inverse Laplace operator ∇⁻² applied to the spectral coefficients `alms` in spherical
coordinates. The radius `R` is omitted in the eigenvalues which are precomputed in `S`.
∇⁻²! is the in-place version which directly stores the output in the first argument `∇⁻²alms`.
The integration constant for Legendre polynomial `l=m=0` is zero. The inverse spherical
Laplace operator is generally

    ∇⁻²alms = alms*R²/(-l(l+1))

with the degree `l` (0-based) of the Legendre polynomial."""
function ∇⁻²!(  ∇⁻²alms::AbstractMatrix{Complex{NF}},   # Output: inverse Laplacian of alms
                alms::AbstractMatrix{Complex{NF}},      # spectral coefficients
                S::SpectralTransform{NF}                # precomputed arrays for spectral space
                ) where {NF<:AbstractFloat}

    @boundscheck size(alms) == size(∇⁻²alms) || throw(BoundsError)
    lmax,mmax = size(alms) .- 1     # degree l, order m of the Legendre polynomials
    
    @unpack eigenvalues⁻¹ = S
    @boundscheck length(eigenvalues⁻¹) >= lmax+1 || throw(BoundsError)

    @inbounds for m in 1:mmax+1     # order m = 0:mmax but 1-based
        for l in m:lmax+1           # degree l = m:lmax but 1-based
            ∇⁻²alms[l,m] = alms[l,m]*eigenvalues⁻¹[l]
        end
    end

    return ∇⁻²alms
end

"""
    ∇²!(    ∇⁻²alms::AbstractMatrix{Complex},
            alms::AbstractMatrix{Complex},
            S::SpectralTransform)

Laplace operator ∇² applied to the spectral coefficients `alms` in spherical
coordinates. The radius `R` is omitted in the eigenvalues which are precomputed in `S`.
∇²! is the in-place version which directly stores the output in the first argument `∇²alms`.
The spherical Laplace operator is generally

    ∇⁻²alms = alms*(-l(l+1))/R²

with the degree `l` (0-based) of the Legendre polynomial."""
function ∇²!(   ∇²alms::AbstractMatrix{Complex{NF}},    # Output: Laplacian of alms
                alms::AbstractMatrix{Complex{NF}},      # spectral coefficients
                S::SpectralTransform{NF}                # precomputed arrays for spectral space
                ) where {NF<:AbstractFloat}

    @boundscheck size(alms) == size(∇²alms) || throw(BoundsError)
    lmax,mmax = size(alms) .- 1     # degree l, order m of the Legendre polynomials
    
    @unpack eigenvalues = S
    @boundscheck length(eigenvalues) >= lmax+1 || throw(BoundsError)

    @inbounds for m in 1:mmax+1     # order m = 0:mmax but 1-based
        for l in m:lmax+1           # degree l = m:lmax but 1-based
            ∇²alms[l,m] = alms[l,m]*eigenvalues[l]
        end
    end

    return ∇²alms
end