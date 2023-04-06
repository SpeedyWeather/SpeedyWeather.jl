"""
    curl!(  curl::LowerTriangularMatrix,
            u::LowerTriangularMatrix,
            v::LowerTriangularMatrix,
            S::SpectralTransform;
            flipsign::Bool=false,
            add::Bool=false,
            )

Curl of a vector `u,v` written into `curl`, `curl = ∇×(u,v)`.
`u,v` are expected to have a 1/coslat-scaling included, then `curl` is not scaled.
`flipsign` option calculates -∇×(u,v) instead. `add` option calculates `curl += ∇×(u,v)` instead.
`flipsign` and `add` can be combined. This functions only creates the kernel and calls the generic
divergence function _divergence! subsequently with flipped u,v -> v,u for the curl."""
function curl!( curl::LowerTriangularMatrix,
                u::LowerTriangularMatrix,
                v::LowerTriangularMatrix,
                S::SpectralTransform;
                flipsign::Bool=false,
                add::Bool=false,
                )

    # = -(∂λ - ∂θ) or (∂λ - ∂θ), adding or overwriting the output curl
    kernel(o,a,b,c) = flipsign ? (add ? o-(a+b-c) : -(a+b-c)) :
                                 (add ? o+(a+b-c) :   a+b-c )    
    _divergence!(kernel,curl,v,u,S)             # flip u,v -> v,u
end

"""
    divergence!(div::LowerTriangularMatrix,
                u::LowerTriangularMatrix,
                v::LowerTriangularMatrix,
                S::SpectralTransform{NF};
                flipsign::Bool=false,
                add::Bool=false,
                )

Divergence of a vector `u,v` written into `div`, `div = ∇⋅(u,v)`. 
`u,v` are expected to have a 1/coslat-scaling included, then `div` is not scaled.
`flipsign` option calculates -∇⋅(u,v) instead. `add` option calculates `div += ∇⋅(u,v)` instead.
`flipsign` and `add` can be combined. This functions only creates the kernel and calls
the generic divergence function _divergence! subsequently."""
function divergence!(   div::LowerTriangularMatrix,
                        u::LowerTriangularMatrix,
                        v::LowerTriangularMatrix,
                        S::SpectralTransform;
                        flipsign::Bool=false,
                        add::Bool=false,
                        )

    # = -(∂λ + ∂θ) or (∂λ + ∂θ), adding or overwriting the output div
    kernel(o,a,b,c) = flipsign ? (add ? o-(a-b+c) : -(a-b+c)) :
                                 (add ? o+(a-b+c) :   a-b+c )                
    _divergence!(kernel,div,u,v,S)
end

"""
    _divergence!(   kernel,
                    div::LowerTriangularMatrix,
                    u::LowerTriangularMatrix,
                    v::LowerTriangularMatrix,
                    S::SpectralTransform)

Generic divergence function of vector `u`,`v` that writes into the output into `div`.
Generic as it uses the kernel `kernel` such that curl, div, add or flipsign
options are provided through `kernel`, but otherwise a single function is used."""
function _divergence!(  kernel,
                        div::LowerTriangularMatrix{Complex{NF}},
                        u::LowerTriangularMatrix{Complex{NF}},
                        v::LowerTriangularMatrix{Complex{NF}},
                        S::SpectralTransform{NF}
                        ) where {NF<:AbstractFloat}

    @boundscheck size(u) == size(div) || throw(BoundsError)
    @boundscheck size(v) == size(div) || throw(BoundsError)

    @unpack grad_y_vordiv1,grad_y_vordiv2 = S
    @boundscheck size(grad_y_vordiv1) == size(div) || throw(BoundsError)
    @boundscheck size(grad_y_vordiv2) == size(div) || throw(BoundsError)
    lmax,mmax = size(div) .- (2,1)              # 0-based lmax,mmax 

    z = zero(Complex{NF})
    lm = 0
    @inbounds for m in 1:mmax+1                 # 1-based l,m
        
        # DIAGONAL (separate to avoid access to v[l-1,m])
        lm += 1                                 
        ∂u∂λ  = ((m-1)*im)*u[lm]
        ∂v∂θ1 = zero(Complex{NF})               # always above the diagonal
        ∂v∂θ2 = grad_y_vordiv2[lm]*v[lm+1]
        div[lm] = kernel(div[lm], ∂u∂λ, ∂v∂θ1, ∂v∂θ2)

        # BELOW DIAGONAL (but skip last row)
        for l in m+1:lmax+1
            lm += 1
            ∂u∂λ  = ((m-1)*im)*u[lm]
            ∂v∂θ1 = grad_y_vordiv1[lm]*v[lm-1]
            ∂v∂θ2 = grad_y_vordiv2[lm]*v[lm+1]  # this pulls in data from the last row though
            div[lm] = kernel(div[lm], ∂u∂λ, ∂v∂θ1, ∂v∂θ2)
        end

        # Last row
        lm += 1
    end

    return nothing
end

"""
    UV_from_vor!(   U::LowerTriangularMatrix,
                    V::LowerTriangularMatrix,
                    vor::LowerTriangularMatrix,
                    S::SpectralTransform)

Get U,V (=(u,v)*coslat) from vorticity ζ spectral space (divergence D=0)
Two operations are combined into a single linear operation. First, invert the
spherical Laplace ∇² operator to get stream function from vorticity. Then
compute zonal and meridional gradients to get U,V."""
function UV_from_vor!(  U::LowerTriangularMatrix{Complex{NF}},
                        V::LowerTriangularMatrix{Complex{NF}},
                        vor::LowerTriangularMatrix{Complex{NF}},
                        S::SpectralTransform{NF}
                        ) where {NF<:AbstractFloat}

    @unpack vordiv_to_uv_x,vordiv_to_uv1,vordiv_to_uv2 = S
    lmax,mmax = size(vor) .- (2,1)                  # 0-based lmax,mmax
    
    @boundscheck lmax == mmax || throw(BoundsError)
    @boundscheck size(U) == size(vor) || throw(BoundsError)
    @boundscheck size(V) == size(vor) || throw(BoundsError)
    @boundscheck size(vordiv_to_uv_x) == size(vor) || throw(BoundsError)
    @boundscheck size(vordiv_to_uv1) == size(vor) || throw(BoundsError)
    @boundscheck size(vordiv_to_uv2) == size(vor) || throw(BoundsError)

    lm = 0
    @inbounds for m in 1:mmax                       # 1-based l,m, exclude last column

        # DIAGONAL (separated to avoid access to l-1,m which is above the diagonal)
        lm += 1

        # U = -∂/∂lat(Ψ) and V = V = ∂/∂λ(Ψ) combined with Laplace inversion ∇⁻², omit radius R scaling
        U[lm] = vordiv_to_uv2[lm]*vor[lm+1]         # - vordiv_to_uv1[lm]*vor[l-1,m] <- is zero
        V[lm] = im*vordiv_to_uv_x[lm]*vor[lm]

        # BELOW DIAGONAL
        for l in m+1:lmax                           # skip last two rows
            lm += 1

            # U = -∂/∂lat(Ψ) and V = V = ∂/∂λ(Ψ) combined with Laplace inversion ∇⁻², omit radius R scaling
            # U[lm] = vordiv_to_uv2[lm]*vor[lm+1] - vordiv_to_uv1[lm]*vor[lm-1]
            U[lm] = muladd(vordiv_to_uv2[lm], vor[lm+1], -vordiv_to_uv1[lm]*vor[lm-1])
            V[lm] = im*vordiv_to_uv_x[lm]*vor[lm]
        end

        # SECOND LAST ROW
        lm += 1
        U[lm] = -vordiv_to_uv1[lm]*vor[lm-1]        # meridional gradient again (but only 2nd term from above)
        V[lm] = im*vordiv_to_uv_x[lm]*vor[lm]       # zonal gradient again (as above)

        # LAST ROW (separated to avoid out-of-bounds access to l+2,m)
        lm += 1
        U[lm] = -vordiv_to_uv1[lm]*vor[lm-1]        # meridional gradient again (but only 2nd term from above)
        V[lm] = zero(Complex{NF})                   # set explicitly to 0 as Ψ does not contribute to last row of V
    end

    # LAST COLUMN
    @inbounds begin
        lm += 1                     # second last row
        U[lm] = zero(Complex{NF})
        V[lm] = im*vordiv_to_uv_x[lm]*vor[lm]

        lm += 1                     # last row
        U[lm] = -vordiv_to_uv1[lm]*vor[lm-1]
        V[lm] = zero(Complex{NF})
    end

    return nothing
end

"""
    UV_from_vordiv!(U::LowerTriangularMatrix,
                    V::LowerTriangularMatrix,
                    vor::LowerTriangularMatrix,
                    div::LowerTriangularMatrix,
                    S::SpectralTransform)

Get U,V (=(u,v)*coslat) from vorticity ζ and divergence D in spectral space.
Two operations are combined into a single linear operation. First, invert the
spherical Laplace ∇² operator to get stream function from vorticity and
velocity potential from divergence. Then compute zonal and meridional gradients
to get U,V."""
function UV_from_vordiv!(   U::LowerTriangularMatrix{Complex{NF}},
                            V::LowerTriangularMatrix{Complex{NF}},
                            vor::LowerTriangularMatrix{Complex{NF}},
                            div::LowerTriangularMatrix{Complex{NF}},
                            S::SpectralTransform{NF}
                            ) where {NF<:AbstractFloat}

    @unpack vordiv_to_uv_x,vordiv_to_uv1,vordiv_to_uv2 = S
    lmax,mmax = size(vor) .- (2,1)                  # 0-based lmax,mmax
    @boundscheck lmax == mmax || throw(BoundsError)
    @boundscheck size(div) == size(vor) || throw(BoundsError)
    @boundscheck size(U) == size(vor) || throw(BoundsError)
    @boundscheck size(V) == size(vor) || throw(BoundsError)
    @boundscheck size(vordiv_to_uv_x) == size(vor) || throw(BoundsError)
    @boundscheck size(vordiv_to_uv1) == size(vor) || throw(BoundsError)
    @boundscheck size(vordiv_to_uv1) == size(vor) || throw(BoundsError)

    lm = 0
    @inbounds for m in 1:mmax                       # 1-based l,m, skip last column

        # DIAGONAL (separated to avoid access to l-1,m which is above the diagonal)
        lm += 1
        
        # div,vor contribution to meridional gradient
        ∂ζθ =  vordiv_to_uv2[lm]*vor[lm+1]          # lm-1 term is zero
        ∂Dθ = -vordiv_to_uv2[lm]*div[lm+1]          # lm-1 term is zero
        
        # the following is moved into the muladd        
        # ∂Dλ = im*vordiv_to_uv_x[lm]*div[lm]       # divergence contribution to zonal gradient
        # ∂ζλ = im*vordiv_to_uv_x[lm]*vor[lm]       # vorticity contribution to zonal gradient

        z = im*vordiv_to_uv_x[lm]
        U[lm] = muladd(z, div[lm], ∂ζθ)             # = ∂Dλ + ∂ζθ
        V[lm] = muladd(z, vor[lm], ∂Dθ)             # = ∂ζλ + ∂Dθ

        # BELOW DIAGONAL (all terms)
        for l in m+1:lmax                               # skip last row (lmax+2)
            lm += 1
            
            # div,vor contribution to meridional gradient
            # ∂ζθ = vordiv_to_uv2[lm]*vor[lm+1] - vordiv_to_uv1[lm]*vor[lm-1]
            # ∂Dθ = vordiv_to_uv1[lm]*div[lm-1] - vordiv_to_uv2[lm]*div[lm+1]
            ∂ζθ = muladd(vordiv_to_uv2[lm], vor[lm+1], -vordiv_to_uv1[lm]*vor[lm-1])
            ∂Dθ = muladd(vordiv_to_uv1[lm], div[lm-1], -vordiv_to_uv2[lm]*div[lm+1])

            # The following is moved into the muladd
            # ∂Dλ = im*vordiv_to_uv_x[lm]*div[lm]   # divergence contribution to zonal gradient
            # ∂ζλ = im*vordiv_to_uv_x[lm]*vor[lm]   # vorticity contribution to zonal gradient

            z = im*vordiv_to_uv_x[lm]
            U[lm] = muladd(z, div[lm], ∂ζθ)         # = ∂Dλ + ∂ζθ
            V[lm] = muladd(z, vor[lm], ∂Dθ)         # = ∂ζλ + ∂Dθ            
        end

        # SECOND LAST ROW (separated to imply that vor,div are zero in last row)
        lm += 1
        U[lm] = im*vordiv_to_uv_x[lm]*div[lm] - vordiv_to_uv1[lm]*vor[lm-1]
        V[lm] = im*vordiv_to_uv_x[lm]*vor[lm] + vordiv_to_uv1[lm]*div[lm-1]

        # LAST ROW (separated to avoid out-of-bounds access to lmax+3
        lm += 1
        U[lm] = -vordiv_to_uv1[lm]*vor[lm-1]        # only last term from 2nd last row
        V[lm] =  vordiv_to_uv1[lm]*div[lm-1]        # only last term from 2nd last row
    end

    # LAST COLUMN 
    @inbounds begin
        lm += 1                                     # second last row
        U[lm] = im*vordiv_to_uv_x[lm]*div[lm]       # other terms are zero
        V[lm] = im*vordiv_to_uv_x[lm]*vor[lm]       # other terms are zero

        lm += 1                                     # last row
        U[lm] = -vordiv_to_uv1[lm]*vor[lm-1]        # other terms are zero
        V[lm] =  vordiv_to_uv1[lm]*div[lm-1]        # other terms are zero
    end
end

"""
    ∇²!(    ∇²alms::LowerTriangularMatrix,
            alms::LowerTriangularMatrix,
            S::SpectralTransform;
            add::Bool=false,
            flipsign::Bool=false,
            inverse::Bool=false)

Laplace operator ∇² applied to the spectral coefficients `alms` in spherical
coordinates. The radius `R` is omitted in the eigenvalues which are precomputed in `S`.
∇²! is the in-place version which directly stores the output in the first argument `∇²alms`.

Options:
    - `add=true` adds the ∇²(alms) to the output
    - `flipsign=true` computes -∇²(alms) instead
    - `inverse=true` computes ∇⁻²(alms) instead

Default is `add=false`, `flipsign=false`, `inverse=false`. These options can be combined."""
function ∇²!(   ∇²alms::LowerTriangularMatrix{Complex{NF}}, # Output: (inverse) Laplacian of alms
                alms::LowerTriangularMatrix{Complex{NF}},   # Input: spectral coefficients
                S::SpectralTransform{NF};                   # precomputed eigenvalues
                add::Bool=false,                            # add to output array or overwrite
                flipsign::Bool=false,                       # -∇² or ∇²
                inverse::Bool=false,                        # ∇⁻² or ∇²
                ) where {NF<:AbstractFloat}

    @boundscheck size(alms) == size(∇²alms) || throw(BoundsError)
    lmax,mmax = size(alms) .- (1,1)     # 0-based degree l, order m of the Legendre polynomials
    
    # use eigenvalues⁻¹/eigenvalues for ∇⁻²/∇² based but name both eigenvalues
    eigenvalues = inverse ? S.eigenvalues⁻¹ : S.eigenvalues
    @boundscheck length(eigenvalues) >= lmax+1 || throw(BoundsError)

    @inline kernel(o,a) = flipsign ? (add ? (o-a) : -a)  :
                                     (add ? (o+a) :  a)

    lm = 0
    @inbounds for m in 1:mmax+1     # order m = 0:mmax but 1-based
        for l in m:lmax+1           # degree l = m:lmax but 1-based
            lm += 1
            ∇²alms[lm] = kernel(∇²alms[lm],alms[lm]*eigenvalues[l])
        end
    end

    return ∇²alms
end

"""
    ∇⁻²!(   ∇⁻²alms::LowerTriangularMatrix,
            alms::LowerTriangularMatrix,
            S::SpectralTransform;
            add::Bool=false,
            flipsign::Bool=false)

Calls `∇²!(∇⁻²alms, alms, S; add, flipsign, inverse=true)`."""
function ∇⁻²!(  ∇⁻²alms::LowerTriangularMatrix{Complex{NF}},# Output: inverse Laplacian of alms
                alms::LowerTriangularMatrix{Complex{NF}},   # Input: spectral coefficients
                S::SpectralTransform{NF};                   # precomputed eigenvalues
                add::Bool=false,                            # add to output array or overwrite
                flipsign::Bool=false,                       # -∇⁻² or ∇⁻²
                ) where {NF<:AbstractFloat}

    inverse = true
    return ∇²!(∇⁻²alms,alms,S;add,flipsign,inverse)
end

function ∇!(dpdx::LowerTriangularMatrix{Complex{NF}},       # Output: zonal gradient
            dpdy::LowerTriangularMatrix{Complex{NF}},       # Output: meridional gradient
            p::LowerTriangularMatrix{Complex{NF}},          # Input: spectral coefficients
            S::SpectralTransform{NF}                        # includes precomputed arrays
            ) where {NF<:AbstractFloat}

    lmax,mmax = size(p) .- (1,1)                            # 0-based, include last row
    @boundscheck size(p) == size(dpdx) || throw(BoundsError)
    @boundscheck size(p) == size(dpdy) || throw(BoundsError)

    @unpack grad_y1, grad_y2 = S

    lm = 0
    @inbounds for m in 1:mmax           # 1-based l,m, skip last column

        # DIAGONAL (separated to avoid access to l-1,m which is above the diagonal)
        lm += 1

        dpdx[lm] = (m-1)*im*p[lm]       # zonal gradient: d/dlon = *i*m
        dpdy[lm] = grad_y2[lm]*p[lm+1]  # meridional gradient: p[lm-1]=0 on diagonal
       
        # BELOW DIAGONAL (all terms)
        for l in m+1:lmax               # skip last row
            lm += 1
            
            dpdx[lm] = (m-1)*im*p[lm]
            dpdy[lm] = grad_y1[lm]*p[lm-1] + grad_y2[lm]*p[lm+1]
        end

        # LAST ROW (separated to avoid out-of-bounds access to lmax+2
        lm += 1
        dpdx[lm] = (m-1)*im*p[lm]
        dpdy[lm] = grad_y1[lm]*p[lm-1]  # only first term from 2nd last row
    end

    # LAST COLUMN
    @inbounds begin
        lm += 1 
        dpdx[lm] = mmax*im*p[lm]
        dpdy[lm] = grad_y2[lm]*p[lm+1]  # only 2nd term

        lm += 1
        dpdx[lm] = mmax*im*p[lm]
        dpdy[lm] = grad_y1[lm]*p[lm-1]  # only 1st term
    end

    return dpdx,dpdy
end