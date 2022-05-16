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
                            flipsign::Bool=false                    # flip the sign to obtain u from Ψ
                            ) where {NF<:AbstractFloat}             # number format NF

    _,mmax = size(Ψ)                # degree l, order m of spherical harmonics
    lmax, mmax = mmax-1, mmax-1     # convert to 0-based l,m, but use mmax for lmax/lmax+1 flexibility
    
    # u needs one more degree/meridional mode l for each m than Ψ due to the recursion
    # Ψ can have size n+1 x n but then the last row is not used in the loop
    size_same = size(coslat_u) == size(Ψ)
    size_compat = size_same || (size(coslat_u) .- (1,0)) == size(Ψ)
    @boundscheck size_compat || throw(BoundsError)

    if flipsign                     # used to get u from streamfunction Ψ (u ~ -∂Ψ/∂lat)
        grad_y1 = S.minus_grad_y1
        grad_y2 = S.minus_grad_y2
    else                            # no sign flip otherwise
        grad_y1 = S.grad_y1
        grad_y2 = S.grad_y2
    end

    coslat_u[1,1] = grad_y2[1,1]*Ψ[2,1]     # l=m=0 mode only with term 2

    @inbounds for m in 1:mmax
        for l in max(2,m):lmax
            coslat_u[l,m] = grad_y1[l,m]*Ψ[l-1,m] + grad_y2[l,m]*Ψ[l+1,m]
        end
        for l in lmax+1:lmax+2-size_same
            coslat_u[l,m] = grad_y1[l,m]*Ψ[l-1,m]
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
                                radius::Real=1                          # radius of the sphere/Earth
                                ) where {NF<:AbstractFloat}             # number format NF

    # Ψ can have size n+1 x n but then the last row is not used in the loop
    size_compat = size(coslat_v) == size(Ψ) || (size(coslat_v) .- (1,0)) == size(Ψ)
    @boundscheck size_compat || throw(BoundsError)
    lmax,mmax = size(Ψ) .- 1    # 0-based max degree l, order m of spherical harmonics

    iradius⁻¹ = convert(Complex{NF},im/radius)      # = imaginary/radius converted to NF

    @inbounds for m in 1:mmax+1                     # loop over all coefficients, order m
        for l in m:lmax+1                           # degree l
            coslat_v[l,m] = (m-1)*iradius⁻¹*Ψ[l,m]  # gradient in lon = *i*m/radius but 1-based order
        end
    end

    # grad in lon does not project onto the last degree l, set explicitly to zero in that case
    if size(Ψ) != size(coslat_v)
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
    lmax,mmax = size(Ψ)
    coslat_v = zeros(NF,lmax+one_more_l,mmax)       # preallocate output array (gradient in longitude)
    return gradient_longitude!(coslat_v,Ψ,R)        # call in-place version
end

"""Divide a gridded field `A` by the cosine of latitude."""
function unscale_coslat!(   A::AbstractMatrix{NF},
                            G::Geometry{NF}) where NF
    
    nlon, nlat = size(A)
    @boundscheck nlat == G.nlat || throw(BoundsError)

    @unpack coslat⁻¹ = G

    @inbounds for j in 1:nlat
        for i in 1:nlon
            A[i,j] *= coslat⁻¹[j]
        end
    end
end

"""Multiply a gridded field `A` by the cosine of latitude."""
function scale_coslat!( A::AbstractMatrix{NF},
                        G::Geometry{NF}) where NF
    
    nlon, nlat = size(A)
    @boundscheck nlat == G.nlat || throw(BoundsError)

    @unpack coslat = G

    @inbounds for j in 1:nlat
        for i in 1:nlon
            A[i,j] *= coslat[j]
        end
    end
end

"""
    ∇⁻²!(   ∇⁻²alms::AbstractMatrix{Complex},
            alms::AbstractMatrix{Complex},
            R::Real=1)

Inverse spherical Laplace operator ∇⁻² applied to the spectral coefficients `alms` on a sphere
of radius `R`. ∇⁻²! is the in-place version which directly stores the output in the argument
∇⁻²alms. The integration constant for Legendre polynomial `l=m=0` is zero. The inverse spherical
Laplace operator is

    ∇⁻²alms = alms*R²/(-l(l+1))

with the degree `l` (0-based) of the Legendre polynomial."""
function ∇⁻²!(  ∇⁻²alms::AbstractMatrix{Complex{NF}},   # Output: inverse Laplacian of alms
                alms::AbstractMatrix{Complex{NF}},      # spectral coefficients
                S::SpectralTransform{NF}                # precomputed arrays for spectral space
                ) where {NF<:AbstractFloat}

    @boundscheck size(alms) == size(∇⁻²alms) || throw(BoundsError)
    lmax,mmax = size(alms) .- 1     # degree l, order m of the Legendre polynomials
    @unpack eigen_values⁻¹ = S

    @inbounds for m in 1:mmax+1     # order m = 0:mmax but 1-based
        for l in m:lmax+1           # degree l = 0:lmax but 1-based
            ∇⁻²alms[l,m] = alms[l,m]*eigen_values⁻¹[l]
        end
    end

    return ∇⁻²alms
end