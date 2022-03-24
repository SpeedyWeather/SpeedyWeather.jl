# function SpectralTransform(P::Parameters,G::Geometry)

#     @unpack nlon, nlat, nlon_half, nlat_half = G
#     @unpack coslat_NH = G
#     @unpack R_earth, trunc = P

#     # SIZE OF SPECTRAL GRID
#     mx = trunc+1
#     nx = trunc+1

#     # PLAN THE FFTs
#     # rfft_plan = plan_rfft(rand(NF,nlon))
#     # irfft_plan = plan_irfft(rand(Complex{NF},nlon_half+1),nlon)

#     # LEGENDRE WEIGHTS from pole to equator (=first half or array)
#     leg_weight = FastGaussQuadrature.gausslegendre(nlat)[2][1:nlat_half]

#     # Spectral packing of speedy is m',n', where m' = m, n' = m+n with m,n being the
#     # conventional wavenumbers. Due to Julia's 1-based indexing subtract two, as
#     # the Legendre polynomials start with
#     nsh2 = zeros(Int, nx)
#     for n in 1:nx
#         for m in 1:mx
#             if m + n - 2 <= trunc + 1
#                 nsh2[n] = nsh2[n] + 1
#             end
#         end
#     end

#     # Epsilon-factors for the recurrence relation of the associated Legendre polynomials.
#     ε,ε⁻¹ = ε_recurrence(mx,nx)

#     # Generate associated Legendre polynomials
#     # get_legendre_poly computes the polynomials at a particular latitiude
#     leg_poly = zeros(mx, nx, nlat_half)
#     for j in 1:nlat_half
#         leg_poly[:,:,j] = legendre_polynomials(j,ε,ε⁻¹,mx,nx,G)
#     end

#     # leg_poly = zeros(mx, nx, nlat_half)
#     # for j in 1:nlat_half
#     #     leg_poly[:,:,j] = AssociatedLegendrePolynomials.λlm(0:mx-1,0:nx-1,coslat_NH[j])
#     # end

#     # LAPLACIANS for harmonic & biharmonic diffusion
#     ∇²,∇⁻²,∇⁴ = Laplacians(mx,nx,R_earth)

#     gradx   = zeros(mx)          #TODO what's this?
#     uvdx    = zeros(mx, nx)
#     uvdym   = zeros(mx, nx)
#     uvdyp   = zeros(mx, nx)
#     gradym  = zeros(mx, nx)
#     gradyp  = zeros(mx, nx)
#     vddym   = zeros(mx, nx)
#     vddyp   = zeros(mx, nx)

#     for m in 1:mx
#         for n in 1:nx
#             m1  = m - 1
#             m2  = m1 + 1
#             el1 = m + n - 2
#             if n == 1
#                 gradx[m]   = m1/R_earth
#                 uvdx[m,1]  = -R_earth/(m1 + 1)
#                 uvdym[m,1] = 0.0
#                 vddym[m,1] = 0.0
#             else
#                 uvdx[m,n]   = -R_earth*m1/(el1*(el1 + 1.0))
#                 gradym[m,n] = (el1 - 1.0)*ε[m2,n]/R_earth
#                 uvdym[m,n]  = -R_earth*ε[m2,n]/el1
#                 vddym[m,n]  = (el1 + 1.0)*ε[m2,n]/R_earth
#             end
#             gradyp[m,n] = (el1 + 2.0)*ε[m2,n+1]/R_earth
#             uvdyp[m,n]  = -R_earth*ε[m2,n+1]/(el1 + 1.0)
#             vddyp[m,n]  = el1*ε[m2,n+1]/R_earth
#         end
#     end

#     SpectralTransform{P.NF}(    trunc,mx,nx,
#                                 # rfft_plan,irfft_plan,
#                                 leg_weight,nsh2,leg_poly,∇²,∇⁻²,∇⁴,
#                                 gradx,uvdx,uvdym,uvdyp,gradym,gradyp,vddym,vddyp)
# end


# function grad!( ψ::Array{NF,2},
#                 psdx::Array{Complex{NF},2},
#                 psdy::Array{NF,2},
#                 G::GeoSpectral{NF}) where {NF<:AbstractFloat}

#     #TODO boundscheck

#     @unpack trunc, mx, nx = G.spectral
#     @unpack gradx, gradyp, gradym = G.spectral

#     for n in 1:nx
#         psdx[:,n] = gradx.*ψ[:,n]*im
#     end

#     for m in 1:mx
#         psdy[m,1]  =  gradyp[m,1]*ψ[m,2]
#         psdy[m,nx] = -gradym[m,nx]*ψ[m,trunc+1]
#     end

#     for n in 2:trunc+1
#         for m in 1:mx
#             psdy[m,n] = -gradym[m,n]*ψ[m,n-1] + gradyp[m,n]*ψ[m,n+1]
#         end
#     end
# end

# function vds!(  ucosm::Array{NF,2},
#                 vcosm::Array{NF,2},
#                 vorm::Array{NF,2},
#                 divm::Array{NF,2},
#                 G::GeoSpectral{NF}) where {NF<:AbstractFloat}

#     #TODO boundscheck

#     @unpack trunc, mx, nx = G.spectral
#     @unpack gradx, vddym, vddyp = G.spectral

#     #TODO preallocate in a diagnosticvars struct
#     zp = zeros(Complex{NF}, mx,nx)
#     zc = zeros(Complex{NF}, mx,nx)

#     for n in 1:nx
#         zp[:,n] = gradx.*ucosm[:,n]*im
#         zc[:,n] = gradx.*vcosm[:,n]*im
#     end

#     for m in 1:mx
#         #TODO this has an implicit conversion to complex{NF}, issue?
#         vorm[m,1]  = zc[m,1] - vddyp[m,1]*ucosm[m,2]
#         vorm[m,nx] = vddym[m,nx]*ucosm[m,trunc+1]
#         divm[m,1]  = zp[m,1] + vddyp[m,1]*vcosm[m,2]
#         divm[m,nx] = -vddym[m,nx]*vcosm[m,trunc+1]
#     end

#     for n in 2:trunc+1
#         for m in 1:mx
#             #TODO same here
#             vorm[m,n] =  vddym[m,n]*ucosm[m,n-1] - vddyp[m,n]*ucosm[m,n+1] + zc[m,n]
#             divm[m,n] = -vddym[m,n]*vcosm[m,n-1] + vddyp[m,n]*vcosm[m,n+1] + zp[m,n]
#         end
#     end
# end

# function uvspec!(   vorm::Array{NF,2},
#                     divm::Array{NF,2},
#                     ucosm::Array{NF,2},
#                     vcosm::Array{NF,2},
#                     G::GeoSpectral{NF}) where {NF<:AbstractFloat}

#     #TODO boundscheck

#     @unpack trunc, mx, nx = G.spectral
#     @unpack uvdx, uvdyp, uvdym = G.spectral

#     #TODO preallocate elsewhere
#     zp = uvdx.*vorm*im
#     zc = uvdx.*divm*im

#     for m in 1:mx
#         ucosm[m,1]  =  zc[m,1] - uvdyp[m,1]*vorm[m,2]
#         ucosm[m,nx] =  uvdym[m,nx]*vorm[m,trunc+1]
#         vcosm[m,1]  =  zp[m,1] + uvdyp[m,1]*divm[m,2]
#         vcosm[m,nx] = -uvdym[m,nx]*divm[m,trunc+1]
#     end

#     for n in 2:trunc+1
#         for m in 1:mx
#           vcosm[m,n] = -uvdym[m,n]*divm[m,n-1] + uvdyp[m,n]*divm[m,n+1] + zp[m,n]
#           ucosm[m,n] =  uvdym[m,n]*vorm[m,n-1] - uvdyp[m,n]*vorm[m,n+1] + zc[m,n]
#         end
#     end
# end

# function vdspec!(   ug::Array{NF,2},
#                     vg::Array{NF,2},
#                     vorm::Array{NF,2},
#                     divm::Array{NF,2},
#                     kcos::Bool,
#                     G::GeoSpectral{NF}) where {NF<:AbstractFloat}

#     #TODO boundscheck

#     @unpack nlat, nlon, cosgr, cosgr2 = G.geometry

#     #TODO preallocate elsewhere
#     ug1 = zeros(NF, nlon, nlat)
#     vg1 = zeros(NF, nlon, nlat)

#     # either cosgr or cosgr2
#     cosgr = kcos ? cosgr : cosgr2

#     for j in 1:nlat
#         for i in 1:nlon
#             ug1[i,j] = ug[i,j]*cosgr[j]
#             vg1[i,j] = vg[i,j]*cosgr[j]
#         end
#     end

#     #TODO add spectral_trans and geometry as arguments
#     specu = spectral(ug1,G)
#     specv = spectral(vg1,G)
#     vds!(specu, specv, vorm, divm)
# end


"""
    ∇²!(∇²alms::AbstractMatrix{Complex{NF}},    # Output: Laplacian of alms
        alms::AbstractMatrix{Complex{NF}},      # spectral coefficients
        R::Real                                 # radius of the Earth
        ) where {NF<:AbstractFloat}             # number format NF

Spherical Laplace operator ∇² applied to the spectral coefficients `alms` on a sphere
of radius `R`. The spherical Laplace operator is

    ∇²alms = -l(l+1)/R²*alms

with the degree `l` of the Legendre polynomial. ∇²! is the in-place version, storing
the result directly in ∇²alms."""
function ∇²!(   ∇²alms::AbstractMatrix{Complex{NF}},    # Output: Laplacian of alms
                alms::AbstractMatrix{Complex{NF}},      # spectral coefficients
                R::Real                                 # radius of the sphere/Earth
                ) where {NF<:AbstractFloat}             # number format NF

    @boundscheck size(alms) == size(∇²alms) || throw(BoundsError)

    lmax,mmax = size(alms) .- 1         # degree l, order m of the Legendre polynomials
    R_inv = convert(Complex{NF},inv(R))  # =1/R, 1 over radius

    @inbounds for m in 1:mmax+1     # order m = 0:mmax but 1-based
        for l in m:lmax+1           # degree l = 0:lmax but 1-based
            # ∇²alms = -l(l+1)/R²*alms, but 1-based (l=>l-1)
            # R⁻² is split to avoid under/overflows
            ∇²alms[l,m] = ((1-l)*R_inv)*(l*R_inv)*alms[l,m]
        end
    end
    return ∇²alms
end

"""
    ∇²!(∇²alms::AbstractMatrix{Complex},    # Output: Laplacian of alms
        alms::AbstractMatrix{Complex})      # spectral coefficients

Same as ∇²!(::AbstractMatrix,::AbstractMatrix,R::Real) but assuming a sphere of radius `R=1`.
"""
function ∇²!(   ∇²alms::AbstractMatrix{Complex{NF}},    # Output: Laplacian of alms
                alms::AbstractMatrix{Complex{NF}}       # spectral coefficients
            ) where NF                                  # number format

    @boundscheck size(alms) == size(∇²alms) || throw(BoundsError)

    lmax,mmax = size(alms) .- 1     # degree l, order m of the Legendre polynomials

    @inbounds for m in 1:mmax+1     # order m = 0:mmax but 1-based
        for l in m:lmax+1           # degree l = 0:lmax but 1-based
            # ∇²alms = -l(l+1)/R²*alms, but 1-based (l=>l-1) and R=1
            ∇²alms[l,m] = (l*(1-l))*alms[l,m]
        end
    end
    return ∇²alms
end

"""
    ∇²alms = ∇²(alms::AbstractMatrix{Complex},R::Real=1)

Spherical Laplace operator ∇² applied to the spectral coefficients `alms` on a sphere
of radius `R`. ∇²(alms,R) is the non-in-place version of ∇²! which first allocates the
output array ∇²alms before calling ∇²!."""
function ∇²(alms::AbstractMatrix{Complex},      # spectral coefficients
            R::Real=1)                          # radius of the sphere/Earth
    ∇²alms = copy(alms)                         # allocate output
    return R == 1 ? ∇²!(∇²alms,alms) : ∇²!(∇²alms,alms,R)
end

"""
    ∇⁴!(∇⁴alms::AbstractMatrix{Complex},    # Output: Laplacian of alms
        alms::AbstractMatrix{Complex},      # spectral coefficients
        R::Real=1)                          # radius of the sphere/Earth

Spherical Bi-Laplace operator ∇⁴ = ∇²(∇²) applied to the spectral coefficients `alms` on a sphere
of radius `R`. ∇⁴! operates by applying ∇²! twice. ∇⁴! is the in-place version, storing
the result directly in the argument ∇⁴alms."""
function ∇⁴!(   ∇⁴alms::AbstractMatrix{Complex},     # Output: Bi-Laplacian of alms
                alms::AbstractMatrix{Complex},       # spectral coefficients
                R::Real=1)
    
    if R == 1                       # execute the non-R version
        ∇²!(∇⁴alms,alms)            # apply first Laplacian
        ∇²!(∇⁴alms,∇⁴alms)          # apply 2nd Laplacian
    else                            # scale by 1/R²
        ∇²!(∇⁴alms,alms,R)          # apply first Laplacian
        ∇²!(∇⁴alms,∇⁴alms,R)        # apply 2nd Laplacian
    end
    return ∇⁴alms
end

"""
    ∇⁴alms = ∇⁴(alms::AbstractMatrix{Complex},R::Real=1)

Spherical Bi-Laplace operator ∇⁴ = ∇²(∇²) applied to the spectral coefficients `alms` on a sphere
of radius `R`. ∇⁴ operates by applying ∇² twice. ∇⁴ is the non-in-place version of ∇⁴! which first
allocates the output array ∇⁴alms before calling ∇⁴!."""
function ∇⁴(alms::AbstractMatrix{Complex},  # spectral coefficients
            R::Real=1)                      # radius of the Earth
    
    ∇⁴alms = copy(alms)                     # allocate output array
    return R == 1 ? ∇⁴!(∇⁴alms,alms) : ∇⁴!(∇⁴alms,alms,R)
end

"""
    ∇⁻²!(   ∇⁻²alms::AbstractMatrix{Complex},       # Out: inverse Laplace of alms
            alms::AbstractMatrix{Complex})          # In: spectral coefficients alms

Inverse spherical Laplace operator ∇⁻² applied to the spectral coefficients `alms` on a sphere
of radius `R=1`. ∇⁻²! is the in-place version which directly stores the output in the argument
∇⁻²alms. The integration constant for Legendre polynomial `l=m=0` is zero."""
function ∇⁻²!(  ∇⁻²alms::AbstractMatrix{Complex},   # Output: inverse Laplacian of alms
                alms::AbstractMatrix{Complex})      # spectral coefficients

    @boundscheck size(alms) == size(∇⁻²alms) || throw(BoundsError)
    lmax,mmax = size(alms) .- 1     # degree l, order m of the Legendre polynomials

    @inbounds for m in 1:mmax+1     # order m = 0:mmax but 1-based
        for l in m:lmax+1           # degree l = 0:lmax but 1-based
            # ∇⁻²alms = R²/(-l(l+1))*alms, but 1-based (l=>l-1) and R=1
            ∇⁻²alms[l,m] = alms[l,m]/(l*(1-l))
        end
    end

    # set the integration constant (l=m=0 polynomial) to zero 
    ∇⁻²alms[1,1] = zero(∇⁻²alms[1,1]) 
    return ∇⁻²alms
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
                R::Real                                 # radius of the sphere/Earth
                ) where {NF<:AbstractFloat}

    @boundscheck size(alms) == size(∇⁻²alms) || throw(BoundsError)
    lmax,mmax = size(alms) .- 1     # degree l, order m of the Legendre polynomials
    R² = convert(NF,R^2)

    @inbounds for m in 1:mmax+1     # order m = 0:mmax but 1-based
        for l in m:lmax+1           # degree l = 0:lmax but 1-based
            # ∇²alms = R²/(-l(l+1))*alms, but 1-based and R=1
            ∇⁻²alms[l,m] = alms[l,m]/(l*(1-l))*R²
        end
    end

    # set integration constant for Legendre polynomial l=m=0 (0-based index) to zero
    ∇⁻²alms[1,1] = zero(Complex{NF})
    return ∇⁻²alms
end

"""
    ∇⁻²alms = ∇⁻²(alms::AbstractMatrix{Complex},R::Real=1)

Inverse spherical Laplace operator ∇⁻² applied to the spectral coefficients `alms` on a sphere
of radius `R`. ∇⁻² is the non-in-place version of ∇⁻²! and therefore first allocates the output array
∇⁻²alms before calling ∇⁻²!. The integration constant for Legendre polynomial `l=m=0` is zero."""
function ∇⁻²(   alms::AbstractMatrix{Complex},  # spectral coefficients
                R::Real=1)                      # radius of the Earth
    ∇⁻²alms = copy(alms)                        # preallocate output
    return R == 1 ? ∇⁻²!(∇⁻²alms,alms) : ∇⁻²!(∇⁻²alms,alms,R)
end