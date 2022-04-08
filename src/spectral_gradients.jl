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
function gradient_latitude!(coslat_u::AbstractMatrix{Complex{NF}},   # output: cos(lat)*zonal velocity u
                            Ψ::AbstractMatrix{Complex{NF}},          # input: streamfunction Ψ
                            ϵlms::AbstractMatrix{NF},                # recursion factors
                            R::Real=1                               # radius of the sphere/Earth
                            ) where {NF<:AbstractFloat}             # number format NF

    _,mmax = size(Ψ)                # degree l, order m of spherical harmonics
    lmax, mmax = mmax-1, mmax-1     # convert to 0-based l,m, but use mmax for lmax/lmax+1 flexibility
    
    # u needs one more degree/meridional mode l for each m than Ψ due to the recursion
    # Ψ can have size n+1 x n but then the last row is not used in the loop
    size_compat = size(coslat_u) == size(Ψ) || (size(coslat_u) .- (1,0)) == size(Ψ)
    @boundscheck size_compat || throw(BoundsError)
    R⁻¹ = convert(Complex{NF},1/R)                              # 1/radius of the sphere

    # for loops implement the recursion formula (0-based degree l, order m)
    # (coslat*u)_lm = -1/R(  -(l-1)*ϵ(l,m)  *Ψ_(l-1,m)          # recursion term 1
    #                       (l+2)*ϵ(l+1,m)*Ψ_(l+1,m))           # recursion term 2

    # convert to 1-based l,m
    @inbounds for m in 1:mmax         # exclude m=mmax+1 as for m=l=mmax+1 term1=term2=0

        # 1. Diagonal (l=m), term 1 = 0 for the l=m modes
        # coslat_u[m,m] = -R⁻¹*(m+1)*ϵlms[m+1,m]*Ψ[m+1,m]         # recursion term 2 only
        coslat_u[m,m] = (1-m)*ϵlms[m+1,m]*Ψ[m+1,m]              # recursion term 2 only

        # 2. Below diagonal modes (l>m, but l<=lmax)
        for l in m+1:lmax
            # coslat_u[l,m] = -R⁻¹*(-(l-2)*ϵlms[l  ,m]*Ψ[l-1,m] + # term 1
            #                        (l+1)*ϵlms[l+1,m]*Ψ[l+1,m])  # term 2
            coslat_u[l,m] = (l-2)*ϵlms[l  ,m]*Ψ[l-1,m] -        # term 1
                            (l+1)*ϵlms[l+1,m]*Ψ[l+1,m]          # term 2
        end

        # 3. Last two rows
        for l in lmax+1:lmax+2                                  # recursion term 1 only
            # coslat_u[l,m] = R⁻¹*(l-2)*ϵlms[l,m]*Ψ[l-1,m]
            coslat_u[l,m] = (l-2)*ϵlms[l,m]*Ψ[l-1,m]
        end
    end
    # coslat_u[end,end] not needed as ϵ(l=m) = 0 in this case (l=lmax+1,m=mmax)

    return coslat_u
end

"""gradient_latitude! but precalculate the recursion factors `ϵlms` in case they are not provided."""
function gradient_latitude!(coslat_u::AbstractArray{Complex{NF}},   # output: cos(lat)*u
                            Ψ::AbstractArray{Complex{NF}},          # input: streamfunction Ψ
                            R::Real=1                               # radius of the sphere/Earth
                            ) where {NF<:AbstractFloat}             # number format NF
    _,mmax = size(Ψ) .- 1                                           # degree l, order m of spherical harmonics   
    ϵlms = get_recursion_factors(NF,mmax,mmax)                      # precalculate recursion factors
    return gradient_latitude!(coslat_u,Ψ,ϵlms,R)                    # call in-place function
end

function gradient_latitude( Ψ::AbstractArray{Complex{NF}},  # input: streamfunction Ψ
                            R::Real=1                       # radius of the sphere/Earth
                            ) where {NF<:AbstractFloat}     # number format NF
    _,mmax = size(Ψ) .- 1                                   # max degree l, order m of spherical harmonics
    coslat_u = zeros(Complex{NF},mmax+2,mmax+1)             # preallocate output, one more l for recursion
    return gradient_latitude!(coslat_u,Ψ,R)                 # call in-place version
end

"""
    coslat_∂alms_∂lon = gradient_longitude!(  coslat_∂alms_∂lon::AbstractMatrix{Complex{NF}},
                                            alms::AbstractMatrix{Complex{NF}};
                                            R::Real=1
                                            ) where {NF<:AbstractFloat}

Zonal gradient in spectral space of spherical harmonic coefficients `alms` on a sphere with radius `R`.
While the zonal gradient is 1/coslat*∂/∂lon in spherical coordinates, this functions omits the 1/coslat scaling
such that the return array is coslat*∂alms/∂lon.
"""
function gradient_longitude!(   coslat_∂alms_∂lon::AbstractMatrix{Complex{NF}}, # output: coslat*zonal gradient
                                alms::AbstractMatrix{Complex{NF}},          # input: spectral coefficients
                                R::Real=1                                   # radius of the sphere/Earth
                                ) where {NF<:AbstractFloat}                 # number format NF

    @boundscheck size(alms) == size(coslat_∂alms_∂lon) || throw(BoundsError)
    lmax,mmax = size(alms) .- 1                         # maximum degree l, order m of spherical harmonics
    iR⁻¹ = convert(Complex{NF},im/R)                    # = imaginary/radius converted to NF

    @inbounds for m in 1:mmax+1                         # loop over all coefficients, order m
        for l in m:lmax+1                               # degree l
            coslat_∂alms_∂lon[l,m] = (m-1)*iR⁻¹*alms[l,m] # gradient in lon = *i*m/R but order m is 1-based
        end
    end

    return coslat_∂alms_∂lon
end

"""Gradient in longitude in spectral space. Input: coefficients `alms` of the spherical harmonics."""
function gradient_longitude(alms::AbstractMatrix{Complex{NF}},  # input array: spectral coefficients
                            R::Real=1                           # radius of the sphere/Earth
                            ) where NF                          # number format NF
    ∂alms_∂lon = zero(alms)                             # preallocate output array (gradient in longitude)
    return gradient_longitude!(∂alms_∂lon,alms,R)       # call in-place version
end

"""Spectral tendency of ∇⋅(uv*ω) from vector uv=(u,v) in grid space and absolute vorticity ω.
Step 1 (grid space): Add Coriolis f to the relative vorticity ζ (=`vor_grid`) to obtain abs vorticity ω.
Step 2 (grid space): Multiply u,v with abs vorticity ω.
Step 3 (grid space): Unscale with coslat, cosine of latitude, as the gradients will include a coslat term.
Step 4 (spectral space): convert uω/coslat, vω/coslat from grid to spectral space
Step 5 (spectral space): Compute gradients ∂/∂lon(uω/coslat) and ∂/∂lat(vω/coslat)
Step 6 (spectral space): Add ∂/∂lon(uω/coslat)+∂/∂θ(vω/coslat) and return.
"""
function divergence_uvω_spectral(   u_grid::AbstractMatrix{NF},     # zonal velocity in grid space
                                    v_grid::AbstractMatrix{NF},     # meridional velocity in grid space
                                    vor_grid::AbstractMatrix{NF},   # relative vorticity in grid space       
                                    G::GeoSpectral{NF}              # struct with geometry and spectral transform
                                    ) where {NF<:AbstractFloat}

    nlon,nlat = size(u_grid)
    @boundscheck size(u_grid) == size(v_grid) || throw(BoundsError)

    @unpack f_coriolis,coslat⁻¹ = G.geometry

    uω_grid_coslat⁻¹ = zero(u_grid)                             # TODO preallocate elsewhere
    vω_grid_coslat⁻¹ = zero(v_grid)

    @inbounds for j in 1:nlat
        for i in 1:nlon
            ω = vor_grid[i,j] + f_coriolis[j]                   # = relative vorticity + coriolis
            uω_grid_coslat⁻¹[i,j] = ω*u_grid[i,j]*coslat⁻¹[j]   # = u(vor+f)/cos(ϕ)
            vω_grid_coslat⁻¹[i,j] = ω*v_grid[i,j]*coslat⁻¹[j]   # = v(vor+f)/cos(ϕ)
        end
    end

    # TODO preallocate returned coefficients elsewhere
    uω_coslat⁻¹ = spectral(uω_grid_coslat⁻¹,G.spectral,one_more_l=true)         
    vω_coslat⁻¹ = spectral(vω_grid_coslat⁻¹,G.spectral,one_more_l=true)
    
    ∂uω_∂lon = gradient_longitude(uω_coslat⁻¹)                  # spectral gradients
    ∂vω_∂lat = gradient_latitude(vω_coslat⁻¹)

    return ∂uω_∂lon + ∂vω_∂lat                                  # add for divergence
end


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

    lmax,mmax = size(alms) .- 1             # degree l, order m of the Legendre polynomials
    R_inv = convert(Complex{NF},inv(R))     # =1/R, 1 over radius

    @inbounds for m in 1:mmax+1             # order m = 0:mmax but 1-based
        for l in m:lmax+1                   # degree l = 0:lmax but 1-based
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
                ) where NF                              # number format NF

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
function ∇²(alms::AbstractMatrix{Complex{NF}},  # spectral coefficients
            R::Real=1) where NF                 # radius of the sphere/Earth
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
function ∇⁴!(   ∇⁴alms::AbstractMatrix{Complex{NF}},    # Output: Bi-Laplacian of alms
                alms::AbstractMatrix{Complex{NF}},      # spectral coefficients
                R::Real=1) where NF
    
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
function ∇⁴(alms::AbstractMatrix{Complex{NF}},  # spectral coefficients
            R::Real=1) where NF                 # radius of the Earth
    
    ∇⁴alms = copy(alms)                     # allocate output array
    return R == 1 ? ∇⁴!(∇⁴alms,alms) : ∇⁴!(∇⁴alms,alms,R)
end

"""
    ∇⁻²!(   ∇⁻²alms::AbstractMatrix{Complex},       # Out: inverse Laplace of alms
            alms::AbstractMatrix{Complex})          # In: spectral coefficients alms

Inverse spherical Laplace operator ∇⁻² applied to the spectral coefficients `alms` on a sphere
of radius `R=1`. ∇⁻²! is the in-place version which directly stores the output in the argument
∇⁻²alms. The integration constant for Legendre polynomial `l=m=0` is zero."""
function ∇⁻²!(  ∇⁻²alms::AbstractMatrix{Complex{NF}},   # Output: inverse Laplacian of alms
                alms::AbstractMatrix{Complex{NF}}       # spectral coefficients
                ) where NF      

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
function ∇⁻²(   alms::AbstractMatrix{Complex{NF}},  # spectral coefficients
                R::Real=1) where NF                 # radius of the Earth
    ∇⁻²alms = copy(alms)                            # preallocate output
    return R == 1 ? ∇⁻²!(∇⁻²alms,alms) : ∇⁻²!(∇⁻²alms,alms,R)
end