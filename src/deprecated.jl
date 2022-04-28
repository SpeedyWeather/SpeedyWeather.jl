# """
#     ∇²!(∇²alms::AbstractMatrix{Complex{NF}},    # Output: Laplacian of alms
#         alms::AbstractMatrix{Complex{NF}},      # spectral coefficients
#         R::Real                                 # radius of the Earth
#         ) where {NF<:AbstractFloat}             # number format NF

# Spherical Laplace operator ∇² applied to the spectral coefficients `alms` on a sphere
# of radius `R`. The spherical Laplace operator is

#     ∇²alms = -l(l+1)/R²*alms

# with the degree `l` of the Legendre polynomial. ∇²! is the in-place version, storing
# the result directly in ∇²alms."""
# function ∇²!(   ∇²alms::AbstractMatrix{Complex{NF}},    # Output: Laplacian of alms
#                 alms::AbstractMatrix{Complex{NF}},      # spectral coefficients
#                 R::Real                                 # radius of the sphere/Earth
#                 ) where {NF<:AbstractFloat}             # number format NF

#     @boundscheck size(alms) == size(∇²alms) || throw(BoundsError)

#     lmax,mmax = size(alms) .- 1             # degree l, order m of the Legendre polynomials
#     R_inv = convert(Complex{NF},inv(R))     # =1/R, 1 over radius

#     @inbounds for m in 1:mmax+1             # order m = 0:mmax but 1-based
#         for l in m:lmax+1                   # degree l = 0:lmax but 1-based
#             # ∇²alms = -l(l+1)/R²*alms, but 1-based (l=>l-1)
#             # R⁻² is split to avoid under/overflows
#             ∇²alms[l,m] = ((1-l)*R_inv)*(l*R_inv)*alms[l,m]
#         end
#     end
#     return ∇²alms
# end

# """
#     ∇²!(∇²alms::AbstractMatrix{Complex},    # Output: Laplacian of alms
#         alms::AbstractMatrix{Complex})      # spectral coefficients

# Same as ∇²!(::AbstractMatrix,::AbstractMatrix,R::Real) but assuming a sphere of radius `R=1`.
# """
# function ∇²!(   ∇²alms::AbstractMatrix{Complex{NF}},    # Output: Laplacian of alms
#                 alms::AbstractMatrix{Complex{NF}}       # spectral coefficients
#                 ) where NF                              # number format NF

#     @boundscheck size(alms) == size(∇²alms) || throw(BoundsError)

#     lmax,mmax = size(alms) .- 1     # degree l, order m of the Legendre polynomials

#     @inbounds for m in 1:mmax+1     # order m = 0:mmax but 1-based
#         for l in m:lmax+1           # degree l = 0:lmax but 1-based
#             # ∇²alms = -l(l+1)/R²*alms, but 1-based (l=>l-1) and R=1
#             ∇²alms[l,m] = (l*(1-l))*alms[l,m]
#         end
#     end
#     return ∇²alms
# end

# """
#     ∇²alms = ∇²(alms::AbstractMatrix{Complex},R::Real=1)

# Spherical Laplace operator ∇² applied to the spectral coefficients `alms` on a sphere
# of radius `R`. ∇²(alms,R) is the non-in-place version of ∇²! which first allocates the
# output array ∇²alms before calling ∇²!."""
# function ∇²(alms::AbstractMatrix{Complex{NF}},  # spectral coefficients
#             R::Real=1) where NF                 # radius of the sphere/Earth
#     ∇²alms = copy(alms)                         # allocate output
#     return R == 1 ? ∇²!(∇²alms,alms) : ∇²!(∇²alms,alms,R)
# end

# """
#     ∇⁴!(∇⁴alms::AbstractMatrix{Complex},    # Output: Laplacian of alms
#         alms::AbstractMatrix{Complex},      # spectral coefficients
#         R::Real=1)                          # radius of the sphere/Earth

# Spherical Bi-Laplace operator ∇⁴ = ∇²(∇²) applied to the spectral coefficients `alms` on a sphere
# of radius `R`. ∇⁴! operates by applying ∇²! twice. ∇⁴! is the in-place version, storing
# the result directly in the argument ∇⁴alms."""
# function ∇⁴!(   ∇⁴alms::AbstractMatrix{Complex{NF}},    # Output: Bi-Laplacian of alms
#                 alms::AbstractMatrix{Complex{NF}},      # spectral coefficients
#                 R::Real=1) where NF
    
#     if R == 1                       # execute the non-R version
#         ∇²!(∇⁴alms,alms)            # apply first Laplacian
#         ∇²!(∇⁴alms,∇⁴alms)          # apply 2nd Laplacian
#     else                            # scale by 1/R²
#         ∇²!(∇⁴alms,alms,R)          # apply first Laplacian
#         ∇²!(∇⁴alms,∇⁴alms,R)        # apply 2nd Laplacian
#     end
#     return ∇⁴alms
# end

# """
#     ∇⁴alms = ∇⁴(alms::AbstractMatrix{Complex},R::Real=1)

# Spherical Bi-Laplace operator ∇⁴ = ∇²(∇²) applied to the spectral coefficients `alms` on a sphere
# of radius `R`. ∇⁴ operates by applying ∇² twice. ∇⁴ is the non-in-place version of ∇⁴! which first
# allocates the output array ∇⁴alms before calling ∇⁴!."""
# function ∇⁴(alms::AbstractMatrix{Complex{NF}},  # spectral coefficients
#             R::Real=1) where NF                 # radius of the Earth
    
#     ∇⁴alms = copy(alms)                     # allocate output array
#     return R == 1 ? ∇⁴!(∇⁴alms,alms) : ∇⁴!(∇⁴alms,alms,R)
# end

# """
#     ∇⁻²!(   ∇⁻²alms::AbstractMatrix{Complex},       # Out: inverse Laplace of alms
#             alms::AbstractMatrix{Complex})          # In: spectral coefficients alms

# Inverse spherical Laplace operator ∇⁻² applied to the spectral coefficients `alms` on a sphere
# of radius `R=1`. ∇⁻²! is the in-place version which directly stores the output in the argument
# ∇⁻²alms. The integration constant for Legendre polynomial `l=m=0` is zero."""
# function ∇⁻²!(  ∇⁻²alms::AbstractMatrix{Complex{NF}},   # Output: inverse Laplacian of alms
#                 alms::AbstractMatrix{Complex{NF}}       # spectral coefficients
#                 ) where NF      

#     @boundscheck size(alms) == size(∇⁻²alms) || throw(BoundsError)
#     lmax,mmax = size(alms) .- 1     # degree l, order m of the Legendre polynomials

#     @inbounds for m in 1:mmax+1     # order m = 0:mmax but 1-based
#         for l in m:lmax+1           # degree l = 0:lmax but 1-based
#             # ∇⁻²alms = R²/(-l(l+1))*alms, but 1-based (l=>l-1) and R=1
#             ∇⁻²alms[l,m] = alms[l,m]/(l*(1-l))
#         end
#     end

#     # set the integration constant (l=m=0 polynomial) to zero 
#     ∇⁻²alms[1,1] = zero(∇⁻²alms[1,1]) 
#     return ∇⁻²alms
# end

# """
#     ∇⁻²alms = ∇⁻²(alms::AbstractMatrix{Complex},R::Real=1)

# Inverse spherical Laplace operator ∇⁻² applied to the spectral coefficients `alms` on a sphere
# of radius `R`. ∇⁻² is the non-in-place version of ∇⁻²! and therefore first allocates the output array
# ∇⁻²alms before calling ∇⁻²!. The integration constant for Legendre polynomial `l=m=0` is zero."""
# function ∇⁻²(   alms::AbstractMatrix{Complex{NF}},  # spectral coefficients
#                 R::Real=1) where NF                 # radius of the Earth
#     ∇⁻²alms = copy(alms)                            # preallocate output
#     return R == 1 ? ∇⁻²!(∇⁻²alms,alms) : ∇⁻²!(∇⁻²alms,alms,R)
# end

# function gradient_latitude!(coslat_u::AbstractMatrix{Complex{NF}},   # output: cos(lat)*zonal velocity u
#                             Ψ::AbstractMatrix{Complex{NF}},          # input: streamfunction Ψ
#                             ϵlms::AbstractMatrix{NF},                # recursion factors
#                             R::Real=1                               # radius of the sphere/Earth
#                             ) where {NF<:AbstractFloat}             # number format NF

#     _,mmax = size(Ψ)                # degree l, order m of spherical harmonics
#     lmax, mmax = mmax-1, mmax-1     # convert to 0-based l,m, but use mmax for lmax/lmax+1 flexibility
    
#     # u needs one more degree/meridional mode l for each m than Ψ due to the recursion
#     # Ψ can have size n+1 x n but then the last row is not used in the loop
#     size_compat = size(coslat_u) == size(Ψ) || (size(coslat_u) .- (1,0)) == size(Ψ)
#     @boundscheck size_compat || throw(BoundsError)
#     R⁻¹ = convert(Complex{NF},1/R)                              # 1/radius of the sphere

#     # for loops implement the recursion formula (0-based degree l, order m)
#     # (coslat*u)_lm = -1/R(  -(l-1)*ϵ(l,m)  *Ψ_(l-1,m)          # recursion term 1
#     #                       (l+2)*ϵ(l+1,m)*Ψ_(l+1,m))           # recursion term 2

#     # convert to 1-based l,m
#     @inbounds for m in 1:mmax         # exclude m=mmax+1 as for m=l=mmax+1 term1=term2=0

#         # 1. Diagonal (l=m), term 1 = 0 for the l=m modes
#         coslat_u[m,m] = -R⁻¹*(m+1)*ϵlms[m+1,m]*Ψ[m+1,m]         # recursion term 2 only
#         # coslat_u[m,m] = (1-m)*ϵlms[m+1,m]*Ψ[m+1,m]              # recursion term 2 only

#         # 2. Below diagonal modes (l>m, but l<=lmax)
#         for l in m+1:lmax
#             coslat_u[l,m] = -R⁻¹*(-(l-2)*ϵlms[l  ,m]*Ψ[l-1,m] + # term 1
#                                    (l+1)*ϵlms[l+1,m]*Ψ[l+1,m])  # term 2
#             # coslat_u[l,m] = (l-2)*ϵlms[l  ,m]*Ψ[l-1,m] -        # term 1
#                             # (l+1)*ϵlms[l+1,m]*Ψ[l+1,m]          # term 2
#         end

#         # 3. Last two rows
#         for l in lmax+1:lmax+2                                  # recursion term 1 only
#             coslat_u[l,m] = R⁻¹*(l-2)*ϵlms[l,m]*Ψ[l-1,m]
#             # coslat_u[l,m] = (l-2)*ϵlms[l,m]*Ψ[l-1,m]
#         end
#     end
#     # coslat_u[end,end] not needed as ϵ(l=m) = 0 in this case (l=lmax+1,m=mmax)

#     return coslat_u
# end

# function ∇⁻²!(  ∇⁻²alms::AbstractArray{Complex{NF},3},  # Output: inverse Laplacian of alms
#                 alms::AbstractArray{Complex{NF},3},     # spectral coefficients
#                 R::Real                                 # radius of the sphere/Earth
#                 ) where {NF<:AbstractFloat}

#     for k in 1:size(alms)[end]
#         ∇⁻²alms_layer = view(∇⁻²alms,:,:,k)
#         alms_layer = view(alms,:,:,k)
#         ∇⁻²!(∇⁻²alms_layer,alms_layer,R)
#     end
# end