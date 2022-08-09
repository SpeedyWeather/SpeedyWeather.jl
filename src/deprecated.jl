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

# """gradient_latitude! but precalculate the recursion factors `ϵlms` in case they are not provided."""
# function gradient_latitude!(coslat_u::AbstractMatrix{Complex{NF}},  # output: cos(lat)*u
#                             Ψ::AbstractMatrix{Complex{NF}},         # input: streamfunction Ψ
#                             R::Real=1                               # radius of the sphere/Earth
#                             ) where {NF<:AbstractFloat}             # number format NF
#     _,mmax = size(Ψ) .- 1                                           # degree l, order m of spherical harmonics   
#     ϵlms = get_recursion_factors(NF,mmax,mmax)                      # precalculate recursion factors
#     return gradient_latitude!(coslat_u,Ψ,ϵlms,R)                    # call in-place function
# end

# function gradient_latitude( Ψ::AbstractMatrix{Complex{NF}}, # input: streamfunction Ψ
#                             R::Real=1                       # radius of the sphere/Earth
#                             ) where {NF<:AbstractFloat}     # number format NF
#     _,mmax = size(Ψ) .- 1                                   # max degree l, order m of spherical harmonics
#     coslat_u = zeros(Complex{NF},mmax+2,mmax+1)             # preallocate output, one more l for recursion
#     return gradient_latitude!(coslat_u,Ψ,R)                 # call in-place version
# end

# create a model hierarchy: BarotropicModel < ShallowWaterModel < PrimitiveEquationModel
# for instances
# Base.:(<)(M1::BarotropicModel,M2::BarotropicModel) = false
# Base.:(<)(M1::ShallowWaterModel,M2::ShallowWaterModel) = false
# Base.:(<)(M1::PrimitiveEquationModel,M2::PrimitiveEquationModel) = false

# Base.:(<)(M1::BarotropicModel,M2::ShallowWaterModel) = true
# Base.:(<)(M1::BarotropicModel,M2::PrimitiveEquationModel) = true
# Base.:(<)(M1::ShallowWaterModel,M2::PrimitiveEquationModel) = true

# Base.:(<)(M1::ShallowWaterModel,M2::BarotropicModel) = false
# Base.:(<)(M1::PrimitiveEquationModel,M2::BarotropicModel) = false
# Base.:(<)(M1::PrimitiveEquationModel,M2::ShallowWaterModel) = false

# # and also for types
# Base.:(<)(::Type{ShallowWaterModel{NF1}},::Type{BarotropicModel{NF2}}) where {NF1,NF2} = false
# Base.:(<)(::Type{PrimitiveEquationModel{NF1}},::Type{BarotropicModel{NF2}}) where {NF1,NF2} = false
# Base.:(<)(::Type{PrimitiveEquationModel{NF1}},::Type{ShallowWaterModel{NF2}}) where {NF1,NF2} = false

# Base.:(<)(::Type{BarotropicModel{NF1}},::Type{ShallowWaterModel{NF2}}) where {NF1,NF2} = true
# Base.:(<)(::Type{BarotropicModel{NF1}},::Type{PrimitiveEquationModel{NF2}}) where {NF1,NF2} = true
# Base.:(<)(::Type{ShallowWaterModel{NF1}},::Type{PrimitiveEquationModel{NF2}}) where {NF1,NF2} = true

# Base.:(<)(::Type{BarotropicModel{NF1}},::Type{BarotropicModel{NF2}}) where {NF1,NF2} = false
# Base.:(<)(::Type{ShallowWaterModel{NF1}},::Type{ShallowWaterModel{NF2}}) where {NF1,NF2} = false
# Base.:(<)(::Type{PrimitiveEquationModel{NF1}},::Type{PrimitiveEquationModel{NF2}}) where {NF1,NF2} = false

# Base.:(<)(M::ModelSetup{NF},::Type{BarotropicModel}) where NF = typeof(M) < BarotropicModel{NF}
# Base.:(<)(M::ModelSetup{NF},::Type{ShallowWaterModel}) where NF = typeof(M) < ShallowWaterModel{NF}
# Base.:(<)(M::ModelSetup{NF},::Type{PrimitiveEquationModel}) where NF = typeof(M) < PrimitiveEquationModel{NF}

# """Compute curl and divergence of vectors u,v in one go. Currently not used."""
# function curl_div!( curl::AbstractMatrix{Complex{NF}},
#                     div::AbstractMatrix{Complex{NF}},
#                     u::AbstractMatrix{Complex{NF}},
#                     v::AbstractMatrix{Complex{NF}},
#                     S::SpectralTransform{NF}
#                     ) where {NF<:AbstractFloat}

#     lmax,mmax = size(div) .- (1,1)                  # 0-based lmax,mmax
#     @boundscheck size(curl) == size(div) || throw(BoundsError)
#     @boundscheck size(u) == (lmax+2,mmax+1) || throw(BoundsError)
#     @boundscheck size(v) == (lmax+2,mmax+1) || throw(BoundsError)

#     @unpack grad_y_vordiv1,grad_y_vordiv2 = S

#     div[1,1]  = zero(Complex{NF})                   # l=m=0 harmonic is zero
#     curl[1,1] = zero(Complex{NF})                   # l=m=0 harmonic is zero
#     @inbounds for m in 1:mmax+1                     # 1-based l,m
#         for l in max(2,m):lmax+1                    # skip l=m=0 harmonic (mean) to avoid access to v[0,1]
#             ∂u∂λ  = ((m-1)*im)*u[l,m]
#             ∂v∂λ  = ((m-1)*im)*v[l,m]
#             ∂v∂θ1 = grad_y_vordiv1[l,m]*v[l-1,m]
#             ∂u∂θ1 = grad_y_vordiv1[l,m]*u[l-1,m]
#             ∂v∂θ2 = grad_y_vordiv2[l,m]*v[l+1,m]
#             ∂u∂θ2 = grad_y_vordiv2[l,m]*u[l+1,m]
#             div[l,m]  = ∂u∂λ - ∂v∂θ1 + ∂v∂θ2
#             curl[l,m] = ∂v∂λ + ∂u∂θ1 - ∂u∂θ2
#         end
#     end

#     return nothing
# end
# """
#     gradient_latitude!( coslat_u::AbstractArray{Complex{NF}},   # output: cos(lat)*zonal velocity u
#                         Ψ::AbstractArray{Complex{NF}},          # input: streamfunction Ψ
#                         ϵlms::AbstractArray{NF},                # recursion factors
#                         R::Real=1                               # radius of the sphere/Earth
#                         ) where {NF<:AbstractFloat}             # number format NF

# Meridional gradient in spectral space of spherical harmonic coefficients `Ψ` on a sphere with
# radius R. Returns `coslat_u`, i.e. the gradient ∂Ψ/∂lat with an additional cosine of latitude scaling.
# This function uses the recursion relation (0-based degree l, order m)

#     (coslat u)_lm = -1/R*(-(l-1)*ϵ_lm*Ψ_l-1,m + (l+2)*ϵ_l+1,m*Ψ_l+1,m ).
    
# As u = -1/R*∂Ψ/∂lat, this function can be generally used to compute the gradient in latitude."""
# function gradient_latitude!(coslat_u::AbstractMatrix{Complex{NF}},  # output: cos(lat)*zonal velocity u
#                             Ψ::AbstractMatrix{Complex{NF}},         # input: streamfunction Ψ
#                             S::SpectralTransform{NF};               # use precomputed recursion factors
#                             flipsign::Bool=false,                   # flip the sign to obtain u from Ψ
#                             add::Bool=false                         # coslat_u += (add) or = (overwrite)
#                             ) where {NF<:AbstractFloat}             # number format NF

#     # u needs one more degree/meridional mode l for each m than Ψ due to the recursion
#     # Ψ can have size n+1 x n but then the last row is not used in the loop
#     lmax_out, mmax_out = size(coslat_u)
#     lmax_in,  mmax_in  = size(Ψ)

#     @boundscheck mmax_out == mmax_in || throw(BoundsError)
#     mmax = mmax_out - 1                 # 0-based max order m of harmonics

#     @boundscheck abs(lmax_out-lmax_in) <= 1 || throw(BoundsError)
#     lmax = lmax_in - 1                  # 0-based max degree l of harmonics
#     output_larger = lmax_out - lmax_in

#     if flipsign                     # used to get u from streamfunction Ψ (u ~ -∂Ψ/∂lat)
#         grad_y1 = S.minus_grad_y1
#         grad_y2 = S.minus_grad_y2
#     else                            # no sign flip otherwise
#         grad_y1 = S.grad_y1
#         grad_y2 = S.grad_y2
#     end

#     g = grad_y2[1,1]*Ψ[2,1]
#     coslat_u[1,1] = add ? coslat_u[1,1] + g : g         # l=m=0 mode only with term 2

#     for m in 1:mmax+1
#         for l in max(2,m):lmax
#             g = grad_y1[l,m]*Ψ[l-1,m] + grad_y2[l,m]*Ψ[l+1,m]
#             coslat_u[l,m] = add ? coslat_u[l,m] + g : g
#         end
#         for l in lmax+1:lmax+1+output_larger            # execute 2x for out > in, 1x for out == in and
#             g = grad_y1[l,m]*Ψ[l-1,m]                   # 0x for out < in
#             coslat_u[l,m] = add ? coslat_u[l,m] + g : g
#         end
#     end

#     return coslat_u
# end

# function gradient_latitude( Ψ::AbstractMatrix{Complex{NF}}, # input: streamfunction Ψ
#                             S::SpectralTransform{NF};       # precomputed gradient arrays
#                             one_more_l::Bool=true,          # allocate output with one more degree l?
#                             flipsign::Bool=false            # flip the sign to obtain u from Ψ?
#                             ) where {NF<:AbstractFloat}     # number format NF
#     _,mmax = size(Ψ) .- 1                                   # max degree l, order m of spherical harmonics
#     coslat_u = zeros(Complex{NF},mmax+one_more_l+1,mmax+1)  # preallocate output, one more l for recursion
#     return gradient_latitude!(coslat_u,Ψ,S;flipsign)        # call in-place version
# end

# """
#     coslat_v = gradient_longitude!( coslat_v::AbstractMatrix{Complex{NF}},
#                                     Ψ::AbstractMatrix{Complex{NF}};
#                                     radius::Real=1
#                                     ) where {NF<:AbstractFloat}

# Zonal gradient in spectral space of spherical harmonic coefficients `Ψ` on a sphere with radius `radius`.
# While the zonal gradient has a 1/cos(lat) scaling in spherical coordinates, this functions omits the scaling
# such that the returned array is scaled with coslat.
# """
# function gradient_longitude!(   coslat_v::AbstractMatrix{Complex{NF}},  # output: cos(latitude)*meridional velocity
#                                 Ψ::AbstractMatrix{Complex{NF}},         # input: spectral coefficients of stream function
#                                 radius::Real=1;                         # radius of the sphere/Earth
#                                 add::Bool=false,                        # coslat_u += (add) or = (overwrite)
#                                 flipsign::Bool=false                    # flip sign of output?
#                                 ) where {NF<:AbstractFloat}             # number format NF

#     lmax_out, mmax_out = size(coslat_v)
#     lmax_in,  mmax_in  = size(Ψ)

#     @boundscheck mmax_out == mmax_in || throw(BoundsError)
#     mmax = mmax_out - 1                 # 0-based max order m of harmonics

#     @boundscheck abs(lmax_out-lmax_in) <= 1 || throw(BoundsError)
#     lmax = min(lmax_out,lmax_in) - 1    # 0-based max degree l of harmonics
#     output_larger = lmax_out > lmax_in

#     iradius⁻¹ = convert(Complex{NF},(-1)^flipsign*im/radius)            # = ±imaginary/radius converted to NF

#     @inbounds for m in 1:mmax+1                     # loop over all coefficients, order m
#         for l in m:lmax+1                           # degree l
#             g = (m-1)*iradius⁻¹*Ψ[l,m]              # gradient in lon = *i*m/radius but 1-based order
#             coslat_v[l,m] = add ? coslat_v[l,m] + g : g 
#         end
#     end

#     # if grad in lon does not project onto the last degree l, set explicitly to zero in that case
#     if output_larger & ~add
#         for m in 1:mmax+1        
#             coslat_v[end,m] = zero(Complex{NF}) 
#         end
#     end

#     return coslat_v
# end

# """Gradient in longitude in spectral space. Input: coefficients `alms` of the spherical harmonics."""
# function gradient_longitude(Ψ::AbstractMatrix{NF},  # input array: spectral coefficients
#                             R::Real=1;              # radius of the sphere/Earth
#                             one_more_l::Bool=true   # allocate output with one more degree l
#                             ) where NF              # number format NF
#     _,mmax = size(Ψ)
#     coslat_v = zeros(NF,mmax+one_more_l,mmax)       # preallocate output array (gradient in longitude)
#     return gradient_longitude!(coslat_v,Ψ,R)        # call in-place version
# end