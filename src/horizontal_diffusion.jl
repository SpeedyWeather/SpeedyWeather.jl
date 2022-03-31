"""
Struct containing all the preallocated arrays for the calculation of horizontal diffusion.
"""
struct HorizontalDiffusion{NF<:AbstractFloat}   # Number format NF
    # Explicit part of the diffusion, precalculated damping coefficients for each spectral mode
    damping::Array{NF,2}            # for temperature and vorticity (explicit)
    damping_div::Array{NF,2}        # for divergence (explicit)
    damping_strat::Array{NF,2}      # for extra diffusion in the stratosphere (explicit)
    
    # Implicit part of the diffusion, precalculated damping coefficients for each spectral mode
    damping_impl::Array{NF,2}       # for temperature and vorticity (implicit)
    damping_div_impl::Array{NF,2}   # for divergence (implicit)
    damping_strat_impl::Array{NF,2} # for extra diffusion in the stratosphere (implicit)
    
    # Vertical component of orographic correction
    tcorv::Array{NF,1}              # for temperature
    qcorv::Array{NF,1}              # for humidity
    
    # Horizontal component of orographic correction
    tcorh::Array{Complex{NF},2}     # for temperature
    qcorh::Array{Complex{NF},2}     # for humidity
end

"""
Generator function for a HorizontalDiffusion struct.
"""
function HorizontalDiffusion(   P::Parameters{NF},      # Parameter struct
                                G::GeoSpectral{NF},     # Geometry and spectral struct
                                B::Boundaries{NF}       # Boundaries struct
                                ) where NF              # number format NF

    # for diffusion
    @unpack lmax,mmax = G.spectral
    @unpack diffusion_power, diffusion_time, diffusion_time_div = P
    @unpack diffusion_time_strat, damping_time_strat = P
    
    # for orographic correction
    @unpack nlev, σ_levels_full = G.geometry
    @unpack g, R, γ, hscale, hshum, rh_ref = P
    @unpack ϕ0trunc = B

    # Diffusion is applied by multiplication of the (absolute) eigenvalues of the Laplacian l*(l+1)
    # normalise by the largest eigenvalue lmax*(lmax+1) such that the highest wavenumber lmax
    # is dampened to 0 at the time scale diffusion_time
    # raise to a power of the Laplacian for hyperdiffusion (=less damping for smaller wavenumbers)
    largest_eigenvalue = lmax*(lmax+1)

    # PREALLOCATE
    # conversion to number format NF later, one more degree l for meridional gradient recursion
    # Damping coefficients for explicit part of the diffusion (=ν∇²ⁿ)
    damping = zeros(lmax+2,mmax+1)              # for temperature and vorticity (explicit)
    damping_div = zeros(lmax+2,mmax+1)          # for divergence (explicit)
    damping_strat = zeros(lmax+2,mmax+1)        # for extra diffusion in the stratosphere (explicit)

    # Damping coefficients for implicit part of the diffusion (= 1/(1+2Δtν∇²ⁿ))
    damping_impl = zeros(lmax+2,mmax+1)         # for temperature and vorticity (implicit)
    damping_div_impl = zeros(lmax+2,mmax+1)     # for divergence (implicit)
    damping_strat_impl = zeros(lmax+2,mmax+1)   # for extra diffusion in the stratosphere (implicit)

    # PRECALCULATE the damping coefficients for every spectral mode
    for m in 1:mmax+1                           # fill only the lower triangle
        for l in m:lmax+2
            # eigenvalue is l*(l+1), but 1-based here l→l-1
            norm_eigenvalue = l*(l-1)/largest_eigenvalue        # normal diffusion ∇²
            norm_eigenvalueⁿ = norm_eigenvalue^diffusion_power  # hyper diffusion ∇²ⁿ

            # Explicit part (=ν∇²ⁿ)
            # convert diffusion time scales to damping frequencies [1/s] times norm. eigenvalue
            damping[l,m] = norm_eigenvalueⁿ/(3600*diffusion_time)               # temperature/vorticity
            damping_div[l,m] = norm_eigenvalueⁿ/(3600*diffusion_time_div)       # divergence
            damping_strat[l,m] = norm_eigenvalue/(3600*diffusion_time_strat)    # stratosphere (no hyperdiff)

            # and implicit part of the diffusion (= 1/(1+2Δtν∇²ⁿ))
            damping_impl[l,m] = 1/(1+2Δt*damping[l,m])              # for temperature/vorticity
            damping_div_impl[l,m] = 1/(1+2Δt*damping_div[l,m])      # for divergence
            damping_strat_impl[l,m] = 1/(1+2Δt*damping_strat[l,m])  # for stratosphere (no hyperdiffusion)
        end
    end

    # Orographic correction terms for temperature and humidity (vertical component)
    γ_g = γ/(1000.0g)
    rgam = R*γ_g
    qexp = hscale/hshum

    # preallocate (high precision, conversion to NF later)
    tcorv = zeros(nlev)    # Vertical component of orographic correction for temperature
    qcorv = zeros(nlev)    # Vertical component of orographic correction for humidity

    for k in 2:nlev
        tcorv[k] = σ_levels_full[k]^rgam
        if k > 2
            qcorv[k] = σ_levels_full[k]^qexp
        end
    end

    # Orographic correction term for temperature (horizontal component)
    corh = zeros(nlon, nlat)            # in grid-point space

    for j in 1:nlat
        for i = 1:nlon
            corh[i,j] = γ_g*ϕ0trunc[i,j]
        end
    end

    tcorh = spectral(corh,G) # correction in spectral space

    # Orographic correction terms for humidity (horizontal component)
    corh .= rh_ref                  # relative humidity reference value
    qcorh = spectral(corh,G)

    # convert to number format NF here
    return HorizontalDiffusion{NF}( damping,damping_div,damping_strat,
                                    damping_impl,damping_div_impl,damping_strat_impl,
                                    tcorv,qcorv,tcorh,qcorh)
end

"""Apply horizontal diffusion to 2D field in spectral space."""
function horizontal_diffusion!( A::AbstractArray{Complex{NF},2},        # spectral horizontal field
                                tendency::AbstractArray{Complex{NF},2}, # its tendency
                                dmp::AbstractArray{NF,2},               # damping coefficients (explicit)
                                dmp1::AbstractArray{NF,2}               # damping coefficients (implicit)
                                ) where {NF<:AbstractFloat}

    @boundscheck size(A) == size(tendency) || throw(BoundsError())
    @boundscheck size(A) == size(dmp) || throw(BoundsError())
    @boundscheck size(A) == size(dmp1) || throw(BoundsError())
    
    @inbounds for i in eachindex(A)
        tendency[i] = (tendency[i] - dmp[i]*A[i])*dmp1[i]
    end
end



"""Apply horizontal diffusion to 3D field layer by layer in spectral space."""
function horizontal_diffusion!( A::AbstractArray{Complex{NF},3},        # spectral horizontal field
                                tendency::AbstractArray{Complex{NF},3}, # its tendency
                                dmp::AbstractArray{NF,2},               # damping coefficients (explicit)
                                dmp1::AbstractArray{NF,2}               # damping coefficients (implicit)
                                ) where {NF<:AbstractFloat}             # number format NF
    _,_,nlev = size(A)
    @boundscheck size(A) == size(tendency) || throw(BoundsError())
    
    for k in 1:nlev
        A_layer = view(A,:,:,k)
        tendency_layer = view(tendency,:,:,k)
        horizontal_diffusion!(A_layer, tendency_layer, dmp, dmp1)
    end
end

"""Zonal drag in the stratosphere (uppermost layer) for vorticity and divergence."""
function stratospheric_zonal_drag!( A::AbstractArray{Complex{NF},4},        # spectral vorticity or divergence
                                    tendency::AbstractArray{Complex{NF},3}, # its tendency
                                    sdrag::NF                               # drag coefficient [1/s]
                                    ) where {NF<:AbstractFloat}             # number format NF
    
    mx,nx,nlev,nleapfrog = size(A)
    @boundscheck (mx,nx,nlev) == size(tendency) || throw(BoundsError())
    @boundscheck nleapfrog == 2 || throw(BoundsError())

    @inbounds for j in 1:nx
        # size(A) = mx x nx x nlev, nlev = 1 is uppermost model level
        # apply drag only to largest zonal wavenumber (mx = 1)
        tendency[1,j,1] = tendency[1,j,1] - sdrag*A[1,j,1,1]
    end
end

"""Orographic temperature correction for absolute temperature to be applied before the horizontal diffusion."""
function orographic_correction!(A_corrected::AbstractArray{Complex{NF},3},  # Corrected variable
                                A::AbstractArray{Complex{NF},4},            # Variable (temperature or humidity)
                                l::Int,                                     # leapfrog index
                                hori_correction::AbstractArray,             # horizontal correction array
                                vert_correction::AbstractArray,             # vertical correction array
                                ) where NF
    
    mx,nx,nlev,nleapfrog = size(A)
    @boundscheck (mx,nx,nlev) == size(A_corrected) || throw(BoundsError())
    @boundscheck (mx,nx) == size(hori_correction) || throw(BoundsError())
    @boundscheck (nlev,) == size(vert_correction) || throw(BoundsError())
    @boundscheck nleapfrog == 2 || throw(BoundsError())
    @boundscheck l in [1,2] || throw(BoundsError())

    @inbounds for k in 1:nlev
        for j in 1:nx
            for i in 1:mx
                A_corrected[i,j,k] = A[i,j,k,l] + hori_correction[i,j]*vert_correction[k]
            end
        end
    end
end