"""
Struct containing all the preallocated arrays for the calculation of horizontal diffusion.
"""
struct HorizontalDiffusion{NF<:AbstractFloat}   # Number format NF
    dmp::Array{NF,2}                # Damping coefficient for temperature and vorticity (explicit)
    dmpd::Array{NF,2}               # Damping coefficient for divergence (explicit)
    dmps::Array{NF,2}               # Damping coefficient for extra diffusion in the stratosphere (explicit)
    
    dmp1::Array{NF,2}               # Damping coefficient for temperature and vorticity (implicit)
    dmp1d::Array{NF,2}              # Damping coefficient for divergence (implicit)
    dmp1s::Array{NF,2}              # Damping coefficient for extra diffusion in the stratosphere (implicit)

    tcorv::Array{NF,1}              # Vertical component of orographic correction for temperature
    qcorv::Array{NF,1}              # Vertical component of orographic correction for humidity
    
    tcorh::Array{Complex{NF},2}     # Horizontal component of orographic correction for temperature
    qcorh::Array{Complex{NF},2}     # Horizontal component of orographic correction for humidity
end

"""
Generator function for a HorizontalDiffusion struct.
"""
function HorizontalDiffusion(   P::Parameters,      # Parameter struct
                                G::GeoSpectral,     # Geometry and spectral struct
                                B::Boundaries)      # Boundaries struct

    @unpack nlev, σ_full = G.geometry
    @unpack trunc, mx, nx = G.spectral
    @unpack g, R, γ, hscale, hshum, rh_ref = P
    @unpack npowhd, thd, thdd, thds, tdrs = P
    @unpack ϕ0trunc = B

    # Damping frequencies [1/s]
    hdiff = 1.0/(3600.0*thd)            # Spectral damping for temperature and vorticity
    hdifd = 1.0/(3600.0*thdd)           # Spectral damping coefficient for divergence
    hdifs = 1.0/(3600.0*thds)           # Spectral damping coefficient for stratosphere
    rlap  = 1.0/(trunc*(trunc + 1))     # ?

    # preallocate (high precision, conversion to NF later)
    dmp = zeros(mx,nx)                  # Damping coefficient for temperature and vorticity (explicit)
    dmpd = zeros(mx,nx)                 # Damping coefficient for divergence (explicit)
    dmps = zeros(mx,nx)                 # Damping coefficient for extra diffusion in the stratosphere (explicit)

    dmp1 = zeros(mx,nx)                 # Damping coefficient for temperature and vorticity (implicit)
    dmp1d = zeros(mx,nx)                # Damping coefficient for divergence (implicit)
    dmp1s = zeros(mx,nx)                # Damping coefficient for extra diffusion in the stratosphere (implicit)

    for j in 1:nx
        for i in 1:mx
            N = i+j-2
            elap = (N*(N + 1.0)*rlap)
            elapn = elap^npowhd
            dmp[i,j]  = hdiff*elapn     
            dmpd[i,j] = hdifd*elapn
            dmps[i,j] = hdifs*elap

            dmp1[i,j]  = 1/(1+dmp*Δt)
            dmp1d[i,j] = 1/(1+dmpd*Δt)
            dmp1s[i,j] = 1/(1+dmps*Δt)
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
        tcorv[k] = σ_full[k]^rgam
        if k > 2
            qcorv[k] = σ_full[k]^qexp
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
    return HorizontalDiffusion{P.NF}(   dmp,dmpd,dmps,
                                        dmp1,dmpd1,dmps1,
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

# """Orographic temperature correction for absolute temperature to be applied before the horizontal diffusion."""
# function temperature_correction!(   Tabs_corrected::AbstractArray{Complex{NF},3},   # Corrected abs temperature T
#                                     Tabs::AbstractArray{Complex{NF},4},             # Absolute temperature
#                                     l::Int,                                         # leapfrog index
#                                     HD::HorizontalDiffusion{NF}                     # struct for correction arrays
#                                     ) where NF
    
#     @unpack tcorh, tcorv = HD   # load horizontal (tcorh) and vertical (tcorv) correction arrays

#     mx,nx,nlev,nleapfrog = size(Tabs)
#     @boundscheck (mx,nx,nlev) == size(Tabs_corrected) || throw(BoundsError())
#     @boundscheck (mx,nx) == size(tcorh) || throw(BoundsError())
#     @boundscheck (nlev,) == size(tcorv) || throw(BoundsError())
#     @boundscheck nleapfrog == 2 || throw(BoundsError())
#     @boundscheck l in [1,2] || throw(BoundsError())

#     @inbounds for k in 1:nlev
#         for j in 1:nx
#             for i in 1:mx
#                 Tabs_corrected[i,j,k] = Tabs[i,j,k,l] + tcorh[i,j]*tcorv[k]
#             end
#         end
#     end
# end

# """Orographic temperature correction for absolute temperature to be applied before the horizontal diffusion."""
# function humidity_correction!(  humid_corrected::AbstractArray{Complex{NF},3},  # Corrected abs temperature T
#                                 humid::AbstractArray{Complex{NF},4},            # Absolute temperature
#                                 l::Int,                                         # leapfrog index
#                                 HD::HorizontalDiffusion{NF}                     # struct for correction arrays
#                                 ) where NF
    
#     @unpack qcorh, qcorv = HD       # load horizontal (qcorh) and vertical (qcorv) correction arrays

#     mx,nx,nlev,nleapfrog = size(humid)
#     @boundscheck (mx,nx,nlev) == size(humid_corrected) || throw(BoundsError())
#     @boundscheck (mx,nx) == size(qcorh) || throw(BoundsError())
#     @boundscheck (nlev,) == size(qcorv) || throw(BoundsError())
#     @boundscheck nleapfrog == 2 || throw(BoundsError())
#     @boundscheck l in [1,2] || throw(BoundsError())

#     @inbounds for k in 1:nlev
#         for j in 1:nx
#             for i in 1:mx
#                 humid_corrected[i,j,k] = humid[i,j,k,l] + qcorh[i,j]*qcorv[k]
#             end
#         end
#     end
# end

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