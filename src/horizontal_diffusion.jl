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
function HorizontalDiffusion(   P::Params,
                                G::GeoSpectral{NF},
                                B::Boundaries{NF}) where {NF<:AbstractFloat}

    @unpack nlev, σ_full = G.geometry
    @unpack trunc, mx, nx = G.spectral
    @unpack g, R, γ, hscale, hshum, rh_ref = P
    # @unpack npowhd, thd, thdd, thds, tdrs = P
    @unpack ϕ0trunc = B

    # TODO load from paramters struct instead
    npowhd = 4.0        # Power of Laplacian in horizontal diffusion
    thd    = 2.4        # Damping time [hrs] for diffusion (del^6) of temperature and vorticity
    thdd   = 2.4        # Damping time [hrs] for diffusion (del^6) of divergence
    thds   = 12.0       # Damping time [hrs] for extra diffusion (del^2) in the stratosphere
    tdrs   = 24.0*30.0  # Damping time [hrs] for drag on zonal-mean wind in the stratosphere

    # Damping frequencies [1/s]
    hdiff = 1.0/(3600.0*thd)            # Spectral damping for temperature and vorticity
    hdifd = 1.0/(3600.0*thdd)           # Spectral damping coefficient for divergence
    hdifs = 1.0/(3600.0*thds)           # Spectral damping coefficient for stratosphere
    rlap  = 1.0/(trunc*(trunc + 1))     # ?

    # preallocate (high precision, conversion to NF later)
    dmp = zeros(mx, nx)                 # Damping coefficient for temperature and vorticity (explicit)
    dmpd = zeros(mx, nx)                # Damping coefficient for divergence (explicit)
    dmps = zeros(mx,nx)                 # Damping coefficient for extra diffusion in the stratosphere (explicit)

    for j in 1:nx
        for k in 1:mx
            N = k+j-2
            elap = (N*(N + 1.0)*rlap)
            elapn = elap^npowhd
            dmp[k,j]  = hdiff*elapn
            dmpd[k,j] = hdifd*elapn
            dmps[k,j] = hdifs*elap
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

    tcorh = convert_to_spectral(corh,G) # correction in spectral space

    # Orographic correction terms for humidity (horizontal component)
    corh .= rh_ref                  # relative humidity reference value
    qcorh = convert_to_spectral(corh)

    return HorizontalDiffusion{NF}( dmp,dmpd,dmps,
                                    dmp1,dmpd1,dmps1,
                                    tcorv,qcorv,tcorh,qcorh)
end

"""Apply horizontal diffusion to 2D field in spectral space."""
function do_horizontal_diffusion!(  A::AbstractArray{Complex{NF},2}         # spectral horizontal field
                                    tendency::AbstractArray{Complex{NF},2}, # its tendency
                                    dmp::AbstractArray{NF,2},               # damping coefficients (explicit)
                                    dmp1::AbsrtactArray{NF,2}               # damping coefficients (implicit)
                                    ) where {NF<:AbstractFloat}

    @boundscheck size(A) == size(tendency) || throw(BoundsError())
    @boundscheck size(A) == size(dmp) || throw(BoundsError())
    @boundscheck size(A) == size(dmp1) || throw(BoundsError())
    
    tendency = (tendency - dmp*field)*dmp1
end

"""Apply horizontal diffusion to 3D field layer by layer in spectral space."""
function do_horizontal_diffusion!(  A::AbstractArray{Complex{NF},3}         # spectral horizontal field
                                    tendency::AbstractArray{Complex{NF},3}, # its tendency
                                    dmp::AbstractArray{NF,2},               # damping coefficients (explicit)
                                    dmp1::AbsrtactArray{NF,2}               # damping coefficients (implicit)
                                    ) where {NF<:AbstractFloat}
    mx,nx,nlev = size(A)
    @boundscheck (mx,nx,nlev) == size(tendency) || throw(BoundsError())
    
    for k in 1:nlev
        A_layer = view(A,:,:,k)
        tendency_layer = view(tendency,:,:,k)
        do_horizontal_diffusion!(A_layer, tendency_layer, dmp, dmp1)
    end
end