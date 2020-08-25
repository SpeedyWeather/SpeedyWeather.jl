"""
Struct containing all the preallocated arrays for the calculation of
horizontal diffusion.
"""
struct HorizontalDiffusion{T<:AbstractFloat}
    tdrs::T
    dmp::Array{T,2}
    dmpd::Array{T,2}
    dmps::Array{T,2}
    tcorv::Array{T,1}
    qcorv::Array{T,1}
    tcorh::Array{T,2}
    qcorh::Array{T,2}
end

"""
Generator function for a HorizontalDiffusion struct.
"""
function HorizontalDiffusion(P::Params,G::GeoSpectral{T}) where {T<:AbstractFloat}

    @unpack nlev, σ_full = G.geometry
    @unpack trunc, mx, nx = G.spectral
    @unpack g, R, γ, hscale, hshum, rh_ref = P

    # Power of Laplacian in horizontal diffusion
    npowhd = 4.0

    thd    = 2.4       # Max damping time (in hours) for horizontal diffusion (del^6) of temperature
                       # and vorticity
    thdd   = 2.4       # Max damping time (in hours) for horizontal diffusion (del^6) of divergence
    thds   = 12.0      # Max damping time (in hours) for extra diffusion (del^2) in the stratosphere
    tdrs   = 24.0*30.0 # Damping time (in hours) for drag on zonal-mean wind in the stratosphere

    # Coefficients for horizontal diffusion
    # Spectral damping coefficients
    hdiff = 1.0/(3600.0thd)
    hdifd = 1.0/(3600.0thdd)
    hdifs = 1.0/(3600.0thds)
    rlap  = 1.0/(trunc*(trunc + 1))

    dmp = zeros(Float64, mx, nx)
    dmpd = zeros(Float64, mx, nx)
    dmps = zeros(Float64, mx,nx)
    for j in 1:nx
        for k in 1:mx
            N = k +j - 2
            elap = (N*(N + 1.0)*rlap)
            elapn = elap^npowhd
            dmp[k,j]  = hdiff*elapn
            dmpd[k,j] = hdifd*elapn
            dmps[k,j] = hdifs*elap
        end
    end

    # 5.2 Orographic correction terms for temperature and humidity
    #     (vertical component)
    γ_g = γ/(1000.0g)
    rgam = R*γ_g
    qexp = hscale/hshum

    tcorv = zeros(Float64, nlev)        # Temperature CORrection Vertical?
    qcorv = zeros(Float64, nlev)        # same for humidity?
    for k in 2:nlev
        tcorv[k] = σ_full[k]^rgam
        if k > 2
            qcorv[k] = σ_full[k]^qexp
        end
    end

    corh = zeros(T, nlon, nlat)
    for j in 1:nlat
        for i = 1:nlon
            corh[i,j] = γ_g*boundaries.ϕ₀ₛ[i,j]
        end
    end
    tcorh = grid_to_spec(corh)

    corh = refrh1
    qcorh = grid_to_spec(corh)

    HorizontalDiffusion(tdrs,dmp,dmpd,dmps,tcorv,qcorv,tcorh,qcorh)
end

function do_horizontal_diffusion_2d!(field, fdt, dmp, dmp1)
    fdt = (fdt - dmp.*field).*dmp1
end

# Add horizontal diffusion tendency of field to spectral tendency fdt at nlev
# levels using damping coefficients dmp and dmp1
function do_horizontal_diffusion_3d!(field, fdt, dmp, dmp1)
    for k in 1:nlev
        do_horizontal_diffusion_2d!(field[:,:,k], @view(fdt[:,:,k]), dmp, dmp1)
    end
end
