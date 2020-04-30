"""Struct holding the prognostic spectral variables."""
struct Prognostics{T<:AbstractFloat}
    ξ::Array{Complex{T},3}       # Vorticity [?]
    D::Array{Complex{T},3}       # Divergence [?]
    Tabs::Array{Complex{T},3}    # Absolute temperature [?]
    logp::Array{Complex{T},2}    # Log of surface pressure [?]
    humid::Array{Complex{T},3}   # Specific humidity [g/kg]
    ϕ::Array{Complex{T},3}       # Atmospheric geopotential [?]
    ϕ0::Array{Complex{T},2}      # Surface geopotential [?]
end

function initialize_from_rest(  ::Type{T},
                                ϕ0_grid::AbstractMatrix,
                                geometry::Geometry,
                                constants::Constants,
                                spectral_trans::SpectralTrans)

    @unpack nlon, nlat, nlev, mx, nx, σ_full = geometry
    @unpack g, R, γ, hscale, hshum, refrh1 = constants

    ξ     = zeros(Complex{T}, mx, nx, nlev)
    D     = zeros(Complex{T}, mx, nx, nlev)
    Tabs  = zeros(Complex{T}, mx, nx, nlev)
    logp  = zeros(Complex{T}, mx, nx)
    humid = zeros(Complex{T}, mx, nx, nlev)
    ϕ     = zeros(Complex{T}, mx, nx, nlev)

    # Compute spectral surface geopotential
    ϕ0    = Complex{T}.(grid_to_spec(ϕ0_grid,spectral_trans, geometry))

    surfg = zeros(T, nlon, nlat)

    γ_g = γ/(1000.0g)   # Lapse rate [K/m] scaled by gravity



    # 2. Start from reference atmosphere (at rest)

    # 2.2 Set reference temperature :
    #     tropos:  T = 288 degK at z = 0, constant lapse rate
    #     stratos: T = 216 degK, lapse rate = 0
    Tabs_ref  = 288.0   #TODO pass on from param
    Tabs_top  = 216.0

    # Surface and stratospheric air temperature
    surfs = -γ_g*ϕₛ

    Tₐ[1,1,1,1] = √(2.0)*Complex{T}(1.0)*Tₐ_top
    Tₐ[1,1,2,1] = √(2.0)*Complex{T}(1.0)*Tₐ_top
    surfs[1,1] = √(2.0)*Complex{T}(1.0)*Tₐ_ref - γ_g*ϕₛ[1,1]

    # Temperature at tropospheric levels
    for k in 3:nlev
        Tₐ[:,:,k,1] = surfs*σ_full[k]^(R*γ_g)
    end

    # 2.3 Set log(ps) consistent with temperature profile
    #     p_ref = 1013 hPa at z = 0
    for j in 1:nlat
        for i in 1:nlon
            surfg[i,j] = log(1.013) + log(1.0 - γ_g*ϕ₀ₛ[i,j]/Tₐ_ref)/(R*γ_g)
        end
    end

    pₛ[:,:,1] = truncate(spectral_trans.trfilt, grid_to_spec(geometry, spectral_trans, surfg))

    # 2.4 Set tropospheric specific humidity in g/kg
    #     Qref = RHref * Qsat(288K, 1013hPa)
    esref = 17.0
    qref = refrh1*0.622*esref
    qexp = hscale/hshum

    # Specific humidity at the surface
    for j in 1:nlat
        for i in 1:nlon
            surfg[i,j] = qref*exp(qexp*surfg[i,j])
        end
    end

    surfs = truncate(spectral_trans.trfilt, grid_to_spec(geometry, spectral_trans, surfg))

    # Specific humidity at tropospheric levels
    for k in 3:nlev
        tr[:,:,k,1,1] = surfs*σ_full[k]^qexp
    end

    # Print diagnostics from initial conditions
    check_diagnostics(geometry, spectral_trans, ξ[:,:,:,1], D[:,:,:,1], Tₐ[:,:,:,1], 0)

    Prognostics(ξ, D, Tabs, logp, humid, ϕ, ϕ0)
end
