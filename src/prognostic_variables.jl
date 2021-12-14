"""Struct holding the prognostic spectral variables."""
struct PrognosticVariables{NF<:AbstractFloat}
    vor     ::Array{Complex{NF},3}      # Vorticity of horizontal wind field
    div     ::Array{Complex{NF},3}      # Divergence of horizontal wind field
    Tabs    ::Array{Complex{NF},3}      # Absolute temperature [K]
    logp    ::Array{Complex{NF},2}      # Log of surface pressure [log(Pa)]
    humid   ::Array{Complex{NF},3}      # Specific humidity [g/kg]
    # ϕ       ::Array{Complex{NF},3}      # Atmospheric geopotential [m²/s²]
    # ϕ0      ::Array{Complex{NF},2}      # Surface geopotential [m²/s²]
end

"""Initialize prognostic variables from rest or restart from file."""
function initial_conditions(    P::Params,                      # Parameter struct
                                B::Boundaries{NF},              # Boundaries struct
                                G::GeoSpectral{NF}              # GeoSpectral struct
                                ) where {NF<:AbstractFloat}     # number format NF

    @unpack initial_conditions = P

    # if initial_conditions == :rest
    Progs = initialize_from_rest(P,B,G)
    # else initial_conditions == :rest
    #TODO allow for restart from file

    return Progs
end


function initialize_from_rest(  P::Params,
                                B::Boundaries{T},
                                G::GeoSpectral{T}
                                ) where {NF<:AbstractFloat}

    @unpack nlev = G.geometry
    @unpack mx, nx = G.spectral
    @unpack ϕ0 = B

    # conversion to type T later when creating a Prognostics struct
    ξ     = zeros(Complex{Float64}, mx, nx, nlev)   # Vorticity
    D     = zeros(Complex{Float64}, mx, nx, nlev)   # Divergence
    Tabs  = zeros(Complex{Float64}, mx, nx, nlev)   # Absolute Temperature
    logp  = zeros(Complex{Float64}, mx, nx)         # logarithm of surface pressure
    humid = zeros(Complex{Float64}, mx, nx, nlev)   # specific humidity
    ϕ     = zeros(Complex{Float64}, mx, nx, nlev)   # geopotential

    # Compute spectral surface geopotential
    ϕ0spectral    = spectral(Float64.(ϕ0),G)

    # TEMPERATURE
    # Tabs_ref: Reference absolute T [K] at surface z = 0, constant lapse rate
    # Tabs_top: Reference absolute T in the stratosphere [K], lapse rate = 0
    @unpack Tabs_ref, Tabs_top, γ, g, R = P
    @unpack σ_full = G.geometry

    γ_g = γ/g/1000                  # Lapse rate [K/m] scaled by gravity
    T0 = -γ_g*ϕ0spectral            # Surface air temperature
    # adjust the mean value (spectral coefficient 1,1) with Tabs_ref
    T0[1,1] += Complex(√2*Tabs_ref)

    # Stratosphere, set the first spectral coefficient (=mean value)
    # in two uppermost levels (i.e. 1,2) for lapse rate = 0
    #TODO why √2?
    Tabs[1,1,1] = Complex(√2*Tabs_top)
    Tabs[1,1,2] = Complex(√2*Tabs_top)

    # Temperature at tropospheric levels
    for k in 3:nlev
        Tabs[:,:,k] = T0*σ_full[k]^(R*γ_g)
    end

    # PRESSURE
    # Set logp - logarithm of surface pressure consistent with temperature profile
    @unpack nlon, nlat, p_ref = P
    logp_ref = log(p_ref)             # logarithm of reference surface pressure
    logp0 = zeros(nlon, nlat)         # logarithm of surface pressure by grid point

    for j in 1:nlat
        for i in 1:nlon
            logp0[i,j] = logp_ref + log(1.0 - γ_g*ϕ0[i,j]/Tabs_ref)/(R*γ_g)
        end
    end

    logp = spectral(T.(logp0),G)    # logarithm of surface pressure in spectral space
    truncate!(logp,G.spectral.trunc)# spectral truncation

    # SPECIFIC HUMIDITY
    @unpack es_ref, rh_ref, hshum, hscale = P
    qref = rh_ref*0.622*es_ref      # reference specific humidity
    qexp = hscale/hshum             # ratio of scale heights
    q = zeros(nlon,nlat)            # specific humidity by grid point

    # Specific humidity at the surface
    for j in 1:nlat
        for i in 1:nlon
            q[i,j] = qref*exp(qexp*logp0[i,j])
        end
    end

    humid0 = spectral(T.(q),G)
    truncate!(humid0,G.spectral.trunc)

    # Specific humidity at tropospheric levels (zero in the stratosphere[?])
    for k in 3:nlev
        humid[:,:,k] = humid0*σ_full[k]^qexp
    end

    # GEOPOTENTIAL [?]
    geopotential!(ϕ,ϕ0spectral,Tabs,G)

    # Print diagnostics from initial conditions
    check_global_mean_temperature(Tabs, P)

    # conversion to T happens here implicitly
    return Prognostics(ξ, D, Tabs, logp, humid, ϕ, ϕ0spectral)
end
