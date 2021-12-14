"""Struct holding the prognostic spectral variables."""
struct PrognosticVariables{NF<:AbstractFloat}
    vor     ::Array{Complex{NF},3}      # Vorticity of horizontal wind field
    div     ::Array{Complex{NF},3}      # Divergence of horizontal wind field
    Tabs    ::Array{Complex{NF},3}      # Absolute temperature [K]
    logp0   ::Array{Complex{NF},2}      # Log of surface pressure [log(Pa)]
    humid   ::Array{Complex{NF},3}      # Specific humidity [g/kg]
    geopot  ::Array{Complex{NF},3}      # Atmospheric geopotential [m²/s²]
end

"""Initialize prognostic variables from rest or restart from file."""
function initial_conditions(    P::Params,                      # Parameter struct
                                B::Boundaries{NF},              # Boundaries struct
                                G::GeoSpectral{NF}              # GeoSpectral struct
                                ) where {NF<:AbstractFloat}     # number format NF

    @unpack initial_conditions = P

    if initial_conditions == :rest
        ProgVars = initialize_from_rest(P,B,G)
    elseif initial_conditions == :restart
        ProgVars = initialize_from_file(P,B,G)
    else
        throw(error("Incorrect initialization option, $initial_conditions given."))
    end

    return ProgVars
end

"""Initialize a PrognosticVariables struct for an atmosphere at rest. No winds,
hence zero vorticity and divergence, but temperature, pressure and humidity are
initialised """
function initialize_from_rest(  P::Params,
                                B::Boundaries{NF},
                                G::GeoSpectral{NF}
                                ) where {NF<:AbstractFloat}

    @unpack nlev = G.geometry
    @unpack mx, nx = G.spectral
    @unpack ϕ0 = B

    # conversion to type NF later when creating a PrognosticVariables struct
    vor    = zeros(Complex{Float64}, mx, nx, nlev)  # Vorticity
    div    = zeros(Complex{Float64}, mx, nx, nlev)  # Divergence
    Tabs   = zeros(Complex{Float64}, mx, nx, nlev)  # Absolute Temperature
    logp0  = zeros(Complex{Float64}, mx, nx)        # logarithm of surface pressure
    humid  = zeros(Complex{Float64}, mx, nx, nlev)  # specific humidity
    geopot = zeros(Complex{Float64}, mx, nx, nlev)  # geopotential

    initialize_temperature!(Tabs,P,B,G)             # temperature from lapse rates    
    logp0_grid = initialize_pressure!(logp0,P,B,G)  # pressure from temperature profile
    initialize_humidity!(humid,logp0_grid,P,G)      # specific humidity from pressure
    geopotential!(geopot,ϕ0spectral,Tabs,G)         # geopotential from surface geopotential

    # conversion to NF happens here implicitly
    return PrognosticVariables{NF}(vor, div, Tabs, logp0, humid, geopot)
end

"""Initialize spectral temperature from surface absolute temperature and constant
lapse rate (troposphere) and zero lapse rate (stratosphere)."""
function initialize_temperature!(   Tabs::AbstractArray{Complex{NF},3}, # spectral temperature in 3D
                                    P::Params,                          # Parameters struct
                                    B::Boundaries{NF},                  # Boundaries struct
                                    G::GeoSpectral{NF}                  # Geospectral struct
                                    ) where {NF<:AbstractFloat}         # number format NF

    # Compute spectral surface geopotential
    @unpack ϕ0 = B
    ϕ0_spectral = spectral(ϕ0,G)

    # Tabs_ref: Reference absolute T [K] at surface z = 0, constant lapse rate
    # Tabs_top: Reference absolute T in the stratosphere [K], lapse rate = 0
    @unpack Tabs_ref, Tabs_top, γ, g, R = P

    γ_g = γ/g/1000                  # Lapse rate [K/m] scaled by gravity
    T0 = -γ_g*ϕ0_spectral            # Surface air temperature
    T0[1,1] += Complex(√2*Tabs_ref) # adjust mean value (spectral coefficient 1,1) with Tabs_ref

    # Stratosphere, set the first spectral coefficient (=mean value)
    # in two uppermost levels (i.e. k=1,2) for lapse rate = 0
    #TODO why √2? (normalisation?)
    Tabs[1,1,1] = Complex(√2*Tabs_top)
    Tabs[1,1,2] = Complex(√2*Tabs_top)

    # Temperature at tropospheric levels
    @unpack σ_full = G.geometry

    for k in 3:nlev
        Tabs[:,:,k] = T0*σ_full[k]^(R*γ_g)
    end
end

"""Initialize the logarithm of surface pressure `logp0` consistent with temperature profile."""
function initiliaze_pressure!(  logp0::AbstractArray{Complex{NF},2},    # logarithm of surface pressure
                                P::Params,                              # Parameters struct
                                B::Boundaries{NF},                      # Boundaries struct
                                G::GeoSpectral{NF}                      # Geospectral struct
                                ) where {NF<:AbstractFloat}             # number format NF
    
    @unpack nlon, nlat, p0_ref = P
    @unpack Tabs_ref, Tabs_top, γ, g, R = P
    @unpack ϕ0 = B                  # load surface geopotential

    γ_g = γ/g/1000                  # Lapse rate [K/m] scaled by gravity
    logp0_ref = log(p0_ref)         # logarithm of reference surface pressure
    logp0_grid = zeros(nlon, nlat)  # logarithm of surface pressure by grid point

    for j in 1:nlat
        for i in 1:nlon
            logp0_grid[i,j] = logp0_ref + log(1 - γ_g*ϕ0[i,j]/Tabs_ref)/(R*γ_g)
        end
    end

    spectral!(logp0,logp0_grid,G)   # convert to spectral space
    spectral_truncation!(logp0,G)   # smooth in spectral space

    return logp0_grid
end

"""Initialize specific humidity in spectral space."""
function initialize_humidity!(  humid::AbstractArray{Complex{NF},3},    # specific humidity
                                logp0_grid::AbstractArray{NF,2},        # logarithm of surface pressure
                                P::Params,                              # Parameters struct
                                G::GeoSpectral{NF}                      # Geospectral struct
                                ) where {NF<:AbstractFloat}             # number format NF

    mx,nx,nlev = size(humid)
    @unpack nlon, nlat = P
    @boundscheck nlon, nlat == size(logp0_grid) || throw(BoundsError())

    @unpack es_ref, rh_ref = P      # reference saturation water vapour pressure [Pa] + relative humidity [1]
    @unpack hscale, hshum = P       # scale height [km], scale height for spec humidity [km]
    q_ref = rh_ref*0.622*es_ref     # reference specific humidity [Pa]
    q_exp = hscale/hshum            # ratio of scale heights [1]

    # Specific humidity at the surface (grid space)
    humid0_grid = q_ref*exp.(q_exp*logp0_grid)

    humid0 = spectral(humid0_grid,G)
    spectral_truncation!(humid0,G)

    # Specific humidity at tropospheric levels (zero in the stratosphere[?])
    for k in 3:nlev
        for j in 1:nx
            for i in 1:mx
                humid[i,j,k] = humid0[i,j]*σ_full[k]^q_exp
            end
        end
    end
end