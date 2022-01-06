"""Struct holding the prognostic spectral variables."""
struct PrognosticVariables{NF<:AbstractFloat}
    vor         ::Array{Complex{NF},3}      # Vorticity of horizontal wind field
    div         ::Array{Complex{NF},3}      # Divergence of horizontal wind field
    temp        ::Array{Complex{NF},3}      # Absolute temperature [K]
    pres_surf   ::Array{Complex{NF},2}      # Logarithm of surface pressure [log(Pa)]
    humid       ::Array{Complex{NF},3}      # Specific humidity [g/kg]
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

    # conversion to type NF later when creating a PrognosticVariables struct
    vor         = zeros(Complex{Float64},mx,nx,nlev)  # vorticity
    div         = zeros(Complex{Float64},mx,nx,nlev)  # divergence
    temp        = zeros(Complex{Float64},mx,nx,nlev)  # absolute Temperature
    pres_surf   = zeros(Complex{Float64},mx,nx)       # logarithm of surface pressure
    humid       = zeros(Complex{Float64},mx,nx,nlev)  # specific humidity

    initialize_temperature!(temp,P,B,G)                     # temperature from lapse rates    
    pres_surf_grid = initialize_pressure!(pres_surf,P,B,G)  # pressure from temperature profile
    initialize_humidity!(humid,pres_surf_grid,P,G)          # specific humidity from pressure

    # conversion to NF happens here implicitly
    return PrognosticVariables{NF}(vor,div,temp,pres_surf,humid)
end

"""Initialize spectral temperature from surface absolute temperature and constant
lapse rate (troposphere) and zero lapse rate (stratosphere)."""
function initialize_temperature!(   temp::AbstractArray{Complex{NF},3}, # spectral temperature in 3D
                                    P::Params,                          # Parameters struct
                                    B::Boundaries{NF},                  # Boundaries struct
                                    G::GeoSpectral{NF}                  # Geospectral struct
                                    ) where {NF<:AbstractFloat}         # number format NF

    @unpack geopot_surf = B     # spectral surface geopotential [m²/s²]

    # temp_ref:     Reference absolute T [K] at surface z = 0, constant lapse rate
    # temp_top:     Reference absolute T in the stratosphere [K], lapse rate = 0
    # lapse_rate:   Reference temperature lapse rate -dT/dz [K/km]
    # gravity:      Gravitational acceleration [m/s^2]
    # R:            Specific gas constant for dry air [J/kg/K]
    @unpack temp_ref, temp_top, lapse_rate, gravity, R = P

    lapse_rate_scaled = lape_rate/gravity/1000      # Lapse rate scaled by gravity [K/m / (m²/s²)]

    temp_surf = -lapse_rate_scaled*geopot_surf      # spectral surface air temperature from orography and lapse rate
    temp_surf[1,1] += Complex(√2*temp_ref)          # adjust mean value (spectral coefficient 1,1) with temp_ref

    # Stratosphere, set the first spectral coefficient (=mean value)
    # in two uppermost levels (i.e. k=1,2) for lapse rate = 0
    #TODO why √2? (normalisation?)
    temp[1,1,1] = Complex(√2*temp_top)
    temp[1,1,2] = Complex(√2*temp_top)

    # Temperature at tropospheric levels
    @unpack σ_full = G.geometry

    for k in 3:nlev
        temp[:,:,k] .= temp_surf*σ_full[k]^(R*lapse_rate_scaled)
    end
end

"""Initialize the logarithm of surface pressure `logp0` consistent with temperature profile."""
function initiliaze_pressure!(  pres_surf::AbstractArray{Complex{NF},2},    # logarithm of surface pressure
                                P::Params,                                  # Parameters struct
                                B::Boundaries{NF},                          # Boundaries struct
                                G::GeoSpectral{NF}                          # Geospectral struct
                                ) where {NF<:AbstractFloat}                 # number format NF
    
    @unpack nlon, nlat = P

    # temp_ref:     Reference absolute T [K] at surface z = 0, constant lapse rate
    # temp_top:     Reference absolute T in the stratosphere [K], lapse rate = 0
    # lapse_rate:   Reference temperature lapse rate -dT/dz [K/km]
    # gravity:      Gravitational acceleration [m/s^2]
    # R:            Specific gas constant for dry air [J/kg/K]
    # pres_ref:     Reference surface pressure [hPa]
    @unpack temp_ref, temp_top, lapse_rate, gravity, pres_ref, R = P
    @unpack geopot_surf = B                     # spectral surface geopotential
    geopot_surf_grid = gridded(geopot_surf,G)   # convert to grid-point space

    lapse_rate_scaled = lape_rate/gravity/1000  # Lapse rate scaled by gravity [K/m / (m²/s²)]
    log_pres_ref = log(pres_ref)                # logarithm of reference surface pressure
    pres_surf_grid = zeros(nlon, nlat)          # logarithm of surface pressure by grid point

    for j in 1:nlat
        for i in 1:nlon
            pres_surf_grid[i,j] = log_pres_ref + 
                log(1 - lapse_rate_scaled*geopot_surf_grid[i,j]/temp_ref)/(R*lapse_rate_scaled)
        end
    end

    spectral!(pres_surf,pres_surf_grid,G)   # convert to spectral space
    spectral_truncation!(pres_surf,G)       # truncate in spectral space

    return pres_surf_grid                   # return grid for use in initialize_humidity!
end

"""Initialize specific humidity in spectral space."""
function initialize_humidity!(  humid::AbstractArray{Complex{NF},3},    # spectral specific humidity
                                pres_surf_grid::AbstractArray{NF,2},    # logarithm of surface pressure (grid space)
                                P::Params,                              # Parameters struct
                                G::GeoSpectral{NF}                      # Geospectral struct
                                ) where {NF<:AbstractFloat}             # number format NF

    mx,nx,nlev = size(humid)
    @unpack nlon, nlat = P

    # reference saturation water vapour pressure [Pa]
    # relative humidity reference [1]
    @unpack water_pres_ref, relhumid_ref = P
    humid_ref = relhumid_ref*0.622*water_pres_ref   # reference specific humidity [Pa]

    @unpack scale_height, scale_height_humid = P            # scale height [km], scale height for spec humidity [km]   
    scale_height_ratio = scale_height/scale_height_humid    # ratio of scale heights [1]

    # Specific humidity at the surface (grid space)
    humid_surf_grid = humid_ref*exp.(scale_height_ratio*pres_surf_grid)

    # Convert to spectral space and truncate
    humid_surf = spectral(humid_surf_grid,G)
    spectral_truncation!(humid_surf,G)

    # stratospheric humidity zero
    humid[:,:,1:2] .= 0

    # Specific humidity at tropospheric levels
    for k in 3:nlev
        for j in 1:nx
            for i in 1:mx
                humid[i,j,k] = humid0[i,j]*σ_full[k]^scale_height_ratio
            end
        end
    end
end