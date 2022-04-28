"""Struct holding the prognostic spectral variables."""
struct PrognosticVariables{NF<:AbstractFloat}
    # variables are lmax x mmax x nleapfrog x nlev
    # except for pres_surf lmax x mmax x nleapfrog
    vor         ::Array{Complex{NF},4}      # Vorticity of horizontal wind field
    div         ::Array{Complex{NF},4}      # Divergence of horizontal wind field
    temp        ::Array{Complex{NF},4}      # Absolute temperature [K]
    pres_surf   ::Array{Complex{NF},3}      # Logarithm of surface pressure [log(Pa)]
    humid       ::Array{Complex{NF},4}      # Specific humidity [g/kg]
end

"""Initialize prognostic variables from rest or restart from file."""
function initial_conditions(    P::Parameters,      # Parameter struct
                                B::Boundaries,      # Boundaries struct
                                G::GeoSpectral)     # GeoSpectral struct

    @unpack initial_conditions = P

    if initial_conditions == :rest
        progn = initialize_from_rest(P,B,G)
    elseif initial_conditions == :restart
        progn = initialize_from_file(P,B,G)
    else
        throw(error("Incorrect initialization option, $initial_conditions given."))
    end

    return progn
end

"""Initialize a PrognosticVariables struct for an atmosphere at rest. No winds,
hence zero vorticity and divergence, but temperature, pressure and humidity are
initialised """
function initialize_from_rest(  P::Parameters,
                                B::Boundaries,
                                G::GeoSpectral)

    @unpack nlev = G.geometry
    @unpack lmax, mmax = G.spectral_transform
    nleapfrog = 2

    # conversion to type NF later when creating a PrognosticVariables struct
    # one more degree l than order m for recursion in meridional gradient
    vor         = zeros(Complex{Float64},lmax+1,mmax+1,nleapfrog,nlev)  # vorticity
    div         = zeros(Complex{Float64},lmax+1,mmax+1,nleapfrog,nlev)  # divergence
    temp        = zeros(Complex{Float64},lmax+1,mmax+1,nleapfrog,nlev)  # absolute Temperature
    pres_surf   = zeros(Complex{Float64},lmax+1,mmax+1,nleapfrog)       # logarithm of surface pressure
    humid       = zeros(Complex{Float64},lmax+1,mmax+1,nleapfrog,nlev)  # specific humidity

    # initialize only the first leapfrog index
    temp_lf1 = view(temp,:,:,1,:)
    pres_surf_lf1 = view(pres_surf,:,:,1)
    humid_lf1 = view(humid,:,:,1,:)

    initialize_temperature!(temp_lf1 ,P,B,G)                    # temperature from lapse rates    
    pres_surf_grid = initialize_pressure!(pres_surf_lf1,P,B,G)  # pressure from temperature profile
    initialize_humidity!(humid_lf1,pres_surf_grid,P,G)          # specific humidity from pressure

    # conversion to NF happens here implicitly
    return PrognosticVariables{P.NF}(vor,div,temp,pres_surf,humid)
end

"""Initialize spectral temperature from surface absolute temperature and constant
lapse rate (troposphere) and zero lapse rate (stratosphere)."""
function initialize_temperature!(   temp::AbstractArray{Complex{NF},3},    # spectral temperature in 3D
                                    P::Parameters,                         # Parameters struct
                                    B::Boundaries,                         # Boundaries struct
                                    G::GeoSpectral                         # Geospectral struct
                                    ) where NF

    lmax,mmax,nlev = size(temp)     # number of vertical levels nlev
    lmax, mmax = lmax-1, mmax-1     # get 0-based max degree l, order m of spherical harmonics
    @unpack geopot_surf = B         # spectral surface geopotential [m²/s²]

    # temp_ref:     Reference absolute T [K] at surface z = 0, constant lapse rate
    # temp_top:     Reference absolute T in the stratosphere [K], lapse rate = 0
    # lapse_rate:   Reference temperature lapse rate -dT/dz [K/km]
    # gravity:      Gravitational acceleration [m/s^2]
    # R_gas:        Specific gas constant for dry air [J/kg/K]
    @unpack temp_ref, temp_top, lapse_rate, gravity, R_gas = P
    @unpack n_stratosphere_levels = P               # number of vertical levels used for stratosphere
    @unpack norm_sphere = G.spectral_transform      # normalization of the l=m=0 spherical harmonic

    lapse_rate_scaled = lapse_rate/gravity/1000     # Lapse rate scaled by gravity [K/m / (m²/s²)]
    temp_surf = -lapse_rate_scaled*geopot_surf      # spectral surface air temperature from orography and lapse rate
    temp_surf[1,1] += norm_sphere*temp_ref          # adjust mean value (spectral coefficient 1,1) with temp_ref

    # Stratosphere, set the first spectral coefficient (=mean value)
    # in uppermost levels (default: k=1,2) for lapse rate = 0
    for k in 1:n_stratosphere_levels
        temp[1,1,k] = norm_sphere*temp_top
    end

    # Temperature at tropospheric levels
    @unpack σ_levels_full = G.geometry

    for k in n_stratosphere_levels+1:nlev
        for m in 1:mmax+1
            for l in m:lmax+1
                temp[l,m,k] = temp_surf[l,m]*σ_levels_full[k]^(R_gas*lapse_rate_scaled)
            end
        end
    end
end

"""Initialize the logarithm of surface pressure `logp0` consistent with temperature profile."""
function initialize_pressure!(  pres_surf::AbstractMatrix{Complex{NF}}, # logarithm of surface pressure
                                P::Parameters,                          # Parameters struct
                                B::Boundaries,                          # Boundaries struct
                                G::GeoSpectral) where NF                # Geospectral struct
    
    @unpack nlon, nlat = P
    S = G.spectral_transform

    # temp_ref:     Reference absolute T [K] at surface z = 0, constant lapse rate
    # temp_top:     Reference absolute T in the stratosphere [K], lapse rate = 0
    # lapse_rate:   Reference temperature lapse rate -dT/dz [K/km]
    # gravity:      Gravitational acceleration [m/s^2]
    # R:            Specific gas constant for dry air [J/kg/K]
    # pres_ref:     Reference surface pressure [hPa]
    @unpack temp_ref, temp_top, lapse_rate, gravity, pres_ref, R_gas = P
    @unpack geopot_surf = B                     # spectral surface geopotential
    geopot_surf_grid = gridded(geopot_surf,S)   # convert to grid-point space

    lapse_rate_scaled = lapse_rate/gravity/1000 # Lapse rate scaled by gravity [K/m / (m²/s²)]
    log_pres_ref = log(pres_ref)                # logarithm of reference surface pressure
    pres_surf_grid = zeros(nlon, nlat)          # logarithm of surface pressure by grid point

    for j in 1:nlat
        for i in 1:nlon
            pres_surf_grid[i,j] = log_pres_ref + 
                log(1 - lapse_rate_scaled*geopot_surf_grid[i,j]/temp_ref)/(R_gas*lapse_rate_scaled)
        end
    end

    # convert to spectral space
    spectral!(pres_surf,pres_surf_grid,SpectralTransform(NF,nlon,nlat,P.trunc,P.radius_earth,true))
    spectral_truncation!(pres_surf,P.trunc)     # set lmax+1 row to zero
    return pres_surf_grid                       # return grid for use in initialize_humidity!
end

"""Initialize specific humidity in spectral space."""
function initialize_humidity!(  humid::AbstractArray{Complex{NF},3},# spectral specific humidity
                                pres_surf_grid::AbstractMatrix,     # log of surf pressure (grid space)
                                P::Parameters,                      # Parameters struct
                                G::GeoSpectral) where NF            # Geospectral struct

    lmax,mmax,nlev = size(humid)    # of size lmax+1, mmax+1, nlev
    lmax, mmax = lmax-1, mmax-1     # hence correct with -1 for 0-based l,m
    @unpack nlon, nlat, n_stratosphere_levels = P
    @unpack σ_levels_full = G.geometry

    # reference saturation water vapour pressure [Pa]
    # relative humidity reference [1]
    @unpack water_pres_ref, relhumid_ref = P
    humid_ref = relhumid_ref*0.622*water_pres_ref   # reference specific humidity [Pa]

    # scale height [km], scale height for spec humidity [km]
    @unpack scale_height, scale_height_humid = P            
    scale_height_ratio = scale_height/scale_height_humid    # ratio of scale heights [1]

    # Specific humidity at the surface (grid space)
    humid_surf_grid = humid_ref*exp.(scale_height_ratio*pres_surf_grid)
    humid_surf = spectral(humid_surf_grid,one_more_l=true)
    spectral_truncation!(humid_surf,P.trunc)                # set the lmax+1 row to zero

    # stratospheric humidity zero
    fill!(view(humid,:,:,1:n_stratosphere_levels),0)

    # Specific humidity at tropospheric levels
    for k in n_stratosphere_levels+1:nlev
        for m in 1:mmax+1
            for l in m:lmax+1
                humid[l,m,k] = humid_surf[l,m]*σ_levels_full[k]^scale_height_ratio
            end
        end
    end
end