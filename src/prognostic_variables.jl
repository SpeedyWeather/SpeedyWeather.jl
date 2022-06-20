"""Struct holding the prognostic spectral variables."""
struct PrognosticVariables{NF<:AbstractFloat}
    # variables are lmax x mmax x nleapfrog x nlev
    # except for pres which is lmax x mmax x nleapfrog
    vor     ::Array{Complex{NF},4}      # Vorticity of horizontal wind field
    div     ::Array{Complex{NF},4}      # Divergence of horizontal wind field
    temp    ::Array{Complex{NF},4}      # Absolute temperature [K]
    pres    ::Array{Complex{NF},3}      # Logarithm of surface pressure [log(Pa)]
    humid   ::Array{Complex{NF},4}      # Specific humidity [g/kg]
end

"""Initialize prognostic variables from rest or restart from file."""
function initial_conditions(M::Union{BarotropicModel,ShallowWaterModel})

    @unpack initial_conditions = M.parameters

    if initial_conditions == :rest
        progn = initialize_from_rest(M)

    elseif initial_conditions == :barotropic_vorticity
        progn = initialize_from_rest(M)

        P = M.parameters    # unpack and rename
        G = M.geospectral

        @unpack nlon, nlat, nlev = G.geometry
        @unpack latd, lon, coslat, sinlat, radius_earth = G.geometry
        @unpack lmax, mmax = G.spectral_transform

        # zonal wind
        u_grid1 = @. (25*coslat - 30*coslat^3 + 300*sinlat^2*coslat^6)/coslat+100
        u_grid = repeat(u_grid1',nlon,1)
        u = zeros(Complex{P.NF},lmax+1,mmax+1)
        u = spectral!(u,u_grid,G.spectral_transform)
        ζ = gradient_latitude(u,G.spectral_transform,one_more_l=false,flipsign=true)
        progn.vor[:,:,1,1] .= ζ/radius_earth

        # zonal wave perturbation
        A = 1e-4
        m = 6
        θ0 = 45
        θw = 10

        ζp = convert.(P.NF,A/2*cos.(m*lon)) * convert.(P.NF,coslat .* exp.(-((latd .- θ0)/θw).^2))'
        progn.vor[:,:,1,1] .+= spectral(ζp,G.spectral_transform)

        for k in 2:nlev
            progn.vor[:,:,1,k] .= progn.vor[:,:,1,1]
        end

        # make it less symmetric
        progn.vor[15,1:14,1,:] .+= 3e-6*randn(Complex{P.NF},14,nlev)
    
    elseif initial_conditions == :barotropic_divergence
        progn = initialize_from_rest(M)

        P = M.parameters    # unpack and rename
        G = M.geospectral

        @unpack nlon, nlat, nlev = G.geometry
        @unpack latd, lon, coslat, sinlat, radius_earth = G.geometry
        @unpack lmax, mmax = G.spectral_transform

        progn.div[5,4,1,1] = 2e-5
        progn.vor[15,1:14,1,:] .+= 3e-6*randn(ComplexF64,14,nlev)

    elseif initial_conditions == :restart
        progn = initialize_from_file(M)         # TODO this is not implemented yet
    else
        throw(error("Incorrect initialization option, $initial_conditions given."))
    end

    # SCALING
    progn.vor .*= M.geospectral.geometry.radius_earth
    progn.div .*= M.geospectral.geometry.radius_earth

    return progn
end

"""Initialize a PrognosticVariables struct for an atmosphere at rest. No winds,
hence zero vorticity and divergence, but temperature, pressure and humidity are
initialised """
function initialize_from_rest(M::BarotropicModel)

    @unpack nlev = M.geospectral.geometry
    @unpack lmax, mmax = M.geospectral.spectral_transform
    nleapfrog = 2

    # conversion to type NF later when creating a PrognosticVariables struct
    vor     = zeros(Complex{Float64},lmax+1,mmax+1,nleapfrog,nlev)  # vorticity

    # dummy arrays for the rest, not used in this ModelSetup
    div     = zeros(Complex{Float64},1,1,1,1)
    temp    = zeros(Complex{Float64},1,1,1,1)
    pres    = zeros(Complex{Float64},1,1,1)
    humid   = zeros(Complex{Float64},1,1,1,1)

    # conversion to NF happens here
    @unpack NF = M.parameters
    return PrognosticVariables{NF}(vor,div,temp,pres,humid)
end

function initialize_from_rest(M::ShallowWaterModel)

    @unpack nlev = M.geospectral.geometry
    @unpack lmax, mmax = M.geospectral.spectral_transform
    nleapfrog = 2

    # conversion to type NF later when creating a PrognosticVariables struct
    vor     = zeros(Complex{Float64},lmax+1,mmax+1,nleapfrog,nlev)  # vorticity
    div     = zeros(Complex{Float64},lmax+1,mmax+1,nleapfrog,nlev)  # divergence
    pres    = zeros(Complex{Float64},lmax+1,mmax+1,nleapfrog)       # interface displacement

    # dummy arrays for the rest, not used in this ModelSetup
    temp    = zeros(Complex{Float64},1,1,1,1)
    humid   = zeros(Complex{Float64},1,1,1,1)

    # conversion to NF happens here
    @unpack NF = M.parameters
    return PrognosticVariables{NF}(vor,div,temp,pres,humid)
end

"""Initialize a PrognosticVariables struct for an atmosphere at rest. No winds,
hence zero vorticity and divergence, but temperature, pressure and humidity are
initialised """
function initialize_from_rest(M::PrimitiveEquationModel)

    P = M.parameters    # unpack and rename
    G = M.geospectral
    B = M.boundaries

    @unpack nlev = G.geometry
    @unpack lmax, mmax = G.spectral_transform
    nleapfrog = 2

    # conversion to type NF later when creating a PrognosticVariables struct
    # one more degree l than order m for recursion in meridional gradient
    vor     = zeros(Complex{Float64},lmax+1,mmax+1,nleapfrog,nlev)  # vorticity
    div     = zeros(Complex{Float64},lmax+1,mmax+1,nleapfrog,nlev)  # divergence
    temp    = zeros(Complex{Float64},lmax+1,mmax+1,nleapfrog,nlev)  # absolute Temperature
    pres    = zeros(Complex{Float64},lmax+1,mmax+1,nleapfrog)       # logarithm of surface pressure
    humid   = zeros(Complex{Float64},lmax+1,mmax+1,nleapfrog,nlev)  # specific humidity

    # initialize only the first leapfrog index
    temp_lf1 = view(temp,:,:,1,:)
    pres_lf1 = view(pres,:,:,1)
    humid_lf1 = view(humid,:,:,1,:)

    initialize_temperature!(temp_lf1,P,B,G)                    # temperature from lapse rates    
    pres_grid = initialize_pressure!(pres_lf1,P,B,G)  # pressure from temperature profile
    initialize_humidity!(humid_lf1,pres_grid,P,G)          # specific humidity from pressure

    # conversion to NF happens here
    @unpack NF = M.parameters
    return PrognosticVariables{NF}(vor,div,temp,pres,humid)
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
function initialize_pressure!(  pres::AbstractMatrix{Complex{NF}},      # logarithm of surface pressure
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
    pres_grid = zeros(nlon, nlat)               # logarithm of surface pressure by grid point

    for j in 1:nlat
        for i in 1:nlon
            pres_grid[i,j] = log_pres_ref + 
                log(1 - lapse_rate_scaled*geopot_surf_grid[i,j]/temp_ref)/(R_gas*lapse_rate_scaled)
        end
    end

    # convert to spectral space
    spectral!(pres,pres_grid,SpectralTransform(NF,nlon,nlat,P.trunc,P.radius_earth,true))
    # spectral_truncation!(pres,P.trunc)      # set lmax+1 row to zero
    return pres_grid                       # return grid for use in initialize_humidity!
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
    # spectral_truncation!(humid_surf,P.trunc)                # set the lmax+1 row to zero

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