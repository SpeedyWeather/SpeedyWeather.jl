"""
    prognostic_variables = initial_conditions(M::ModelSetup)

Initialize the prognostic variables from rest, with some initial vorticity,
or restart from file."""
function initial_conditions(M::ModelSetup)

    progn = initialize_from_rest(M)             # allocate variables in any case
    @unpack initial_conditions = M.parameters
    @unpack radius_earth = M.geometry

    if initial_conditions == :rest
        nothing 

    elseif initial_conditions == :barotropic_vorticity              # zonal wind with perturbation
        for progn_layer in progn.layers
            progn_layer.leapfrog[1].vor[4,1]  =  80/radius_earth    # zonal wind
            progn_layer.leapfrog[1].vor[6,1]  = -160/radius_earth
            progn_layer.leapfrog[1].vor[8,1]  =  80/radius_earth

            # perturbation
            for m in 2:min(15,progn.mmax+1)
                for l in m:min(15,progn.lmax+1)
                    progn_layer.leapfrog[1].vor[l,m] = 10/radius_earth*randn(ComplexF64)
                end
            end
        end

    elseif initial_conditions == :restart
        initialize_from_file!(progn,M)
    else
        throw(error("Incorrect initialization option, $initial_conditions given."))
    end

    # SCALING
    @unpack radius_earth = M.geometry
    scale!(progn,:vor,radius_earth)
    scale!(progn,:div,radius_earth)

    return progn
end

function initialize_from_rest(M::ModelSetup) 

    @unpack NF = M.parameters
    @unpack nlev = M.geometry
    @unpack lmax, mmax = M.spectral_transform

    return zeros(PrognosticVariables{NF},M,lmax,mmax,nlev)
end

"""Initialize a PrognosticVariables struct for an atmosphere at rest. No winds,
hence zero vorticity and divergence, but temperature, pressure and humidity are
initialised """
function initialize_from_rest(M::PrimitiveEquationModel)

    @unpack NF = M.parameters
    @unpack nlev = M.geometry
    @unpack lmax, mmax = M.spectral_transform

    progn = zeros(PrognosticVariables{NF},M,lmax,mmax,nlev)

    # TODO adjust temp and pres to rest
    # initialize_temperature!(temp_lf1,P,B,G)                    # temperature from lapse rates    
    # pres_grid = initialize_pressure!(pres_lf1,P,B,G)  # pressure from temperature profile
    # initialize_humidity!(humid_lf1,pres_grid,P,G)          # specific humidity from pressure

    return progn
end

function initialize_from_file!(progn_new::PrognosticVariables{NF},M::ModelSetup) where NF
    @unpack restart_path, restart_id = M.parameters
    restart_file = jldopen(joinpath(restart_path,@sprintf("run%04d",restart_id),"restart.jld2"))
    progn_old = restart_file["prognostic_variables"]
    version = restart_file["version"]   # currently unused, TODO check for compat with version
    time = restart_file["time"]         # currently unused

    var_names = [propertynames(progn_old.layers[1].leapfrog[1])]

    for var_name in var_names 
        var = get_var(progn_old, var_name) 
        set_var!(progn_new, var_name, var)
    end 

    return progn_new
end

# """Initialize spectral temperature from surface absolute temperature and constant
# lapse rate (troposphere) and zero lapse rate (stratosphere)."""
# function initialize_temperature!(   temp::AbstractArray{Complex{NF},3},    # spectral temperature in 3D
#                                     P::Parameters,                         # Parameters struct
#                                     B::Boundaries,                         # Boundaries struct
#                                     G::Geometry                            # Geospectral struct
#                                     S::SpectralTransform
#                                     ) where NF

#     lmax,mmax,nlev = size(temp)     # number of vertical levels nlev
#     lmax, mmax = lmax-1, mmax-1     # get 0-based max degree l, order m of spherical harmonics
#     @unpack geopot_surf = B         # spectral surface geopotential [m²/s²]

#     # temp_ref:     Reference absolute T [K] at surface z = 0, constant lapse rate
#     # temp_top:     Reference absolute T in the stratosphere [K], lapse rate = 0
#     # lapse_rate:   Reference temperature lapse rate -dT/dz [K/km]
#     # gravity:      Gravitational acceleration [m/s^2]
#     # R_gas:        Specific gas constant for dry air [J/kg/K]
#     @unpack temp_ref, temp_top, lapse_rate, gravity, R_gas = P
#     @unpack n_stratosphere_levels = P               # number of vertical levels used for stratosphere
#     @unpack norm_sphere = G.spectral_transform      # normalization of the l=m=0 spherical harmonic

#     lapse_rate_scaled = lapse_rate/gravity/1000     # Lapse rate scaled by gravity [K/m / (m²/s²)]
#     temp_surf = -lapse_rate_scaled*geopot_surf      # spectral surface air temperature from orography and lapse rate
#     temp_surf[1,1] += norm_sphere*temp_ref          # adjust mean value (spectral coefficient 1,1) with temp_ref

#     # Stratosphere, set the first spectral coefficient (=mean value)
#     # in uppermost levels (default: k=1,2) for lapse rate = 0
#     for k in 1:n_stratosphere_levels
#         temp[1,1,k] = norm_sphere*temp_top
#     end

#     # Temperature at tropospheric levels
#     @unpack σ_levels_full = G.geometry

#     for k in n_stratosphere_levels+1:nlev
#         for m in 1:mmax+1
#             for l in m:lmax+1
#                 temp[l,m,k] = temp_surf[l,m]*σ_levels_full[k]^(R_gas*lapse_rate_scaled)
#             end
#         end
#     end
# end

# """Initialize the logarithm of surface pressure `logp0` consistent with temperature profile."""
# function initialize_pressure!(  pres::AbstractMatrix{Complex{NF}},      # logarithm of surface pressure
#                                 P::Parameters,                          # Parameters struct
#                                 B::Boundaries,                          # Boundaries struct
#                                 G::Geomtry,
#                                 S::SpectralTransform
#                                 ) where NF                # Geospectral struct
    
#     @unpack nlon, nlat = P
#     S = G.spectral_transform

#     # temp_ref:     Reference absolute T [K] at surface z = 0, constant lapse rate
#     # temp_top:     Reference absolute T in the stratosphere [K], lapse rate = 0
#     # lapse_rate:   Reference temperature lapse rate -dT/dz [K/km]
#     # gravity:      Gravitational acceleration [m/s^2]
#     # R:            Specific gas constant for dry air [J/kg/K]
#     # pres_ref:     Reference surface pressure [hPa]
#     @unpack temp_ref, temp_top, lapse_rate, gravity, pres_ref, R_gas = P
#     @unpack geopot_surf = B                     # spectral surface geopotential
#     geopot_surf_grid = gridded(geopot_surf,S)   # convert to grid-point space

#     lapse_rate_scaled = lapse_rate/gravity/1000 # Lapse rate scaled by gravity [K/m / (m²/s²)]
#     log_pres_ref = log(pres_ref)                # logarithm of reference surface pressure
#     pres_grid = zeros(nlon, nlat)               # logarithm of surface pressure by grid point

#     for j in 1:nlat
#         for i in 1:nlon
#             pres_grid[i,j] = log_pres_ref + 
#                 log(1 - lapse_rate_scaled*geopot_surf_grid[i,j]/temp_ref)/(R_gas*lapse_rate_scaled)
#         end
#     end

#     # convert to spectral space
#     spectral!(pres,pres_grid,SpectralTransform(NF,nlon,nlat,P.trunc,P.radius_earth,true))
#     # spectral_truncation!(pres,P.trunc)      # set lmax+1 row to zero
#     return pres_grid                       # return grid for use in initialize_humidity!
# end

# """Initialize specific humidity in spectral space."""
# function initialize_humidity!(  humid::AbstractArray{Complex{NF},3},# spectral specific humidity
#                                 pres_surf_grid::AbstractMatrix,     # log of surf pressure (grid space)
#                                 P::Parameters,                      # Parameters struct
#                                 G::GeoSpectral) where NF            # Geospectral struct

#     lmax,mmax,nlev = size(humid)    # of size lmax+1, mmax+1, nlev
#     lmax, mmax = lmax-1, mmax-1     # hence correct with -1 for 0-based l,m
#     @unpack nlon, nlat, n_stratosphere_levels = P
#     @unpack σ_levels_full = G.geometry

#     # reference saturation water vapour pressure [Pa]
#     # relative humidity reference [1]
#     @unpack water_pres_ref, relhumid_ref = P
#     humid_ref = relhumid_ref*0.622*water_pres_ref   # reference specific humidity [Pa]

#     # scale height [km], scale height for spec humidity [km]
#     @unpack scale_height, scale_height_humid = P            
#     scale_height_ratio = scale_height/scale_height_humid    # ratio of scale heights [1]

#     # Specific humidity at the surface (grid space)
#     humid_surf_grid = humid_ref*exp.(scale_height_ratio*pres_surf_grid)
#     humid_surf = spectral(humid_surf_grid,one_more_l=true)
#     # spectral_truncation!(humid_surf,P.trunc)                # set the lmax+1 row to zero

#     # stratospheric humidity zero
#     fill!(view(humid,:,:,1:n_stratosphere_levels),0)

#     # Specific humidity at tropospheric levels
#     for k in n_stratosphere_levels+1:nlev
#         for m in 1:mmax+1
#             for l in m:lmax+1
#                 humid[l,m,k] = humid_surf[l,m]*σ_levels_full[k]^scale_height_ratio
#             end
#         end
#     end
# end