abstract type ZonalWind <: InitialConditions end
abstract type StartFromRest <: InitialConditions end
abstract type StartFromFile <: InitialConditions end
abstract type StartWithVorticity <: InitialConditions end

"""
    prognostic_variables = initial_conditions(M::ModelSetup)

Initialize the prognostic variables from rest, with some initial vorticity,
or restart from file."""
function initial_conditions(model::ModelSetup)

    progn = allocate_prognostic_variables(model)    # allocate variables in any case
    IC = model.parameters.initial_conditions        # type of initial conditions

    initial_conditions!(IC,progn,model)     # dispatch to the type of initial conditions

    # SCALING: we use vorticity*radius,divergence*radius in the dynamical core
    @unpack radius_earth = model.geometry
    scale!(progn,:vor,radius_earth)
    scale!(progn,:div,radius_earth)

    return progn
end

function initial_conditions!(   ::Type{StartFromRest},
                                progn::PrognosticVariables,
                                model::ModelSetup)
    return nothing
end

function initial_conditions!(   ::Type{StartFromRest},
                                progn::PrognosticVariables,
                                model::PrimitiveEquation)
    homogenuous_temperature!(progn,model)
    pres_grid = pressure_on_orography!(progn,model)
    # TODO initialise humidity
end
    
function initial_conditions!(   ::Type{StartWithVorticity},
                                progn::PrognosticVariables,
                                model::ModelSetup)
    @unpack radius_earth = model.geometry

    mmax = min(15,progn.mmax+1)     # perturb only larger modes
    lmax = min(15,progn.lmax+1)

    ξ = randn(Complex{model.parameters.NF},mmax,lmax)

    for progn_layer in progn.layers
        
        # zonal wind
        progn_layer.leapfrog[1].vor[4,1]  =  80/radius_earth    
        progn_layer.leapfrog[1].vor[6,1]  = -160/radius_earth
        progn_layer.leapfrog[1].vor[8,1]  =  80/radius_earth

        # perturbation
        for m in 2:min(15,progn.mmax+1)
            for l in m:min(15,progn.lmax+1)
                progn_layer.leapfrog[1].vor[l,m] = 10/radius_earth*ξ[l,m]
            end
        end
    end
end

function initial_conditions!(   ::Type{ZonalWind},
                                progn::PrognosticVariables,
                                model::ModelSetup)

#         φ = M.geometry.latds

#         # the -4u₀R⁻¹ sinφ cosφ (2 - 5sin²φ) factor
#         A = zero(ϕ)


#         for k in layers
#             η = σ_levels_full[k]
#             σᵥ = (σ .- 0.252)*π/2


#         ηₛ = 1
#         η₀ = 0.252
#         ηᵥ = (ηₛ-η₀)*π/2
#         u₀ = 35
#         s = u₀*cos(ηᵥ)^(3/2)
#         aΩ = radius_earth*rotation_earth
#         g⁻¹ = inv(gravity)
#         ϕ = G.latds
    
#         for ij in eachindex(ϕ,orography.data)
#             sinϕ = sind(ϕ[ij])
#             cosϕ = cosd(ϕ[ij])
#             orography[ij] = g⁻¹*s*(s*(-2*sinϕ^6*(cosϕ^2 + 1/3) + 10/63) + (8/5*cosϕ^3*(sinϕ^2 + 2/3) - π/4)*aΩ)
#         end
end

function allocate_prognostic_variables(model::ModelSetup) 

    @unpack NF = model.parameters
    @unpack nlev = model.geometry
    @unpack lmax, mmax = model.spectral_transform

    return zeros(PrognosticVariables{NF},model,lmax,mmax,nlev)
end

function initial_conditions!(   ::Type{StartFromFile},
                                progn_new::PrognosticVariables,
                                model::ModelSetup)

    @unpack restart_path, restart_id = model.parameters

    restart_file = jldopen(joinpath(restart_path,string("run-",run_id_string(restart_id)),"restart.jld2"))
    progn_old = restart_file["prognostic_variables"]
    version = restart_file["version"]   # currently unused, TODO check for compat with version
    time = restart_file["time"]         # currently unused

    var_names = propertynames(progn_old.layers[1].leapfrog[1])

    for var_name in var_names
        if has(progn_new, var_name) 
            var = get_var(progn_old, var_name) 
            set_var!(progn_new, var_name, var)
        end
    end 
    pres = get_pressure(progn_old)
    set_pressure!(progn_new, pres)

    return progn_new
end

function homogenous_temperature!(   progn::PrognosticVariables,
                                    model::PrimitiveEquation)

    P = model.parameters
    B = model.boundaries
    G = model.geometry
    S = model.spectral_transform

    @unpack geopot_surf = B.orography       # spectral surface geopotential [m²/s²] (orography*gravity)

    # temp_ref:     Reference absolute T [K] at surface z = 0, constant lapse rate
    # temp_top:     Reference absolute T in the stratosphere [K], lapse rate = 0
    # lapse_rate:   Reference temperature lapse rate -dT/dz [K/km]
    # gravity:      Gravitational acceleration [m/s^2]
    # R_dry:        Specific gas constant for dry air [J/kg/K]
    @unpack temp_ref, temp_top, lapse_rate, gravity, R_dry = P
    @unpack n_stratosphere_levels, nlev = G     # number of vertical levels used for stratosphere
    @unpack norm_sphere = S                     # normalization of the l=m=0 spherical harmonic

    Γg⁻¹ = lapse_rate/gravity/1000              # Lapse rate scaled by gravity [K/m / (m²/s²)]

    # SURFACE TEMPERATURE
    temp_surf = progn.layers[end].leapfrog[1].temp  # spectral temperature at k=nlev=surface
    temp_surf[1] = norm_sphere*temp_ref             # set global mean surface temperature
    for lm in eachharmonic(geopot_surf,temp_surf)
        temp_surf[lm] -= Γg⁻¹*geopot_surf[lm]       # lower temperature for higher mountains
    end

    # TROPOPAUSE/STRATOSPHERE set the l=m=0 spectral coefficient (=mean value) only
    # in uppermost levels (default: k=1,2) for lapse rate = 0
    for k in 1:n_stratosphere_levels
        temp = progn.layers[k].leapfrog[1].temp
        temp[1] = norm_sphere*temp_top
    end

    # TROPOSPHERE use lapserate and vertical coordinate σ for profile
    for k in n_stratosphere_levels+1:nlev
        temp = progn.layers[k].leapfrog[1].temp
        σₖᴿ = G.σ_levels_full[k]^(R_dry*Γg⁻¹)   # TODO reference

        for lm in eachharmonic(temp,temp_surf)
            temp[lm] = temp_surf[lm]*σₖᴿ
        end
    end
end

function initialize_pressure!(  progn::PrognosticVariables,
                                model::PrimitiveEquation)
    
    P = model.parameters
    B = model.boundaries
    G = model.geometry
    S = model.spectral_transform

    @unpack Grid,nlat_half = G
    @unpack lmax,mmax = S
    
    # temp_ref:     Reference absolute T [K] at surface z = 0, constant lapse rate
    # temp_top:     Reference absolute T in the stratosphere [K], lapse rate = 0
    # lapse_rate:   Reference temperature lapse rate -dT/dz [K/km]
    # gravity:      Gravitational acceleration [m/s^2]
    # R:            Specific gas constant for dry air [J/kg/K]
    # pres_ref:     Reference surface pressure [hPa]
    @unpack temp_ref, temp_top, lapse_rate, gravity, pres_ref, R_dry = P
    @unpack orography = B.orography # orography on the grid

    Γ = lapse_rate/1000             # Lapse rate [K/km] -> [K/m]
    lnp₀ = log(pres_ref)            # logarithm of reference surface pressure
    lnp_grid = zero(orography)      # allocate log surface pressure on grid

    RΓg⁻¹ = R_dry*Γ/gravity         # for convenience
    ΓT⁻¹ = Γ/temp_ref           

    for ij in eachgridpoint(lnp_grid,orography)
        lnp_grid[ij] = lnp₀ + log(1 - ΓT⁻¹*orography[ij])/RΓg⁻¹
    end

    lnp = progn.pres.leapfrog[1]
    spectral!(lnp,lnp_grid,S)
    spectral_truncation!(lnp,lmax)  # set lmax+1 row to zero

    return lnp_grid                 # return grid for use in initialize_humidity!
end

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