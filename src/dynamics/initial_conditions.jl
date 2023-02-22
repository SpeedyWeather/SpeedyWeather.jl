abstract type ZonalJet <: InitialConditions end             # Galewsky, Scott, Palvani, 2004
abstract type ZonalWind <: InitialConditions end            # Jablonowski & Williamson, 2006
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

    initial_conditions!(IC,progn,model)     # dispatch to the type of initial conditions IC

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
    homogeneous_temperature!(progn,model)
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

"""Initial conditions from Galewsky, 2004, Tellus"""
function initial_conditions!(   ::Type{ZonalJet},
                                progn::PrognosticVariables,
                                model::ShallowWater)

    @unpack latitude, width, umax = model.parameters.zonal_jet_coefs    # for jet
    @unpack perturb_lat, perturb_lon, perturb_xwidth,                   # for perturbation
        perturb_ywidth, perturb_height = model.parameters.zonal_jet_coefs

    θ₀ = (latitude-width)/360*2π    # southern boundary of jet [radians]
    θ₁ = (latitude+width)/360*2π    # northern boundary of jet
    eₙ = exp(-4/(θ₁-θ₀)^2)          # normalisation
    
    θ₂ = perturb_lat*2π/360         # perturbation latitude [radians]
    α = perturb_xwidth*2π/360       # zonal extent of interface perturbation [radians]
    β = perturb_ywidth*2π/360       # meridional extent of interface perturbation [radians]
    λ = perturb_lon*2π/360          # perturbation longitude [radians]

    @unpack radius_earth, rotation_earth, gravity = model.parameters

    # always create on F64 grid then convert to spectral and interpolate there
    Grid = FullGaussianGrid
    nlat_half = 64
    u_grid = zeros(Grid,nlat_half)
    η_grid = zeros(Grid,nlat_half)
    colats = RingGrids.get_colat(Grid,nlat_half)
    _,lons = RingGrids.get_colatlons(Grid,nlat_half)
    weights = FastGaussQuadrature.gausslegendre(2nlat_half)[2]
    η_sum = 0

    for (j,ring) in enumerate(eachring(u_grid,η_grid))
        θ = π/2 - colats[j]             # latitude in radians
        coslat⁻¹j = 1/cos(θ)
        f = 2rotation_earth*sin(θ)
        
        # velocity per latitude
        if θ₀ < θ < θ₁
            u_θ = umax/eₙ*exp(1/(θ-θ₀)/(θ-θ₁))  # u as in Galewsky, 2004
        else
            u_θ = 0
        end

        # integration for layer thickness h / interface height η
        w = weights[j]
        η_sum += 2w*(radius_earth*u_θ/gravity * (f + tan(θ)/radius_earth*u_θ))

        # lon-constant part of perturbation
        ηθ = perturb_height*cos(θ)*exp(-((θ₂-θ)/β)^2)

        # store in all longitudes
        for ij in ring
            u_grid[ij] = u_θ/radius_earth*coslat⁻¹j   # include scaling for curl!
            
            # calculate perturbation (possibly shifted in lon compared to Galewsky 2004)
            ϕ = lons[ij] - λ
            η_grid[ij] = η_sum + exp(-(ϕ/α)^2)*ηθ
        end
    end

    u = spectral(u_grid)
    η = spectral(η_grid)

    # interpolate in spectral space to desired resolution
    @unpack lmax,mmax = model.spectral_transform
    @unpack NF = model.parameters
    u = spectral_truncation(complex(NF),u,lmax+1,mmax)
    
    # get vorticity initial conditions from curl of u,v
    v = zero(u)     # meridional velocity zero for these initial conditions
    @unpack vor = progn.layers[end].leapfrog[1]
    curl!(vor,u,v,model.spectral_transform)

    # transform interface height η (use pres as prognostic variable) in spectral
    pres = progn.pres.leapfrog[1]
    copyto!(pres,η)
    spectral_truncation!(pres)
end

"""Initial conditions from Jablonowski and Williamson, 2006, QJR Meteorol. Soc"""
function initial_conditions!(   ::Type{ZonalWind},
                                progn::PrognosticVariables,
                                model::PrimitiveEquation)

    @unpack u₀, η₀, ΔT = model.parameters.zonal_wind_coefs
    @unpack perturb_lat, perturb_lon, perturb_uₚ, perturb_radius = model.parameters.zonal_wind_coefs
    @unpack temp_ref, R_dry, lapse_rate, gravity, pres_ref = model.parameters
    @unpack radius_earth, rotation_earth = model.parameters
    @unpack σ_tropopause = model.parameters
    @unpack σ_levels_full, Grid, nlat_half, nlev = model.geometry
    @unpack norm_sphere = model.spectral_transform

    φ, λ = model.geometry.latds, model.geometry.londs
    S = model.spectral_transform

    # VORTICITY
    ζ = zeros(Grid,nlat_half)   # relative vorticity
    D = zeros(Grid,nlat_half)   # divergence (perturbation only)

    for (k,layer) in enumerate(progn.layers)

        η = σ_levels_full[k]    # Jablonowski and Williamson use η for σ coordinates
        ηᵥ = (η - η₀)*π/2       # auxiliary variable for vertical coordinate

        # amplitude with height
        cos_ηᵥ = cos(ηᵥ)^(3/2)  # wind increases with height: 1 at model top, ~0.4 at surface

        for (ij,(φij,λij)) in enumerate(zip(φ,λ))
            sinφ = sind(φij)
            cosφ = cosd(φij)
            tanφ = tand(φij)

            # Jablonowski and Williamson, eq. (3) 
            ζ[ij] = -4u₀/radius_earth*cos_ηᵥ*sinφ*cosφ*(2 - 5sinφ^2) 

            # PERTURBATION
            sinφc = sind(perturb_lat)       # location of centre
            cosφc = cosd(perturb_lat)
            λc = perturb_lon
            R = radius_earth*perturb_radius # spatial extent of perturbation

            # great circle distance to perturbation
            X = sinφc*sinφ + cosφc*cosφ*cosd(λij-λc)
            X_norm = 1/sqrt(1-X^2)
            r = radius_earth*acos(X)
            exp_decay = exp(-(r/R)^2)

            # Jablonowski and Williamson, eq. (12) 
            ζ[ij] += perturb_uₚ/radius_earth*exp_decay*
                (tanφ - 2*(radius_earth/R)^2*acos(X)*X_norm*(sinφc*cosφ - cosφc*sinφ*cosd(λij-λc)))
            
            # Jablonowski and Williamson, eq. (13) 
            D[ij] -= 2perturb_uₚ*radius_earth/R^2 * exp_decay * acos(X) * X_norm * cosφc*sind(λij-λc)
        end

        @unpack vor = layer.leapfrog[1]
        spectral!(vor,ζ,S)
        spectral_truncation!(vor)
    end

    # TEMPERATURE
    Tη = zero(σ_levels_full)
    Γ = lapse_rate/1000                         # from [K/km] to [K/m]

    # vertical profile
    for k in 1:nlev
        σ = σ_levels_full[k]
        Tη[k] = temp_ref*σ^(R_dry*Γ/gravity)    # Jablonowski and Williamson eq. 4

        if σ < σ_tropopause
            Tη[k] += ΔT*(σ_tropopause-σ)^5      # Jablonowski and Williamson eq. 5
        end
    end

    T = zeros(Grid,nlat_half)   # temperature
    aΩ = radius_earth*rotation_earth

    for (k,layer) in enumerate(progn.layers)

        η = σ_levels_full[k]    # Jablonowski and Williamson use η for σ coordinates
        ηᵥ = (η - η₀)*π/2       # auxiliary variable for vertical coordinate

        # amplitudes with height
        A1 = 3/4*η*π*u₀/R_dry*sin(ηᵥ)* sqrt(cos(ηᵥ))
        A2 = 2u₀*cos(ηᵥ)^(3/2)

        for (ij,φij) in enumerate(φ)
            sinφ = sind(φij)
            cosφ = cosd(φij)

             # Jablonowski and Williamson, eq. (6) 
            T[ij] = Tη[k] + A1*((-2sinφ^6*(cosφ^2 + 1/3) + 10/63)*A2 + (8/5*cosφ^3*(sinφ^2 + 2/3) - π/4)*aΩ)
        end

        @unpack temp = layer.leapfrog[1]
        spectral!(temp,T,S)
        spectral_truncation!(temp)
    end

    # PRESSURE (constant everywhere)
    lnp₀ = log(pres_ref*100)        # logarithm of reference surface pressure, *100 for [hPa] to [Pa]
    progn.pres.leapfrog[1][1] = norm_sphere*lnp₀
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

function homogeneous_temperature!(  progn::PrognosticVariables,
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

function pressure_on_orography!(progn::PrognosticVariables,
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
    lnp₀ = log(pres_ref*100)        # logarithm of reference surface pressure, *100 for [hPa] to [Pa]
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