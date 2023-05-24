"""
    prognostic_variables = initial_conditions(M::ModelSetup)

Allocate the prognostic variables and then set to initial conditions."""
function initial_conditions(model::ModelSetup)

    progn = allocate_prognostic_variables(model)    # allocate variables in any case
    IC = model.parameters.initial_conditions        # initial conditions struct
    initial_conditions!(IC,progn,model)             # dispatch to initial conditions
    return progn
end

function allocate_prognostic_variables(model::ModelSetup) 

    (;NF) = model.parameters
    (;nlev) = model.geometry
    (;lmax, mmax) = model.spectral_transform

    return zeros(PrognosticVariables{NF},model,lmax,mmax,nlev)
end

struct StartFromRest <: InitialConditions end

function initial_conditions!(   ::StartFromRest,
                                progn::PrognosticVariables,
                                model::ModelSetup)
    return nothing
end

function initial_conditions!(   ::StartFromRest,
                                progn::PrognosticVariables,
                                model::PrimitiveEquation)
    homogeneous_temperature!(progn,model)
    model.parameters.pressure_on_orography && pressure_on_orography!(progn,model)
    # TODO initialise humidity
end

struct StartWithVorticity <: InitialConditions end

function initial_conditions!(   ::StartWithVorticity,
                                progn::PrognosticVariables,
                                model::ModelSetup)

    (;radius) = model.geometry

    mmax = min(15,progn.mmax+1)     # perturb only larger modes
    lmax = min(15,progn.lmax+1)

    ξ = randn(Complex{model.parameters.NF},mmax,lmax)

    for progn_layer in progn.layers
        
        # zonal wind
        progn_layer.timesteps[1].vor[4,1]  =  80/radius
        progn_layer.timesteps[1].vor[6,1]  = -160/radius
        progn_layer.timesteps[1].vor[8,1]  =  80/radius

        # perturbation
        for m in 2:min(15,progn.mmax+1)
            for l in m:min(15,progn.lmax+1)
                progn_layer.timesteps[1].vor[l,m] = 10/radius*ξ[l,m]
            end
        end
    end
end

"""
    Z = ZonalJet(;kwargs...) <: InitialConditions

Create a struct that contains all parameters for the Galewsky et al, 2004 zonal jet
intitial conditions for the shallow water model. Default values as in Galewsky."""
Base.@kwdef struct ZonalJet <: InitialConditions
    # jet
    latitude = 45               # degrees north [˚N]
    width = (1/4-1/7)*180       # ≈ 19.29˚ as in Galewsky et al, 2004 
    umax = 80                   # [m/s]
    
    # perturbation
    perturb_lat = latitude          # [˚N], position in jet by default
    perturb_lon = 270               # [˚E]
    perturb_xwidth = 1/3*360/2π     # ≈ 19.1˚E zonal extent [˚E]
    perturb_ywidth = 1/15*360/2π    # ≈ 3.8˚N meridional extent [˚N]
    perturb_height = 120            # amplitude [m]
end

"""Initial conditions from Galewsky, 2004, Tellus"""
function initial_conditions!(   coefs::ZonalJet,
                                progn::PrognosticVariables,
                                model::ShallowWater)

    (;latitude, width, umax) = coefs               # for jet
    (;perturb_lat, perturb_lon, perturb_xwidth,   # for perturbation
        perturb_ywidth, perturb_height) = coefs

    θ₀ = (latitude-width)/360*2π    # southern boundary of jet [radians]
    θ₁ = (latitude+width)/360*2π    # northern boundary of jet
    eₙ = exp(-4/(θ₁-θ₀)^2)          # normalisation
    
    θ₂ = perturb_lat*2π/360         # perturbation latitude [radians]
    α = perturb_xwidth*2π/360       # zonal extent of interface perturbation [radians]
    β = perturb_ywidth*2π/360       # meridional extent of interface perturbation [radians]
    λ = perturb_lon*2π/360          # perturbation longitude [radians]

    (;radius, rotation, gravity) = model.parameters.planet

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
        f = 2rotation*sin(θ)
        
        # velocity per latitude
        if θ₀ < θ < θ₁
            u_θ = umax/eₙ*exp(1/(θ-θ₀)/(θ-θ₁))  # u as in Galewsky, 2004
        else
            u_θ = 0
        end

        # integration for layer thickness h / interface height η
        w = weights[j]
        η_sum += 2w*(radius*u_θ/gravity * (f + tan(θ)/radius*u_θ))

        # lon-constant part of perturbation
        ηθ = perturb_height*cos(θ)*exp(-((θ₂-θ)/β)^2)

        # store in all longitudes
        for ij in ring
            u_grid[ij] = u_θ/radius*coslat⁻¹j   # include scaling for curl!
            
            # calculate perturbation (possibly shifted in lon compared to Galewsky 2004)
            ϕ = lons[ij] - λ
            η_grid[ij] = η_sum + exp(-(ϕ/α)^2)*ηθ
        end
    end

    u = spectral(u_grid)
    η = spectral(η_grid)

    # interpolate in spectral space to desired resolution
    (;lmax,mmax) = model.spectral_transform
    (;NF) = model.parameters
    u = spectral_truncation(complex(NF),u,lmax+1,mmax)
    
    # get vorticity initial conditions from curl of u,v
    v = zero(u)     # meridional velocity zero for these initial conditions
    (;vor) = progn.layers[end].timesteps[1]
    curl!(vor,u,v,model.spectral_transform)

    # transform interface height η (use pres as prognostic variable) in spectral
    (;pres) = progn.surface.timesteps[1]
    copyto!(pres,η)
    spectral_truncation!(pres)
end

"""
    Z = ZonalWind(;kwargs...) <: InitialConditions

Create a struct that contains all parameters for the Jablonowski and Williamson, 2006
intitial conditions for the primitive equation model. Default values as in Jablonowski."""
Base.@kwdef struct ZonalWind <: InitialConditions
    
    # vertical
    η₀ = 0.252                  # conversion from σ to Jablonowski's ηᵥ-coordinates
    u₀ = 35                     # max amplitude of zonal wind [m/s]
    
    # perturbation
    perturb_lat = 40            # Gaussian profile perturbation centred at [˚N]
    perturb_lon = 20            # and [˚E]
    perturb_uₚ = 1              # strength of perturbation [m/s]
    perturb_radius = 1/10       # radius of Gaussian perturbation in units of Earth's radius [1]

    # temperature
    ΔT = 0                      # temperature difference used for stratospheric lapse rate [K]
                                # Jablonowski uses ΔT = 4.8e5 [K]
    Tmin = 200                  # minimum temperature [K] of profile
end

"""Initial conditions from Jablonowski and Williamson, 2006, QJR Meteorol. Soc"""
function initial_conditions!(   coefs::ZonalWind,
                                progn::PrognosticVariables{NF},
                                model::PrimitiveEquation) where NF

    (;u₀, η₀, ΔT, Tmin) = coefs
    (;perturb_lat, perturb_lon, perturb_uₚ, perturb_radius) = coefs
    (;temp_ref, R_dry, lapse_rate, pres_ref) = model.parameters
    (;radius, rotation, gravity) = model.parameters.planet
    (;σ_tropopause, pressure_on_orography) = model.parameters
    (;σ_levels_full, Grid, nlat_half, nlev) = model.geometry
    (;norm_sphere) = model.spectral_transform

    φ, λ = model.geometry.latds, model.geometry.londs
    S = model.spectral_transform

    # VORTICITY
    ζ = zeros(Grid{NF},nlat_half)   # relative vorticity
    D = zeros(Grid{NF},nlat_half)   # divergence (perturbation only)

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
            ζ[ij] = -4u₀/radius*cos_ηᵥ*sinφ*cosφ*(2 - 5sinφ^2) 

            # PERTURBATION
            sinφc = sind(perturb_lat)       # location of centre
            cosφc = cosd(perturb_lat)
            λc = perturb_lon
            R = radius*perturb_radius # spatial extent of perturbation

            # great circle distance to perturbation
            X = sinφc*sinφ + cosφc*cosφ*cosd(λij-λc)
            X_norm = 1/sqrt(1-X^2)
            r = radius*acos(X)
            exp_decay = exp(-(r/R)^2)

            # Jablonowski and Williamson, eq. (12) 
            ζ[ij] += perturb_uₚ/radius*exp_decay*
                (tanφ - 2*(radius/R)^2*acos(X)*X_norm*(sinφc*cosφ - cosφc*sinφ*cosd(λij-λc)))
            
            # Jablonowski and Williamson, eq. (13) 
            D[ij] = -2perturb_uₚ*radius/R^2 * exp_decay * acos(X) * X_norm * cosφc*sind(λij-λc)
        end

        (;vor,div) = layer.timesteps[1]
        spectral!(vor,ζ,S)
        spectral!(div,D,S)
        spectral_truncation!(vor)
        spectral_truncation!(div)
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

    Tη .= max.(Tη,Tmin)

    T = zeros(Grid{NF},nlat_half)   # temperature
    aΩ = radius*rotation

    for (k,layer) in enumerate(progn.layers)

        η = σ_levels_full[k]    # Jablonowski and Williamson use η for σ coordinates
        ηᵥ = (η - η₀)*π/2       # auxiliary variable for vertical coordinate

        # amplitudes with height
        A1 = 3/4*η*π*u₀/R_dry*sin(ηᵥ)*sqrt(cos(ηᵥ))
        A2 = 2u₀*cos(ηᵥ)^(3/2)

        for (ij,φij) in enumerate(φ)
            sinφ = sind(φij)
            cosφ = cosd(φij)

             # Jablonowski and Williamson, eq. (6) 
            T[ij] = Tη[k] + A1*((-2sinφ^6*(cosφ^2 + 1/3) + 10/63)*A2 + (8/5*cosφ^3*(sinφ^2 + 2/3) - π/4)*aΩ)
        end

        (;temp) = layer.timesteps[1]
        spectral!(temp,T,S)
        spectral_truncation!(temp)
    end

    # PRESSURE (constant everywhere)
    lnp₀ = log(pres_ref*100)        # logarithm of reference surface pressure, *100 for [hPa] to [Pa]
    progn.surface.timesteps[1].pres[1] = norm_sphere*lnp₀
    pressure_on_orography && pressure_on_orography!(progn,model)
end

struct StartFromFile <: InitialConditions end

function initial_conditions!(   ::StartFromFile,
                                progn_new::PrognosticVariables,
                                model::ModelSetup)

    (; restart_path, restart_id ) = model.parameters

    restart_file = jldopen(joinpath(restart_path,string("run-",run_id_string(restart_id)),"restart.jld2"))
    progn_old = restart_file["prognostic_variables"]
    version = restart_file["version"]   # currently unused, TODO check for compat with version
    time = restart_file["time"]         # currently unused

    return copy!(progn_new, progn_old)
end

function homogeneous_temperature!(  progn::PrognosticVariables,
                                    model::PrimitiveEquation)

    P = model.parameters
    B = model.boundaries
    G = model.geometry
    S = model.spectral_transform

    (; geopot_surf ) = B.orography       # spectral surface geopotential [m²/s²] (orography*gravity)

    # temp_ref:     Reference absolute T [K] at surface z = 0, constant lapse rate
    # temp_top:     Reference absolute T in the stratosphere [K], lapse rate = 0
    # lapse_rate:   Reference temperature lapse rate -dT/dz [K/km]
    # gravity:      Gravitational acceleration [m/s^2]
    # R_dry:        Specific gas constant for dry air [J/kg/K]
    (; temp_ref, temp_top, lapse_rate, R_dry ) = P
    (; gravity ) = P.planet
    (; n_stratosphere_levels, nlev ) = G     # number of vertical levels used for stratosphere
    (; norm_sphere ) = S                     # normalization of the l=m=0 spherical harmonic

    Γg⁻¹ = lapse_rate/gravity/1000              # Lapse rate scaled by gravity [K/m / (m²/s²)]

    # SURFACE TEMPERATURE
    temp_surf = progn.layers[end].timesteps[1].temp  # spectral temperature at k=nlev=surface
    temp_surf[1] = norm_sphere*temp_ref             # set global mean surface temperature
    for lm in eachharmonic(geopot_surf,temp_surf)
        temp_surf[lm] -= Γg⁻¹*geopot_surf[lm]       # lower temperature for higher mountains
    end

    # TROPOPAUSE/STRATOSPHERE set the l=m=0 spectral coefficient (=mean value) only
    # in uppermost levels (default: k=1,2) for lapse rate = 0
    for k in 1:n_stratosphere_levels
        temp = progn.layers[k].timesteps[1].temp
        temp[1] = norm_sphere*temp_top
    end

    # TROPOSPHERE use lapserate and vertical coordinate σ for profile
    for k in n_stratosphere_levels+1:nlev
        temp = progn.layers[k].timesteps[1].temp
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

    (; Grid,nlat_half ) = G
    (; lmax,mmax ) = S
    
    # temp_ref:     Reference absolute T [K] at surface z = 0, constant lapse rate
    # temp_top:     Reference absolute T in the stratosphere [K], lapse rate = 0
    # lapse_rate:   Reference temperature lapse rate -dT/dz [K/km]
    # gravity:      Gravitational acceleration [m/s^2]
    # R:            Specific gas constant for dry air [J/kg/K]
    # pres_ref:     Reference surface pressure [hPa]
    (; temp_ref, temp_top, lapse_rate, pres_ref, R_dry ) = P
    (; gravity ) = P.planet
    (; orography ) = B.orography # orography on the grid

    Γ = lapse_rate/1000             # Lapse rate [K/km] -> [K/m]
    lnp₀ = log(pres_ref*100)        # logarithm of reference surface pressure, *100 for [hPa] to [Pa]
    lnp_grid = zero(orography)      # allocate log surface pressure on grid

    RΓg⁻¹ = R_dry*Γ/gravity         # for convenience
    ΓT⁻¹ = Γ/temp_ref           

    for ij in eachgridpoint(lnp_grid,orography)
        lnp_grid[ij] = lnp₀ + log(1 - ΓT⁻¹*orography[ij])/RΓg⁻¹
    end

    lnp = progn.pres.timesteps[1]
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
#     (; nlon, nlat, n_stratosphere_levels ) = P
#     (; σ_levels_full ) = G.geometry

#     # reference saturation water vapour pressure [Pa]
#     # relative humidity reference [1]
#     (; water_pres_ref, relhumid_ref ) = P
#     humid_ref = relhumid_ref*0.622*water_pres_ref   # reference specific humidity [Pa]

#     # scale height [km], scale height for spec humidity [km]
#     (; scale_height, scale_height_humid ) = P            
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