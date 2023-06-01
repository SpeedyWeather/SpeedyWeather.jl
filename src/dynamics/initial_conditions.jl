# default initial conditions by model
initial_conditions_default(::Type{<:Barotropic}) = StartWithRandomVorticity()
initial_conditions_default(::Type{<:ShallowWater}) = ZonalJet()
initial_conditions_default(::Type{<:PrimitiveEquation}) = ZonalWind()

"""
$(TYPEDSIGNATURES)
Allocate the prognostic variables and then set to initial conditions."""
function initial_conditions(model::Model) where Model
    (;spectral_grid) = model
    progn = allocate(PrognosticVariables,spectral_grid,Model)   # allocate variables in any case
    IC = model.initial_conditions                               # initial conditions struct
    initial_conditions!(progn,IC,model)                         # dispatch to initial conditions
    return progn
end

"""
$(TYPEDSIGNATURES)"""
function allocate(
    ::Type{PrognosticVariables},
    spectral_grid::SpectralGrid,
    ::Type{Model},
) where {Model<:ModelSetup}

    (;NF,trunc,nlev) = spectral_grid
    return zeros(PrognosticVariables{NF},Model,trunc,nlev)
end

Base.@kwdef struct StartFromRest <: InitialConditions 
    pressure_on_orography::Bool = false
end

function initial_conditions!(   progn::PrognosticVariables,
                                initial_conditions::StartFromRest,
                                model::ModelSetup)
    return nothing          # everything remains zero 
end

function initial_conditions!(   progn::PrognosticVariables,
                                initial_conditions::StartFromRest,
                                model::PrimitiveEquation)
    homogeneous_temperature!(progn,model)
    initial_conditions.pressure_on_orography && pressure_on_orography!(progn,model)
    # TODO initialise humidity
end

"""Start with random vorticity as initial conditions
$(TYPEDFIELDS)"""
Base.@kwdef struct StartWithVorticity <: InitialConditions
    "Power of the spectral distribution k^power"
    power::Float64 = -3

    "(approximate) amplitude in [1/s], used as standard deviation of spherical harmonic coefficients"
    amplitude::Float64 = 1e-5
end

"""
$(TYPEDSIGNATURES)
Start with random vorticity as initial conditions"""
function initial_conditions!(   progn::PrognosticVariables{NF},
                                initial_conditions::StartWithRandomVorticity,
                                model::ModelSetup) where NF

    lmax = progn.trunc+1
    power = initial_conditions.power + 1    # +1 as power is summed of orders m
    ξ = randn(Complex{NF},lmax,lmax)*convert(NF,initial_conditions.amplitude)

    for progn_layer in progn.layers
        for m in 1:lmax
            for l in m:lmax
                progn_layer.timesteps[1].vor[l,m] = ξ[l,m]*l^power
            end
        end
        # don't perturb l=m=0 mode to have zero mean
        progn_layer.timesteps[1].vor[1] = 0
    end
end

"""
$(TYPEDSIGNATURES)
Create a struct that contains all parameters for the Galewsky et al, 2004 zonal jet
intitial conditions for the shallow water model. Default values as in Galewsky.
$(TYPEDFIELDS)"""
Base.@kwdef struct ZonalJet <: InitialConditions
    "jet latitude [˚N]"
    latitude::Float64 = 45
    
    "jet width [˚], default ≈ 19.29˚"
    width::Float64 = (1/4-1/7)*180

    "jet maximum velocity [m/s]"
    umax::Float64 = 80
    
    "perturbation latitude [˚N], position in jet by default"
    perturb_lat::Float64 = latitude
    
    "perturbation longitude [˚E]"
    perturb_lon::Float64 = 270
    
    "perturbation zonal extent [˚], default ≈ 19.1˚"
    perturb_xwidth::Float64 = 1/3*360/2π

    "perturbation meridinoal extent [˚], default ≈ 3.8˚"
    perturb_ywidth::Float64 = 1/15*360/2π
    
    "perturbation amplitude [m]"
    perturb_height::Float64 = 120
end

function Base.show(io::IO,IC::InitialConditions)
    print(io,"$(typeof(IC)) <: InitialConditions:")
    for key in propertynames(IC)
        val = getfield(IC,key)
        print(io,"\n $key::$(typeof(val)) = $val")
    end
end

"""
$(TYPEDSIGNATURES)
Initial conditions from Galewsky, 2004, Tellus"""
function initial_conditions!(   progn::PrognosticVariables,
                                initial_conditions::ZonalJet,
                                model::ShallowWater)

    (;latitude, width, umax) = initial_conditions               # for jet
    (;perturb_lat, perturb_lon, perturb_xwidth,                 # for perturbation
        perturb_ywidth, perturb_height) = initial_conditions

    θ₀ = (latitude-width)/360*2π    # southern boundary of jet [radians]
    θ₁ = (latitude+width)/360*2π    # northern boundary of jet
    eₙ = exp(-4/(θ₁-θ₀)^2)          # normalisation
    
    θ₂ = perturb_lat*2π/360         # perturbation latitude [radians]
    α = perturb_xwidth*2π/360       # zonal extent of interface perturbation [radians]
    β = perturb_ywidth*2π/360       # meridional extent of interface perturbation [radians]
    λ = perturb_lon*2π/360          # perturbation longitude [radians]

    (;rotation, gravity) = model.planet
    (;radius) = model.spectral_grid

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
    (;NF) = model.spectral_grid
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
$(TYPEDSIGNATURES)
Create a struct that contains all parameters for the Jablonowski and Williamson, 2006
intitial conditions for the primitive equation model. Default values as in Jablonowski.
$(TYPEDFIELDS)"""
Base.@kwdef struct ZonalWind <: InitialConditions
    "conversion from σ to Jablonowski's ηᵥ-coordinates"
    η₀::Float64 = 0.252
    
    "max amplitude of zonal wind [m/s]"
    u₀::Float64 = 35

    # PERTURBATION
    "perturbation centred at [˚N]"
    perturb_lat::Float64 = 40

    "perturbation centred at [˚E]"
    perturb_lon::Float64 = 20

    "perturbation strength [m/s]"
    perturb_uₚ::Float64 = 1
    
    "radius of Gaussian perturbation in units of Earth's radius [1]"
    perturb_radius::Float64 = 1/10

    # TERMPERATURE
    "temperature difference used for stratospheric lapse rate [K], Jablonowski uses ΔT = 4.8e5 [K]"
    ΔT::Float64 = 0                 

    "minimum temperature [K] of profile"
    Tmin::Float64 = 200

    # PRESSURE
    "initialize pressure given the `atmosphere.lapse_rate` on orography?"
    pressure_on_orography::Bool = false
end

"""
$(TYPEDSIGNATURES)
Initial conditions from Jablonowski and Williamson, 2006, QJR Meteorol. Soc"""
function initial_conditions!(   progn::PrognosticVariables{NF},
                                initial_conditions::ZonalWind,
                                model::PrimitiveEquation) where NF

    (;u₀, η₀, ΔT, Tmin, pressure_on_orography) = initial_conditions
    (;perturb_lat, perturb_lon, perturb_uₚ, perturb_radius) = initial_conditions
    (;temp_ref, R_dry, lapse_rate, pres_ref, σ_tropopause) = model.atmosphere
    (;radius, Grid, nlat_half, nlev) = model.spectral_grid
    (;rotation, gravity) = model.planet
    (;σ_levels_full) = model.geometry
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

    lnpₛ = ones(Grid{NF},nlat_half)
    lnpₛ .= pressure_on_orography ? pressure_on_orography!(progn,model) : lnp₀
    
    # HUMIDITY
    initialize_humidity!(progn,lnpₛ,model)
end

"""
Restart from a previous SpeedyWeather.jl simulation via the restart file restart.jld2
Applies interpolation in the horizontal but not in the vertical. restart.jld2 is
identified by
$(TYPEDFIELDS)"""
Base.@kwdef struct StartFromFile <: InitialConditions
    "path for restart file"
    path::String = pwd()

    "`run_id` of restart file in `run_????/restart.jld2`"
    id::Union{String,Int} = 1
end

"""
$(TYPEDSIGNATURES)
Restart from a previous SpeedyWeather.jl simulation via the restart file restart.jld2
Applies interpolation in the horizontal but not in the vertical."""
function initial_conditions!(   progn_new::PrognosticVariables,
                                initial_conditions::StartFromFile,
                                model::ModelSetup)

    (; path, id ) = initial_conditions

    restart_file = jldopen(joinpath(path,string("run_",run_id_to_string(id)),"restart.jld2"))
    progn_old = restart_file["prognostic_variables"]
    # version = restart_file["version"]             # currently unused
    model.clock.time = restart_file["time"]         # synchronize clocks
    return copy!(progn_new, progn_old)
end

function homogeneous_temperature!(  progn::PrognosticVariables,
                                    model::PrimitiveEquation)
    (; geopot_surf ) = model.orography       # spectral surface geopotential [m²/s²] (orography*gravity)

    # temp_ref:     Reference absolute T [K] at surface z = 0, constant lapse rate
    # temp_top:     Reference absolute T in the stratosphere [K], lapse rate = 0
    # lapse_rate:   Reference temperature lapse rate -dT/dz [K/km]
    # gravity:      Gravitational acceleration [m/s^2]
    # R_dry:        Specific gas constant for dry air [J/kg/K]
    (;temp_ref, temp_top, lapse_rate, R_dry, σ_tropopause) = model.atmosphere
    (;gravity) = model.planet
    (;nlev, σ_levels_full) = model.geometry         
    (; norm_sphere ) = model.spectral_transform # normalization of the l=m=0 spherical harmonic
    n_stratosphere_levels = findfirst(σ->σ>=σ_tropopause,σ_levels_full)

    # Lapse rate scaled by gravity [K/m / (m²/s²)]
    Γg⁻¹ = lapse_rate/gravity/1000                      # /1000 for lapse rate [K/km] → [K/m]

    # SURFACE TEMPERATURE (store in k = nlev, but it's actually surface, i.e. k=nlev+1/2)
    # overwrite with lowermost layer further down
    temp_surf = progn.layers[end].timesteps[1].temp     # spectral temperature at k=nlev+1/2
    temp_surf[1] = norm_sphere*temp_ref                 # set global mean surface temperature
    for lm in eachharmonic(geopot_surf,temp_surf)
        temp_surf[lm] -= Γg⁻¹*geopot_surf[lm]           # lower temperature for higher mountains
    end

    # TROPOPAUSE/STRATOSPHERE set the l=m=0 spectral coefficient (=mean value) only
    # in uppermost levels (default: k=1,2) for lapse rate = 0
    for k in 1:n_stratosphere_levels
        temp = progn.layers[k].timesteps[1].temp
        temp[1] = norm_sphere*temp_top
    end

    # TROPOSPHERE use lapserate and vertical coordinate σ for profile
    for k in n_stratosphere_levels+1:nlev               # k=nlev overwrites the surface temperature
                                                        # with lowermost layer temperature
        temp = progn.layers[k].timesteps[1].temp
        σₖᴿ = σ_levels_full[k]^(R_dry*Γg⁻¹)             # from hydrostatic equation

        for lm in eachharmonic(temp,temp_surf)
            temp[lm] = temp_surf[lm]*σₖᴿ
        end
    end
end

"""
$(TYPEDSIGNATURES)
Initialize surface pressure on orography by integrating the
hydrostatic equation with the reference temperature lapse rate."""
function pressure_on_orography!(progn::PrognosticVariables,
                                model::PrimitiveEquation)
    # temp_ref:     Reference absolute T [K] at surface z = 0, constant lapse rate
    # lapse_rate:   Reference temperature lapse rate -dT/dz [K/km]
    # gravity:      Gravitational acceleration [m/s^2]
    # R:            Specific gas constant for dry air [J/kg/K]
    # pres_ref:     Reference surface pressure [hPa]
    (; temp_ref, lapse_rate, pres_ref, R_dry ) = model.atmosphere
    (; gravity ) = model.planet
    (; orography ) = model.orography # orography on the grid

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
    spectral_truncation!(lnp)       # set lmax+1 row to zero
    return lnp_grid                 # return grid for use in initialize_humidity!
end

function initialize_humidity!(  progn::PrognosticVariables,
                                pres_surf_grid::AbstractGrid,
                                model::PrimitiveDry)
    return nothing
end

function initialize_humidity!(  progn::PrognosticVariables,
                                pres_surf_grid::AbstractGrid,
                                model::PrimitiveWet)

    # reference saturation water vapour pressure [Pa]
    # relative humidity reference [1]
    (;water_pres_ref, relhumid_ref, R_dry, R_vapour, pres_ref) = model.atmosphere
    gas_ratio = R_dry/R_vapour
    humid_ref = relhumid_ref*gas_ratio*water_pres_ref   # reference specific humidity [Pa]

    # ratio of scale heights [1], scale height [km], scale height for spec humidity [km]     
    (;scale_height, scale_height_humid, σ_tropopause) = model.atmosphere
    scale_height_ratio = scale_height/scale_height_humid

    (;nlev, σ_levels_full) = model.geometry
    n_stratosphere_levels = findfirst(σ->σ>=σ_tropopause,σ_levels_full)

    # Specific humidity at the surface (grid space)
    humid_surf_grid = zero(pres_surf_grid)
    # @. humid_surf_grid = humid_ref*(exp(pres_surf_grid)/(pres_ref*100))^scale_height_ratio
    q_ref = 20e-3       # kg/kg at the surface
    @. humid_surf_grid .= q_ref
    scale_coslat²!(humid_surf_grid,model.geometry)

    humid_surf = spectral(humid_surf_grid,model.spectral_transform)
    spectral_truncation!(humid_surf)

    # Specific humidity at tropospheric levels (stratospheric humidity remains zero)
    a = model.spectral_transform.norm_sphere
    for k in n_stratosphere_levels+1:nlev
        for lm in eachharmonic(humid_surf,progn.layers[1].timesteps[1].humid)
            progn.layers[k].timesteps[1].humid[lm] = humid_surf[lm]*σ_levels_full[k]^scale_height_ratio
        end
    end
end