abstract type InitialConditions end

# EXPORT INITIAL CONDITIONS
export  StartFromFile,
        StartFromRest,
        ZonalJet,
        ZonalWind,
        StartWithRandomVorticity

Base.@kwdef struct StartFromRest <: InitialConditions 
    pressure_on_orography::Bool = false
end

function initialize!(   progn::PrognosticVariables,
                        initial_conditions::StartFromRest,
                        model::ModelSetup)
    return nothing          # everything remains zero 
end

function initialize!(   progn::PrognosticVariables,
                        initial_conditions::StartFromRest,
                        model::PrimitiveEquation)
    homogeneous_temperature!(progn,model)
    initial_conditions.pressure_on_orography && pressure_on_orography!(progn,model)
    # TODO initialise humidity
end

"""Start with random vorticity as initial conditions
$(TYPEDFIELDS)"""
Base.@kwdef mutable struct StartWithRandomVorticity <: InitialConditions
    "Power of the spectral distribution k^power"
    power::Float64 = -3

    "(approximate) amplitude in [1/s], used as standard deviation of spherical harmonic coefficients"
    amplitude::Float64 = 1e-5
end

"""
$(TYPEDSIGNATURES)
Start with random vorticity as initial conditions"""
function initialize!(   progn::PrognosticVariables{NF},
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
Base.@kwdef mutable struct ZonalJet <: InitialConditions
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
    println(io,"$(typeof(IC)) <: InitialConditions")
    keys = propertynames(IC)
    print_fields(io,IC,keys)
end

"""
$(TYPEDSIGNATURES)
Initial conditions from Galewsky, 2004, Tellus"""
function initialize!(   progn::PrognosticVariables,
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
    u = spectral_truncation(complex(NF),u,lmax,mmax)
    
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

    "Sigma coordinates of the tropopause [1]"
    σ_tropopause::Float64 = 0.2
    
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
    "reference dry-adiabtic lapse rate [K/m]"
    lapse_rate::Float64 = 5/1000

    "temperature difference used for stratospheric lapse rate [K], Jablonowski uses ΔT = 4.8e5 [K]"
    ΔT::Float64 = 0                 

    "minimum temperature [K] of profile"
    Tmin::Float64 = 200

    # PRESSURE
    "initialize pressure given the `atmosphere.lapse_rate` on orography?"
    pressure_on_orography::Bool = true
end

"""
$(TYPEDSIGNATURES)
Initial conditions from Jablonowski and Williamson, 2006, QJR Meteorol. Soc"""
function initialize!(   progn::PrognosticVariables{NF},
                        initial_conditions::ZonalWind,
                        model::PrimitiveEquation) where NF

    (;u₀, η₀, ΔT, Tmin, pressure_on_orography) = initial_conditions
    (;perturb_lat, perturb_lon, perturb_uₚ, perturb_radius) = initial_conditions
    (;lapse_rate, σ_tropopause) = initial_conditions
    (;temp_ref, R_dry, pres_ref) = model.atmosphere
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

    # vertical profile
    for k in 1:nlev
        σ = σ_levels_full[k]
        Tη[k] = temp_ref*σ^(R_dry*lapse_rate/gravity)   # Jablonowski and Williamson eq. 4

        if σ < σ_tropopause
            Tη[k] += ΔT*(σ_tropopause-σ)^5              # Jablonowski and Williamson eq. 5
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
    lnp₀ = log(pres_ref)            # logarithm of reference surface pressure [log(Pa)]
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
function initialize!(   progn_new::PrognosticVariables,
                        initial_conditions::StartFromFile,
                        model::ModelSetup)

    (; path, id ) = initial_conditions

    restart_file = jldopen(joinpath(path,string("run_",run_id_to_string(id)),"restart.jld2"))
    progn_old = restart_file["prognostic_variables"]
    # version = restart_file["version"]             # currently unused
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
    (;temp_ref, lapse_rate, R_dry) = model.atmosphere
    (;gravity) = model.planet
    (;nlev, σ_levels_full) = model.geometry         
    (;norm_sphere) = model.spectral_transform # normalization of the l=m=0 spherical harmonic

    # Lapse rate scaled by gravity [K/m / (m²/s²)]
    Γg⁻¹ = lapse_rate/gravity

    # SURFACE TEMPERATURE (store in k = nlev, but it's actually surface, i.e. k=nlev+1/2)
    # overwrite with lowermost layer further down
    temp_surf = progn.layers[end].timesteps[1].temp     # spectral temperature at k=nlev+1/2
    temp_surf[1] = norm_sphere*temp_ref                 # set global mean surface temperature
    for lm in eachharmonic(geopot_surf,temp_surf)
        temp_surf[lm] -= Γg⁻¹*geopot_surf[lm]           # lower temperature for higher mountains
    end

    # Use lapserate and vertical coordinate σ for profile
    for k in 1:nlev                                     # k=nlev overwrites the surface temperature
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
    # lapse_rate:   Reference temperature lapse rate -dT/dz [K/m]
    # gravity:      Gravitational acceleration [m/s^2]
    # R:            Specific gas constant for dry air [J/kg/K]
    # pres_ref:     Reference surface pressure [hPa]
    (; temp_ref, lapse_rate, pres_ref, R_dry ) = model.atmosphere
    (; gravity ) = model.planet
    (; orography ) = model.orography # orography on the grid

    lnp₀ = log(pres_ref)            # logarithm of reference surface pressure [log(Pa)]
    lnp_grid = zero(orography)      # allocate log surface pressure on grid

    RΓg⁻¹ = R_dry*lapse_rate/gravity         # for convenience
    ΓT⁻¹ = lapse_rate/temp_ref           

    for ij in eachgridpoint(lnp_grid,orography)
        lnp_grid[ij] = lnp₀ + log(1 - ΓT⁻¹*orography[ij])/RΓg⁻¹
    end

    lnp = progn.surface.timesteps[1].pres
    spectral!(lnp,lnp_grid,model.spectral_transform)
    spectral_truncation!(lnp)       # set lmax+1 row to zero
    return lnp_grid                 # return grid for use in initialize_humidity!
end

function initialize_humidity!(  progn::PrognosticVariables,
                                pres_surf_grid::AbstractGrid,
                                model::PrimitiveDry)
    return nothing
end

function initialize_humidity!(  
    progn::PrognosticVariables{NF},
    pres_surf_grid::AbstractGrid,
    model::PrimitiveWet
) where NF

    # TODO create a InitialHumidity <: InitialConditions to hold these parameters
    relhumid_ref::NF = 0.7
    scale_height_humid::NF = 2.5        # scale height for specific humidity [km]
    scale_height::NF = 7.5              # scale height for pressure [km]
    scale_height_ratio = scale_height/scale_height_humid

    (;nlev, σ_levels_full) = model.geometry

    # Specific humidity at the surface (grid space)
    temp_grid = gridded(progn.layers[end].timesteps[1].temp,model.spectral_transform)
    humid_surf_grid = zero(pres_surf_grid)
    for ij in eachgridpoint(humid_surf_grid)
        q_sat = saturation_humidity(temp_grid[ij],exp(pres_surf_grid[ij]),model.clausius_clapeyron)
        humid_surf_grid[ij] = relhumid_ref*q_sat
    end

    humid_surf = spectral(humid_surf_grid,model.spectral_transform)
    spectral_truncation!(humid_surf)

    # Specific humidity at levels above
    for k in 1:nlev
        for lm in eachharmonic(humid_surf)
            progn.layers[k].timesteps[1].humid[lm] = humid_surf[lm]*σ_levels_full[k]^scale_height_ratio
        end
    end
end

"""Parameters for random initial conditions for the interface displacement η
in the shallow water equations.
$(TYPEDFIELDS)"""
Base.@kwdef struct RandomWaves <: InitialConditions  
    # random interface displacement field
    A::Float64 = 2000       # amplitude [m]
    lmin::Int64 = 10        # minimum wavenumber
    lmax::Int64 = 20        # maximum wavenumber
end

"""
$(TYPEDSIGNATURES)
Random initial conditions for the interface displacement η
in the shallow water equations. The flow (u,v) is zero initially.
This kicks off gravity waves that will interact with orography."""
function initialize!(   progn::PrognosticVariables{NF},
                        initial_conditions::RandomWaves,
                        model::ShallowWater) where NF
        
    (;A, lmin, lmax) = initial_conditions
    (;trunc) = progn

    η = progn.surface.timesteps[1].pres
    η .= randn(LowerTriangularMatrix{Complex{NF}},trunc+2,trunc+1)

    # zero out other wavenumbers
    η[1:min(lmin,trunc+2),:] .= 0
    η[min(lmax+2,trunc+2):trunc+2,:] .= 0

    # scale to amplitude
    η_grid = gridded(η,model.spectral_transform)
    η_min,η_max = extrema(η_grid)
    η .*= (A/max(abs(η_min),abs(η_max)))

    return nothing
end