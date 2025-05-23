abstract type AbstractInitialConditions <: AbstractModelComponent end

export InitialConditions
@kwdef struct InitialConditions{V, P, T, H} <: AbstractInitialConditions
    vordiv::V = ZeroInitially()
    pres::P = ZeroInitially()
    temp::T = ZeroInitially()
    humid::H = ZeroInitially()
end

function initialize!(
    progn::PrognosticVariables,
    IC::InitialConditions,
    model::AbstractModel
)
    has(model, :vor)   && initialize!(progn, IC.vordiv, model)
    has(model, :pres)  && initialize!(progn, IC.pres,   model)
    has(model, :temp)  && initialize!(progn, IC.temp,   model)
    has(model, :humid) && initialize!(progn, IC.humid,  model)
end

InitialConditions(::Type{<:Barotropic}) = InitialConditions(; vordiv = RandomVelocity())
InitialConditions(::Type{<:ShallowWater}) = InitialConditions(; vordiv = ZonalJet())
function InitialConditions(::Type{<:PrimitiveDry})
    vordiv = ZonalWind()
    pres = PressureOnOrography()
    temp = JablonowskiTemperature()
    return InitialConditions(;vordiv, pres, temp)
end

function InitialConditions(::Type{<:PrimitiveWet})
    vordiv = ZonalWind()
    pres = PressureOnOrography()
    temp = JablonowskiTemperature()
    humid = ConstantRelativeHumidity()
    return InitialConditions(;vordiv, pres, temp, humid)
end

export ZeroInitially
struct ZeroInitially <: AbstractInitialConditions end
initialize!(::PrognosticVariables,::ZeroInitially,::AbstractModel) = nothing

# to avoid a breaking change, like ZeroInitially
export StartFromRest
@kwdef struct StartFromRest{P, T, H} <: AbstractInitialConditions
    pres::P = ConstantPressure()
    temp::T = JablonowskiTemperature()
    humid::H = ZeroInitially()
end

initialize!(::PrognosticVariables, ::StartFromRest, ::Barotropic) = nothing

function initialize!(
    progn::PrognosticVariables,
    IC::StartFromRest,
    model::AbstractModel,
)
    has(model, :pres)  && initialize!(progn, IC.pres, model)
    has(model, :temp)  && initialize!(progn, IC.temp, model)
    has(model, :humid) && initialize!(progn, IC.humid, model)
end

export RandomVorticity

"""Start with random vorticity as initial conditions
$(TYPEDFIELDS)"""
@kwdef mutable struct RandomVorticity <: AbstractInitialConditions
    "[OPTION] Power of the spectral distribution k^power"
    power::Float64 = -3

    "[OPTION] the (approximate) amplitude in [1/s], used as standard deviation of spherical harmonic coefficients"
    amplitude::Float64 = 1e-4

    "[OPTION] Maximum wavenumber"
    max_wavenumber::Int = 20

    "[OPTION] Random number generator seed, 0=randomly seed from Julia's GLOBAL_RNG"
    seed::Int = 123

    "Independent random number generator for this random process"
    random_number_generator::Random.Xoshiro = Random.Xoshiro(seed)
end

"""$(TYPEDSIGNATURES)
Start with random vorticity as initial conditions"""
function initialize!(   progn::PrognosticVariables{NF},
                        initial_conditions::RandomVorticity,
                        model::Barotropic) where NF

    # reseed the random number generator, for seed=0 randomly seed from Julia's global RNG
    seed = initial_conditions.seed == 0 ? rand(UInt) : initial_conditions.seed
    RNG = initial_conditions.random_number_generator
    Random.seed!(RNG, seed)

    lmax = progn.trunc + 1
    power = initial_conditions.power + 1    # +1 as power is summed of orders m

    A = initial_conditions.amplitude
    ξ = zeros(LowerTriangularArray{Complex{NF}}, lmax+1, lmax, model.spectral_grid.nlayers)

    lm = 0
    for m in 1:lmax
        for l in m:lmax
            lm += 1
            ξ[lm, :] .= A*l^power*(2rand(RNG, Complex{NF}) - (1 + 1im))
        end
    end

    spectral_truncation!(ξ, initial_conditions.max_wavenumber)
    ξ[1:lmax, :] .= 0       # remove zonal modes and mean
    set!(progn, model; vor=ξ, lf=1)
end

export RandomVelocity

"""Start with random velocity as initial conditions
$(TYPEDFIELDS)"""
@kwdef mutable struct RandomVelocity <: AbstractInitialConditions
    "[OPTION] maximum speed [ms⁻¹]"
    max_speed::Float64 = 60

    "[OPTION] Maximum wavenumber after truncation"
    truncation::Int = 15

    "[OPTION] Random number generator seed, 0=randomly seed from Julia's GLOBAL_RNG"
    seed::Int = 0

    "Independent random number generator for this random process"
    random_number_generator::Random.Xoshiro = Random.Xoshiro(seed)
end

"""$(TYPEDSIGNATURES)
Start with random vorticity as initial conditions"""
function initialize!(   progn::PrognosticVariables{NF},
                        initial_conditions::RandomVelocity,
                        model::Barotropic) where NF

    # reseed the random number generator, for seed=0 randomly seed from Julia's global RNG
    seed = initial_conditions.seed == 0 ? rand(UInt) : initial_conditions.seed
    RNG = initial_conditions.random_number_generator
    Random.seed!(RNG, seed)

    (; GridVariable2D, nlat_half, radius, trunc, nlayers) = model.spectral_grid
    A = convert(NF, initial_conditions.max_speed) * 2
    u = 2A*rand(GridVariable2D, nlat_half) .- A
    v = 2A*rand(GridVariable2D, nlat_half) .- A

    u_spectral = transform(u, model.spectral_transform)
    v_spectral = transform(v, model.spectral_transform)

    spectral_truncation!(u_spectral, initial_conditions.truncation)
    spectral_truncation!(v_spectral, initial_conditions.truncation)

    ξ = curl(u_spectral, v_spectral, model.spectral_transform; radius)
    ξ[1] = 0        # remove mean

    # repeat over vertical layers
    ξks = repeat(ξ, 1, nlayers)
    set!(progn, model; vor=ξks, lf=1)

    return nothing
end

export ZonalJet

"""A struct that contains all parameters for the Galewsky et al, 2004 zonal jet
intitial conditions for the ShallowWaterModel. Default values as in Galewsky.
$(TYPEDFIELDS)"""
@kwdef mutable struct ZonalJet <: AbstractInitialConditions
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

"""
$(TYPEDSIGNATURES)
Initial conditions from Galewsky, 2004, Tellus"""
function initialize!(   progn::PrognosticVariables,
                        initial_conditions::ZonalJet,
                        model::AbstractModel)

    (; latitude, width, umax) = initial_conditions               # for jet
    (; perturb_lat, perturb_lon, perturb_xwidth,                 # for perturbation
        perturb_ywidth, perturb_height) = initial_conditions

    θ₀ = (latitude-width)/360*2π    # southern boundary of jet [radians]
    θ₁ = (latitude+width)/360*2π    # northern boundary of jet
    eₙ = exp(-4/(θ₁-θ₀)^2)          # normalisation

    θ₂ = perturb_lat*2π/360         # perturbation latitude [radians]
    α = perturb_xwidth*2π/360       # zonal extent of interface perturbation [radians]
    β = perturb_ywidth*2π/360       # meridional extent of interface perturbation [radians]
    λ = perturb_lon*2π/360          # perturbation longitude [radians]

    (; rotation, gravity) = model.planet
    (; Grid, NF, nlat_half) = model.spectral_grid
    (; coslat⁻¹, radius) = model.geometry

    u_grid = zeros(Grid{NF}, nlat_half, 1)
    η_perturb_grid = zeros(Grid{NF}, nlat_half)
    lat = RingGrids.get_lat(Grid, nlat_half)
    lons, _ = RingGrids.get_lonlats(Grid, nlat_half)

    for (j, ring) in enumerate(eachring(u_grid))
        θ = lat[j]             # latitude in radians

        # velocity per latitude
        if θ₀ < θ < θ₁
            u_θ = umax/eₙ*exp(1/(θ-θ₀)/(θ-θ₁))  # u as in Galewsky, 2004
        else
            u_θ = 0
        end

        # lon-constant part of perturbation
        ηθ = perturb_height*cos(θ)*exp(-((θ₂-θ)/β)^2)

        # store in all longitudes
        for ij in ring
            u_grid[ij] = u_θ/radius*coslat⁻¹[j]   # include scaling for curl!

            # calculate perturbation (possibly shifted in lon compared to Galewsky 2004)
            ϕ = lons[ij] - λ
            η_perturb_grid[ij] = exp(-(ϕ/α)^2)*ηθ
        end
    end

    # the following obtain initial conditions for η from u, v=0 via
    # 0 = -∇⋅((ζ+f)*(-v, u)) - ∇²((u^2 + v^2)/2), i.e.
    # invert the Laplacian for
    # ∇²(gη) = -∇⋅(0, (ζ+f)*u) - ∇²(u^2/2)
    u = transform(u_grid, model.spectral_transform)

    # get vorticity initial conditions from curl of u, v
    v = zero(u)     # meridional velocity zero for these initial conditions
    vor = progn.vor[1]
    curl!(vor, u, v, model.spectral_transform)

    # compute the div = -∇⋅(0,(ζ+f)*u) = ∇×((ζ+f)*u, 0) term, v=0
    vor_grid = transform(vor, model.spectral_transform)
    f = coriolis(vor_grid; rotation)

    # includes 1/coslat/radius from above for curl!
    # but *radius^2 for the ∇⁻²! operation below! 
    vor_flux_grid = @. (vor_grid + f) * u_grid * radius^2
    vor_flux = transform(vor_flux_grid, model.spectral_transform)
    div = curl(vor_flux, v, model.spectral_transform)

    # compute the -∇²(u^2/2) term, add to div, divide by gravity
    RingGrids.scale_coslat!(u_grid)     # remove coslat scaling
    u_grid .*= radius                   # no radius scaling as we'll apply ∇⁻²(∇²) (would cancel)
    @. u_grid = convert(NF,1/2) * u_grid^2
    u²_half = transform!(u, u_grid, model.spectral_transform)
    ∇²!(div, u²_half, model.spectral_transform, flipsign=true, add=true)
    div .*= inv(gravity)

    # invert Laplacian to obtain η
    pres = progn.pres[1]
    ∇⁻²!(pres, div, model.spectral_transform)

    # add perturbation (reuse u array)
    η_perturb = transform(η_perturb_grid, model.spectral_transform)
    pres .+= η_perturb
    spectral_truncation!(pres)
    return nothing
end

export ZonalWind

"""
$(TYPEDSIGNATURES)
Create a struct that contains all parameters for the Jablonowski and Williamson, 2006
intitial conditions for the primitive equation model. Default values as in Jablonowski.
$(TYPEDFIELDS)"""
@kwdef struct ZonalWind <: AbstractInitialConditions
    "conversion from σ to Jablonowski's ηᵥ-coordinates"
    η₀::Float64 = 0.252

    "[OPTION] max amplitude of zonal wind [m/s]"
    u₀::Float64 = 35

    # PERTURBATION
    "[OPTION] perturbation centred at [˚N]"
    perturb_lat::Float64 = 40

    "[OPTION] perturbation centred at [˚E]"
    perturb_lon::Float64 = 20

    "[OPTION] perturbation strength [m/s]"
    perturb_uₚ::Float64 = 1

    "[OPTION] radius of Gaussian perturbation in units of Earth's radius [1]"
    perturb_radius::Float64 = 1/10
end

"""
$(TYPEDSIGNATURES)
Initial conditions from Jablonowski and Williamson, 2006, QJR Meteorol. Soc"""
function initialize!(   progn::PrognosticVariables{NF},
                        initial_conditions::ZonalWind,
                        model::PrimitiveEquation) where NF

    (; u₀, η₀) = initial_conditions
    (; perturb_uₚ, perturb_radius) = initial_conditions
    λc = initial_conditions.perturb_lon
    sinφc, cosφc = sind(initial_conditions.perturb_lat), cosd(initial_conditions.perturb_lat)
    (; radius) = model.spectral_grid
    R = radius*perturb_radius           # spatial extent of perturbation

    function vor(λ, φ, η)               # longitude, latitude (degree), sigma coordinate

        # great circle distance to perturbation
        X = clamp(sinφc*sind(φ) + cosφc*cosd(φ)*cosd(λ-λc), 0, 1)
        r = radius*acos(X)
        
        # Eq (3), the unperturbed zonal wind
        ζ = -4u₀/radius*cos((η - η₀)*π/2)^(3/2)*sind(φ)*cosd(φ)*(2 - 5sind(φ)^2)

        # Eq (12), the perturbation
        perturbation = perturb_uₚ/radius * exp(-(r/R)^2) *
            (tand(φ) - 2*(radius/R)^2*acos(X)*(sinφc*cosd(φ) - cosφc*sind(φ)*cosd(λ-λc))/sqrt(1-X^2))

        return ζ+perturbation
    end

    function div(λ, φ, η)               # longitude, latitude (degree), sigma coordinate

        # great circle distance to perturbation
        X = clamp(sinφc*sind(φ) + cosφc*cosd(φ)*cosd(λ-λc), 0, 1)
        r = radius*acos(X)

        # Eq (13)
        return -2perturb_uₚ*radius/R^2 * exp(-(r/R)^2) * acos(X)/sqrt(1-X^2) * cosφc * sind(λ-λc)
    end

    # apply those to set the initial conditions for vor, div
    set!(progn, model; vor=vor, div=div, lf=1)
    return nothing
end

export RossbyHaurwitzWave

"""Rossby-Haurwitz wave initial conditions as in Williamson et al. 1992, J Computational Physics
with an additional cut-off amplitude `c` to filter out tiny harmonics in the vorticity field.
Parameters are $(TYPEDFIELDS)"""
@kwdef struct RossbyHaurwitzWave <: AbstractInitialConditions
    m::Int = 4
    ω::Float64 = 7.848e-6
    K::Float64 = 7.848e-6
    c::Float64 = 1e-10
end

"""$(TYPEDSIGNATURES)
Rossby-Haurwitz wave initial conditions as in Williamson et al. 1992, J Computational Physics
with an additional cut-off amplitude `c` to filter out tiny harmonics in the vorticity field."""
function initialize!(
    progn::PrognosticVariables,
    initial_conditions::RossbyHaurwitzWave,
    model::AbstractModel,
)
    (; m, ω, K, c) = initial_conditions
    (; geometry) = model
    Ω = model.planet.rotation
    R = model.spectral_grid.radius
    g = model.planet.gravity

    # Rossby-Haurwitz wave defined through vorticity ζ as a function of
    # longitude λ, latitude θ (in degrees), sigma level σ (vertically constant though)
    # see Williamson et al. 1992, J Computational Physics, eq 145
    ζ(λ, θ, σ) = 2ω*sind(θ) - K*sind(θ)*cosd(θ)^m*(m^2 + 3m + 2)*cosd(m*λ)
    # see Williamson et al. 1992, J Computational Physics, eq 147 - 149, 146
    A(λ, θ) = ω/2 * (2Ω+ω)*cosd(θ)^2 + K^2/4*cosd(θ)^(2m)*((m+1)*cosd(θ)^2 + (2m^2-m-2) - 2m^2/(cosd(θ)^2))
    B(λ, θ) = (2(Ω+ω)*K)/((m+1)*(m+2))*cosd(θ)^m*( (m^2+2m+2) - (m+1)^2*cosd(θ)^2 )
    C(λ, θ) = K^2/4*cosd(θ)^(2m)*((m+1) * cosd(θ)^2 - (m + 2))

    η(λ, θ) = R^2/g*(A(λ,θ) + B(λ,θ)*cosd(m*λ) + C(λ,θ)*cosd(2m*λ))

    set!(progn, geometry, vor = ζ)
    model isa ShallowWater && set!(progn, geometry, pres = η)
    set!(progn, geometry, div = 0)  # technically not needed, but set to zero for completeness

    # filter low values below cutoff amplitude c
    vor = progn.vor[1]  # 1 = first leapfrog timestep
    low_values = abs.(vor) .< c
    vor[low_values] .= 0
    if model isa ShallowWater
        pres = progn.pres[1]
        low_value = abs.(pres) .< c
        pres[low_value] .= 0
    end

    return nothing
end

export JablonowskiTemperature

"""
$(TYPEDSIGNATURES)
Create a struct that contains all parameters for the Jablonowski and Williamson, 2006
intitial conditions for the primitive equation model. Default values as in Jablonowski.
$(TYPEDFIELDS)"""
@kwdef struct JablonowskiTemperature <: AbstractInitialConditions
    "conversion from σ to Jablonowski's ηᵥ-coordinates"
    η₀::Float64 = 0.252

    "Sigma coordinates of the tropopause [1]"
    σ_tropopause::Float64 = 0.2

    "max amplitude of zonal wind [m/s]"
    u₀::Float64 = 35

    "temperature difference used for stratospheric lapse rate [K], Jablonowski uses ΔT = 4.8e5 [K]"
    ΔT::Float64 = 0

    "minimum temperature [K] of profile"
    Tmin::Float64 = 200
end

"""
$(TYPEDSIGNATURES)
Initial conditions from Jablonowski and Williamson, 2006, QJR Meteorol. Soc"""
function initialize!(   progn::PrognosticVariables{NF},
                        initial_conditions::JablonowskiTemperature,
                        model::PrimitiveEquation) where NF

    (;u₀, η₀, ΔT, Tmin) = initial_conditions
    (;σ_tropopause) = initial_conditions
    lapse_rate = model.atmosphere.moist_lapse_rate
    (;temp_ref, R_dry) = model.atmosphere
    (;radius, Grid, nlat_half, nlayers) = model.spectral_grid
    (;rotation, gravity) = model.planet
    (;σ_levels_full) = model.geometry

    φ, λ = model.geometry.latds, model.geometry.londs
    S = model.spectral_transform

    # vertical profile
    Tη = zero(σ_levels_full)
    for k in 1:nlayers
        σ = σ_levels_full[k]
        Tη[k] = temp_ref*σ^(R_dry*lapse_rate/gravity)   # Jablonowski and Williamson eq. 4

        if σ < σ_tropopause
            Tη[k] += ΔT*(σ_tropopause-σ)^5              # Jablonowski and Williamson eq. 5
        end
    end

    Tη .= max.(Tη, Tmin)
    temp_grid = zeros(Grid{NF}, nlat_half, nlayers)     # temperature
    aΩ = radius*rotation

    for k in eachlayer(temp_grid)
        η = σ_levels_full[k]    # Jablonowski and Williamson use η for σ coordinates
        ηᵥ = (η - η₀)*π/2       # auxiliary variable for vertical coordinate

        # amplitudes with height
        A1 = 3/4*η*π*u₀/R_dry*sin(ηᵥ)*sqrt(cos(ηᵥ))
        A2 = 2u₀*cos(ηᵥ)^(3/2)

        for (ij, φij) in enumerate(φ)
            sinφ = sind(φij)
            cosφ = cosd(φij)

             # Jablonowski and Williamson, eq. (6)
            temp_grid[ij, k] = Tη[k] + A1*((-2sinφ^6*(cosφ^2 + 1/3) + 10/63)*A2 + (8/5*cosφ^3*(sinφ^2 + 2/3) - π/4)*aΩ)
        end
    end

    set!(progn, model; temp=temp_grid, lf=1)

    return nothing
end

export StartFromFile

"""
Restart from a previous SpeedyWeather.jl simulation via the restart file restart.jld2
Applies interpolation in the horizontal but not in the vertical. restart.jld2 is
identified by
$(TYPEDFIELDS)"""
@kwdef struct StartFromFile <: AbstractInitialConditions
    "path for restart file"
    path::String = pwd()

    "`run_id` of restart file in `run_????/restart.jld2`"
    id::Union{String, Int} = 1
end

"""
$(TYPEDSIGNATURES)
Restart from a previous SpeedyWeather.jl simulation via the restart file restart.jld2
Applies interpolation in the horizontal but not in the vertical."""
function initialize!(   progn_new::PrognosticVariables,
                        initial_conditions::StartFromFile,
                        model::AbstractModel)

    (; path, id ) = initial_conditions

    restart_file = jldopen(joinpath(path, string("run_", run_id_to_string(id)), "restart.jld2"))
    progn_old = restart_file["prognostic_variables"]
    version = restart_file["version"]
    if version != pkgversion(SpeedyWeather)
        @info "Restart file created with SpeedyWeather $version loaded"*
                "but currently used is $(pkgversion(SpeedyWeather))"
    end
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
    (; temp_ref, lapse_rate, R_dry) = model.atmosphere
    (; gravity) = model.planet
    (; nlayers, σ_levels_full) = model.geometry
    (; norm_sphere) = model.spectral_transform # normalization of the l=m=0 spherical harmonic

    # Lapse rate scaled by gravity [K/m / (m²/s²)]
    Γg⁻¹ = lapse_rate/gravity

    # SURFACE TEMPERATURE (store in k = nlayers, but it's actually surface, i.e. k=nlayers+1/2)
    # overwrite with lowermost layer further down
    temp_surf = progn.temp[1][:,end]     # spectral temperature at k=nlayers+1/2

    temp_surf[1] = norm_sphere*temp_ref                 # set global mean surface temperature
    for lm in eachharmonic(geopot_surf, temp_surf)
        temp_surf[lm] -= Γg⁻¹*geopot_surf[lm]           # lower temperature for higher mountains
    end

    # Use lapserate and vertical coordinate σ for profile
    for k in 1:nlayers                                     # k=nlayers overwrites the surface temperature
                                                        # with lowermost layer temperature
        temp = progn.layers[k].timesteps[1].temp
        σₖᴿ = σ_levels_full[k]^(R_dry*Γg⁻¹)             # from hydrostatic equation

        for lm in eachharmonic(temp, temp_surf)
            temp[lm] = temp_surf[lm]*σₖᴿ
        end
    end
end

export PressureOnOrography
struct PressureOnOrography <: AbstractInitialConditions end

"""
$(TYPEDSIGNATURES)
Initialize surface pressure on orography by integrating the
hydrostatic equation with the reference temperature lapse rate."""
function initialize!(   progn::PrognosticVariables,
                        ::PressureOnOrography,
                        model::PrimitiveEquation)

    # temp_ref:     Reference absolute T [K] at surface z = 0
    # lapse_rate:   Reference temperature lapse rate -dT/dz [K/m]
    # gravity:      Gravitational acceleration [m/s^2]
    # R_dry:        Specific gas constant for dry air [J/kg/K]
    # pres_ref:     Reference surface pressure [hPa]

    (;atmosphere) = model
    lapse_rate = model isa PrimitiveDry ? atmosphere.dry_lapse_rate : atmosphere.moist_lapse_rate

    (; temp_ref, pres_ref, R_dry ) = model.atmosphere
    (; gravity ) = model.planet
    (; orography ) = model.orography # orography on the grid

    lnp₀ = log(pres_ref)            # logarithm of reference surface pressure [log(Pa)]
    lnp_grid = zero(orography)      # allocate log surface pressure on grid

    RΓg⁻¹ = R_dry*lapse_rate/gravity         # for convenience
    ΓT⁻¹ = lapse_rate/temp_ref

    for ij in eachgridpoint(lnp_grid, orography)
        lnp_grid[ij] = lnp₀ + log(1 - ΓT⁻¹*orography[ij])/RΓg⁻¹
    end

    set!(progn, model; pres=lnp_grid, lf=1)

    return nothing
end

export ConstantPressure
struct ConstantPressure <: AbstractInitialConditions end

function initialize!(   progn::PrognosticVariables,
                        ::ConstantPressure,
                        model::PrimitiveEquation)

    # logarithm of reference surface pressure [log(Pa)]
    set!(progn, model; pres=log(model.atmosphere.pres_ref))
    return nothing
end

# for shallow water constant pressure = 0 as pres=interface displacement here
initialize!(::PrognosticVariables, ::ConstantPressure, ::ShallowWater) = nothing

export ConstantRelativeHumidity
@kwdef struct ConstantRelativeHumidity <: AbstractInitialConditions
    relhumid_ref::Float64 = 0.7
end

function initialize!(
    progn::PrognosticVariables,
    IC::ConstantRelativeHumidity,
    model::PrimitiveEquation,
)
    (; relhumid_ref) = IC
    (; nlayers, σ_levels_full) = model.geometry

    # get pressure [Pa] on grid
    lnpₛ = progn.pres[1]
    pres_grid = transform(lnpₛ, model.spectral_transform)
    pres_grid .= exp.(pres_grid)

    temp_grid = transform(progn.temp[1], model.spectral_transform)
    humid_grid = zero(temp_grid)

    for k in eachlayer(temp_grid, humid_grid)
        for ij in eachgridpoint(humid_grid)
            pₖ = σ_levels_full[k] * pres_grid[ij]
            q_sat = saturation_humidity(temp_grid[ij, k], pₖ, model.clausius_clapeyron)
            humid_grid[ij, k] = relhumid_ref*q_sat
        end
    end

    set!(progn, model; humid=humid_grid, lf=1)

    return nothing
end

export RandomWaves

"""Parameters for random initial conditions for the interface displacement η
in the shallow water equations.
$(TYPEDFIELDS)"""
@kwdef struct RandomWaves <: AbstractInitialConditions
    # random interface displacement field
    A::Float64 = 2000       # amplitude [m]
    lmin::Int64 = 10        # minimum wavenumber
    lmax::Int64 = 30        # maximum wavenumber
end

RandomWaves(S::SpectralGrid; kwargs...) = RandomWaves(;kwargs...)

"""
$(TYPEDSIGNATURES)
Random initial conditions for the interface displacement η
in the shallow water equations. The flow (u, v) is zero initially.
This kicks off gravity waves that will interact with orography."""
function initialize!(   progn::PrognosticVariables{NF},
                        initial_conditions::RandomWaves,
                        model::ShallowWater) where NF

    (; A, lmin, lmax) = initial_conditions
    (; trunc) = progn

    # start with matrix to have matrix indexing
    ηm = randn(Complex{NF}, trunc+2, trunc+1)

    # zero out other wavenumbers
    ηm[1:min(lmin, trunc+2), :] .= 0
    ηm[min(lmax+2, trunc+2):trunc+2, :] .= 0

    # convert to LowerTriangularMatrix with vector indexing
    η = LowerTriangularArray(ηm)

    # scale to amplitude
    η_grid = transform(η, model.spectral_transform)
    η_min, η_max = extrema(η_grid)
    η .*= (A/max(abs(η_min), abs(η_max)))

    set!(progn, model; pres=η)

    return nothing
end
