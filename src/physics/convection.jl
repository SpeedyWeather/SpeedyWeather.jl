abstract type AbstractConvection <: AbstractParameterization end

export SimplifiedBettsMiller

"""The simplified Betts-Miller convection scheme from Frierson, 2007,
https://doi.org/10.1175/JAS3935.1. This implements the qref-formulation
in their paper. Fields and options are $(TYPEDFIELDS)"""
@kwdef struct SimplifiedBettsMiller{NF} <: AbstractConvection
    "[OPTION] Relaxation time for profile adjustment"
    time_scale::Second = Hour(4)

    "[OPTION] Relative humidity for reference profile [1]"
    relative_humidity::NF = 0.7
end

Adapt.@adapt_structure SimplifiedBettsMiller

# generator function 
SimplifiedBettsMiller(SG::SpectralGrid; kwargs...) = SimplifiedBettsMiller{SG.NF}(; kwargs...)
initialize!(::SimplifiedBettsMiller, ::PrimitiveEquation) = nothing

# function barrier
function parameterization!(ij, diagn, progn, convection_scheme::SimplifiedBettsMiller, model)
    convection!(ij, diagn, convection_scheme, model)
end

"""
$(TYPEDSIGNATURES)
calculates temperature and humidity tendencies for the convection scheme following the
simplified Betts-Miller convection. Starts with a first-guess relaxation to determine
the convective criteria (none, dry/shallow or deep), then adjusts reference profiles
for thermodynamic consistency (e.g. in dry convection the humidity profile is non-precipitating),
and relaxes current vertical profiles to the adjusted references."""
function convection!(ij, diagn, SBM::SimplifiedBettsMiller, model)

    (; geometry, clausius_clapeyron, planet, atmosphere, time_stepping) = model
    NF = eltype(diagn.grid.temp_grid)
    σ = geometry.σ_levels_full
    σ_half = geometry.σ_levels_half
    Δσ = geometry.σ_levels_thick
    nlayers = length(σ)
    
    # use previous time step for more stable calculations
    temp = diagn.grid.temp_grid_prev
    temp_virt = diagn.grid.temp_virt_grid
    humid = diagn.grid.humid_grid_prev
    geopot = diagn.dynamics.geopot
    temp_tend = diagn.tendencies.temp_tend_grid
    humid_tend = diagn.tendencies.humid_tend_grid
    pₛ = diagn.grid.pres_grid_prev[ij]

    # thermodynamics
    Lᵥ = clausius_clapeyron.latent_heat_condensation  # latent heat of vaporization
    cₚ = clausius_clapeyron.heat_capacity             # heat capacity

    # use scratch arrays for temp_ref_profile, humid_ref_profile
    temp_ref_profile =  diagn.dynamics.a_grid               # temperature [K] reference profile to adjust to
    humid_ref_profile = diagn.dynamics.b_grid               # specific humidity [kg/kg] profile to adjust to
    geopot = diagn.dynamics.uv∇lnp                          # geopotential [m²/s²] on full levels

    # TODO move this to its own parameterization?
    geopotential!(ij, geopot, temp_virt, model.orography.orography, planet.gravity, model.geopotential)

    # CONVECTIVE CRITERIA AND FIRST GUESS RELAXATION
    # Create pseudo column for surface_temp_humid function (this needs to be updated later)
    # For now, use surface values directly
    temp_parcel = temp[ij, nlayers]                         # TODO use surface or skin temperature?
    humid_parcel = humid[ij, nlayers]

    level_zero_buoyancy = pseudo_adiabat!(ij, temp_ref_profile,
                                            temp_parcel, humid_parcel,
                                            temp_virt, geopot, pₛ, σ,
                                            clausius_clapeyron)
            
    for k in level_zero_buoyancy:nlayers
        qsat = saturation_humidity(temp_ref_profile[ij, k], pₛ*σ[k], clausius_clapeyron)
        humid_ref_profile[ij, k] = qsat*SBM.relative_humidity
    end

    local Pq::NF = 0        # precipitation due to drying
    local PT::NF = 0        # precipitation due to cooling
    local ΔT::NF = 0        # vertically uniform temperature profile adjustment
    local Qref::NF = 0      # = ∫_pₛ^p_LZB -humid_ref_profile dp

    # skip constants compared to Frierson 2007, i.e. no /τ, /gravity, *cₚ/Lᵥ
    for k in level_zero_buoyancy:nlayers
        # Frierson's equation (1)
        # δq = -(humid[ij, k] - humid_ref_profile[ij, k])/SBM.time_scale.value
        # Pq -= δq*Δσ[k]/gravity
        #
        # δT = -(temp[ij, k] - temp_ref_profile[ij, k])/SBM.time_scale.value
        # PT += δT*Δσ[k]/gravity*cₚ/Lᵥ

        # shorter form with same sign (τ, gravity, cₚ, Lᵥ all positive) to be reused
        Pq += (humid[ij, k] - humid_ref_profile[ij, k])*Δσ[k]
        PT -= (temp[ij, k] - temp_ref_profile[ij, k])*Δσ[k]
    end

    # ADJUST PROFILES FOLLOWING FRIERSON 2007
    deep_convection = Pq > 0 && PT > 0
    shallow_convection = Pq <= 0 && PT > 0

    # escape immediately for no convection
    no_convection = !(deep_convection || shallow_convection)
    no_convection && return nothing

    # height of zero buoyancy level in σ coordinates
    Δσ_lzb = σ_half[nlayers+1] - σ_half[level_zero_buoyancy]   

    if deep_convection

        ΔT = (PT - Pq*Lᵥ/cₚ)/Δσ_lzb         # eq (5) but reusing PT, Pq, and /cₚ already included

        for k in level_zero_buoyancy:nlayers
            temp_ref_profile[ij, k] -= ΔT   # equation (6)
        end
    
    elseif shallow_convection
        
        # FRIERSON'S QREF SCHEME
        # "changing the reference profiles for both temperature and humidity so the
        # precipitation is zero.

        for k in level_zero_buoyancy:nlayers
            Qref -= humid_ref_profile[ij, k]*Δσ[k]  # eq (11) but in σ coordinates
        end
        fq = 1 - Pq/Qref                    # = 1 - Δq/Qref in eq (12) but we reuse Pq

        ΔT = PT/Δσ_lzb                      # equation (14), reuse PT and in σ coordinates
        for k in level_zero_buoyancy:nlayers
            humid_ref_profile[ij, k] *= fq      # update humidity profile, eq (13)
            temp_ref_profile[ij, k] -= ΔT       # update temperature profile, eq (15)
        end
    end

	(; gravity) = planet
	(; water_density) = atmosphere
	(; Δt_sec) = time_stepping

    # Initialize rain accumulation for this grid point
    rain_convection = zero(NF)

    # GET TENDENCIES FROM ADJUSTED PROFILES
    for k in level_zero_buoyancy:nlayers
        temp_tend[ij, k] -= (temp[ij, k] - temp_ref_profile[ij,k]) / SBM.time_scale.value
        δq = (humid[ij, k] - humid_ref_profile[ij,k]) / SBM.time_scale.value
        humid_tend[ij, k] -= δq

        # convective precipitation (rain), integrate dq\dt [(kg/kg)/s] vertically
        rain = max(δq * Δσ[k], zero(δq))        # only integrate excess humidity for precip (no reevaporation)
        rain_convection += rain                 # integrate vertically, Formula 25, unit [m]
    end
 
    pₛΔt_gρ = (pₛ * Δt_sec / gravity / water_density) * deep_convection # enforce no precip for shallow conv
	rain_convection *= pₛΔt_gρ                                           # convert to [m] of rain during Δt
    
    # Store precipitation in diagnostic arrays
    diagn.physics.rain_convection[ij] += rain_convection

    # TODO: double check, because preivously this was rain_rate_convection, which didn't exist
    diagn.physics.total_precipitation_rate[ij] = rain_convection / Δt_sec    # rate: convert to [m/s] of rain

    # Update cloud top
    diagn.physics.cloud_top[ij] = min(diagn.physics.cloud_top[ij], level_zero_buoyancy)       # clouds reach to top of convection
    return nothing
end

"""
$(TYPEDSIGNATURES)
Calculates the moist pseudo adiabat given temperature and humidity of surface parcel.
Follows the dry adiabat till condensation and then continues on the pseudo moist-adiabat
with immediate condensation to the level of zero buoyancy. Levels above are skipped,
set to NaN instead and should be skipped in the relaxation."""
function pseudo_adiabat!(
    ij,
    temp_ref_profile,
    temp_parcel,
    humid_parcel,
    temp_virt_environment,
    geopot,
    pres,
    σ,
    clausius_clapeyron,
)
    NF = eltype(temp_ref_profile)

    # thermodynamics
    (; R_dry, R_vapour) = clausius_clapeyron
    Lᵥ = clausius_clapeyron.latent_heat_condensation    # latent heat of vaporization
    cₚ = clausius_clapeyron.heat_capacity               # heat capacity

    κ = R_dry/cₚ                         
    ε = clausius_clapeyron.mol_ratio
    μ = (1-ε)/ε                             # for virtual temperature

    nlayers = length(σ)                         # number of vertical levels
    for i in 1:nlayers
        temp_ref_profile[ij, i] = NF(NaN)           # reset profile from any previous calculation, TODO necessary?
    end
    temp_ref_profile[ij, nlayers] = temp_parcel     # start profile at surface with parcel temperature

    local saturated::Bool = false           # did the parcel reach saturation yet?
    local buoyant::Bool = true              # is the parcel still buoyant?
    local k::Int = nlayers                  # layer index top to surface
    local temp_virt_parcel::NF = temp_parcel * (1 + μ*humid_parcel)

    while buoyant && k > 1                  # calculate moist adiabat while buoyant till top
        k -= 1                              # one level up
            
        if !saturated                       # if not saturated yet follow dry adiabat
            # dry adiabatic ascent and saturation humidity of that temperature 
            temp_parcel_dry = temp_parcel*(σ[k]/σ[k+1])^κ
            sat_humid = saturation_humidity(temp_parcel_dry, σ[k]*pres, clausius_clapeyron)
                    
            # set to saturated when the dry adiabatic ascent would reach saturation
            # then follow moist adiabat instead (see below)
            saturated = humid_parcel >= sat_humid
        end
    
        if saturated            
            # calculate moist/pseudo adiabatic lapse rate, dT/dΦ = -Γ/cp
            T, Tᵥ, q = temp_parcel, temp_virt_parcel, humid_parcel  # for brevity
            A = q*Lᵥ / ((1-q)^2 * R_dry)
            B = q*Lᵥ^2 / ((1-q)^2 * cₚ * R_vapour)
            Γ = (1 + A/Tᵥ) / (1 + B/T^2)
                
            ΔΦ = geopot[ij, k] - geopot[ij, k+1]                    # vertical gradient in geopotential
            temp_parcel = temp_parcel - ΔΦ/cₚ*Γ                     # new temperature of parcel at k
                
            # at new (lower) temperature condensation occurs immediately
            # new humidity equals to that saturation humidity
            humid_parcel = saturation_humidity(temp_parcel, σ[k]*pres, clausius_clapeyron)
        else
            temp_parcel = temp_parcel_dry       # else parcel temperature following dry adiabat
        end
    
        # use dry/moist adiabatic ascent for reference profile
        temp_ref_profile[ij, k] = temp_parcel

        # check whether parcel is still buoyant wrt to environment
        # use virtual temperature as it's equivalent to density
        temp_virt_parcel = temp_parcel*(1 + μ*humid_parcel)         # virtual temperature of parcel
        buoyant = temp_virt_parcel > temp_virt_environment[ij, k]     
    end
    
    # if parcel isn't buoyant anymore set last temperature (with negative buoyancy) back to NaN
    temp_ref_profile[ij, k] = !buoyant ? NF(NaN) : temp_ref_profile[ij, k]
    
    # level of zero buoyancy is reached when the loop stops, but in case it's at the top it's still buoyant
    level_zero_buoyancy = k + (1-buoyant)
    return level_zero_buoyancy
end

export DryBettsMiller

"""
The simplified Betts-Miller convection scheme from Frierson, 2007,
https://doi.org/10.1175/JAS3935.1 but with humidity set to zero.
Fields and options are
$(TYPEDFIELDS)"""
@kwdef struct DryBettsMiller{NF} <: AbstractConvection
    "[OPTION] Relaxation time for profile adjustment"
    time_scale::Second = Hour(4)
end

function Adapt.adapt_structure(to, DBM::DryBettsMiller{NF}) where NF
    return DryBettsMiller{NF}(adapt_structure(to, DBM.time_scale))
end

# generator function
DryBettsMiller(SG::SpectralGrid; kwargs...) = DryBettsMiller{SG.NF}(; kwargs...)
initialize!(::DryBettsMiller, ::PrimitiveEquation) = nothing

# function barrier
function parameterization!(ij, diagn, progn, convection_scheme::DryBettsMiller, model)
    convection!(ij, diagn, convection_scheme, model)
end

"""
$(TYPEDSIGNATURES)
calculates temperature tendency for the dry convection scheme following the
simplified Betts-Miller convection from Frierson 2007 but with zero humidity.
Starts with a first-guess relaxation to determine the convective criterion,
then adjusts the reference profiles
for thermodynamic consistency (e.g. in dry convection the humidity profile is non-precipitating),
and relaxes current vertical profiles to the adjusted references."""
function convection!(ij, diagn, DBM::DryBettsMiller, model)

    (; geometry, atmosphere) = model
    NF = eltype(diagn.grid.temp_grid)
    σ = geometry.σ_levels_full
    σ_half = geometry.σ_levels_half
    Δσ = geometry.σ_levels_thick
    nlayers = length(σ)

    # use previous time step for more stable calculations
    temp = diagn.grid.temp_grid_prev
    temp_tend = diagn.tendencies.temp_tend_grid

    # use work arrays for temp_ref_profile
    temp_ref_profile = diagn.dynamics.a_grid     # temperature [K] reference profile to adjust to

    # CONVECTIVE CRITERIA AND FIRST GUESS RELAXATION
    # Use surface temperature directly (simplified for now)
    temp_parcel = temp[ij, nlayers]
    level_zero_buoyancy = dry_adiabat!(ij, temp_ref_profile,
                                        temp,
                                        temp_parcel,
                                        σ,
                                        atmosphere)

    local PT::NF = 0        # precipitation due to cooling
    local ΔT::NF = 0        # vertically uniform temperature profile adjustment

    # skip constants compared to Frierson 2007, i.e. no /τ, /gravity, *cₚ/Lᵥ
    for k in level_zero_buoyancy:nlayers
        # Frierson's equation (1)
        # δT = -(temp[ij, k] - temp_ref_profile[ij, k])/DBM.time_scale.value
        # PT += δT*Δσ[k]/gravity*cₚ/Lᵥ

        # shorter form with same sign (τ, gravity, cₚ, Lᵥ all positive) to be reused
        PT -= (temp[ij, k] - temp_ref_profile[ij, k])*Δσ[k]
    end

    # ADJUST PROFILES FOLLOWING FRIERSON 2007
    convection = PT > 0
    convection || return nothing            # escape immediately for no convection

    # height of zero buoyancy level in σ coordinates
    Δσ_lzb = σ_half[nlayers+1] - σ_half[level_zero_buoyancy]   
    ΔT = PT/Δσ_lzb                          # eq (5) or (14) but reusing PT
    for k in level_zero_buoyancy:nlayers
        temp_ref_profile[ij, k] -= ΔT           # equation (6) or equation (15)
        temp_tend[ij, k] -= (temp[ij, k] - temp_ref_profile[ij, k]) / DBM.time_scale.value
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
Calculates the moist pseudo adiabat given temperature and humidity of surface parcel.
Follows the dry adiabat till condensation and then continues on the pseudo moist-adiabat
with immediate condensation to the level of zero buoyancy. Levels above are skipped,
set to NaN instead and should be skipped in the relaxation."""
function dry_adiabat!(
    ij,
    temp_ref_profile,
    temp_environment,
    temp_parcel,
    σ,
    atmosphere,
)
    κ = atmosphere.heat_capacity

    @boundscheck length(temp_ref_profile) ==
        length(σ) == length(temp_environment) || throw(BoundsError)

    nlayers = length(temp_ref_profile)      # number of vertical levels

    for i in 1:nlayers
        temp_ref_profile[ij, i] = NF(NaN)          # reset profile from any previous calculation
    end
    temp_ref_profile[ij, nlayers] = temp_parcel    # start profile at surface with parcel temperature

    local buoyant::Bool = true              # is the parcel still buoyant?
    local k::Int = nlayers                  # layer index top to surface

    while buoyant && k > 1                  # calculate moist adiabat while buoyant till top
        k -= 1                              # one level up
            
        # dry adiabatic ascent
        temp_parcel = temp_parcel*(σ[k]/σ[k+1])^κ
        temp_ref_profile[ij, k] = temp_parcel

        # check whether parcel is still buoyant wrt to environment
        buoyant = temp_parcel > temp_environment[ij, k] 
    end
    
    # if parcel isn't buoyant anymore set last temperature (with negative buoyancy) back to NaN
    temp_ref_profile[ij, k] = !buoyant ? NF(NaN) : temp_ref_profile[ij, k]    
    
    # level of zero buoyancy is reached when the loop stops, but in case it's at the top it's still buoyant
    level_zero_buoyancy = k + (1-buoyant)
    return level_zero_buoyancy
end

export ConvectiveHeating

"""Convective heating as defined by Lee and Kim, 2003, JAS
implemented as convection parameterization. Fields are
$(TYPEDFIELDS)"""
@kwdef struct ConvectiveHeating{NF, VectorType} <: AbstractConvection
    "[OPTION] Q_max heating strength as 1K/time_scale"
    time_scale::Second = Hour(12)

    "[OPTION] Pressure of maximum heating [hPa]"
    p₀::NF = 525

    "[OPTION] Vertical extent of heating [hPa]"
    σₚ::NF = 200

    "[OPTION] Latitude of heating [˚N]"
    θ₀::NF = 0

    "[OPTION] Latitudinal width of heating [˚]"
    σθ::NF = 20

    "[DERIVED] Latitudinal mask"
    lat_mask::VectorType
end

Adapt.@adapt_structure ConvectiveHeating

# generator
ConvectiveHeating(SG::SpectralGrid; kwargs...) = ConvectiveHeating{SG.NF, SG.VectorType}(lat_mask=zeros(SG.nlat); kwargs...)

# precompute latitudinal mask
function initialize!(C::ConvectiveHeating, model::PrimitiveEquation)
    
    θ = model.geometry.latd
    (; θ₀, σθ) = C
    
    # Lee and Kim, 2003, eq. 2
    @. C.lat_mask .= cosd((θ-θ₀)/σθ)^2
end

# function barrier
parameterization!(ij, diagn, progn, convection_scheme::ConvectiveHeating, model) = convection!(ij, diagn, convection_scheme, model)

@inline function convection!(
    ij,
    diagn::DiagnosticVariables,
    scheme::ConvectiveHeating,
    model::PrimitiveEquation,
)
    # Get latitude ring index and latitude
    j = whichring(diagn.grid.temp_grid.grid, ij)
    latd = model.geometry.latd[j]
    σ = model.geometry.σ_levels_full
    
    # escape immediately if not in the tropics
    abs(latd) >= scheme.σθ && return nothing

    p₀ = scheme.p₀*100      # hPa -> Pa
    σₚ = scheme.σₚ*100      # hPa -> Pa
    cos²θ_term = scheme.lat_mask[j]
    Qmax = 1/Second(scheme.time_scale).value

    # get pressure levels for this grid point
    pres = diagn.grid.pres_grid_prev
    temp_tend = diagn.tendencies.temp_tend_grid
    nlayers = size(temp_tend, 2)
    
    @inbounds for k in 1:nlayers
        p = pres[ij] * σ[k]      # Pressure in Pa

        # Lee and Kim, 2003, eq. 2
        temp_tend[ij, k] += Qmax*exp(-((p-p₀)/σₚ)^2 / 2)*cos²θ_term
    end
    return nothing
end

function variables(::AbstractConvection)
    return (
        DiagnosticVariable(name=:rain_convection, dims=Grid2D(), desc="Convective precipitation", units="m"),
        DiagnosticVariable(name=:cloud_top, dims=Grid2D(), desc="Cloud top level", units="1"),
        DiagnosticVariable(name=:total_precipitation_rate, dims=Grid2D(), desc="Total precipitation rate", units="m/s"),
    )
end