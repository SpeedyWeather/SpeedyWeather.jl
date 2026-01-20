abstract type AbstractConvection <: AbstractParameterization end

export BettsMillerConvection

"""The simplified Betts-Miller convection scheme from Frierson, 2007,
https://doi.org/10.1175/JAS3935.1. This implements the qref-formulation
in their paper. Fields and options are $(TYPEDFIELDS)"""
@kwdef struct BettsMillerConvection{NF} <: AbstractConvection
    "[OPTION] Relaxation time for profile adjustment"
    time_scale::Second = Hour(4)

    "[OPTION] Relative humidity for reference profile [1]"
    relative_humidity::NF = 0.7
end

Adapt.@adapt_structure BettsMillerConvection

# generator function 
BettsMillerConvection(SG::SpectralGrid; kwargs...) = BettsMillerConvection{SG.NF}(; kwargs...)
initialize!(::BettsMillerConvection, ::PrimitiveEquation) = nothing

# function barrier
@propagate_inbounds parameterization!(ij, diagn, progn, convection_scheme::BettsMillerConvection, model) =
    convection!(ij, diagn, convection_scheme, model)

"""
$(TYPEDSIGNATURES)
calculates temperature and humidity tendencies for the convection scheme following the
simplified Betts-Miller convection. Starts with a first-guess relaxation to determine
the convective criteria (none, dry/shallow or deep), then adjusts reference profiles
for thermodynamic consistency (e.g. in dry convection the humidity profile is non-precipitating),
and relaxes current vertical profiles to the adjusted references."""
@propagate_inbounds function convection!(ij, diagn, convection::BettsMillerConvection, model)

    (; geometry, clausius_clapeyron, planet, atmosphere, time_stepping) = model
    σ = geometry.σ_levels_full
    σ_half = geometry.σ_levels_half
    Δσ = geometry.σ_levels_thick
    nlayers = length(σ)
    Δt = time_stepping.Δt_sec
    
    # use previous time step for more stable calculations
    temp = diagn.grid.temp_grid_prev
    humid = diagn.grid.humid_grid_prev
    geopotential = diagn.grid.geopotential
    temp_tend = diagn.tendencies.temp_tend_grid
    humid_tend = diagn.tendencies.humid_tend_grid
    pₛ = diagn.grid.pres_grid_prev[ij]                  # surface pressure [Pa]
    NF = eltype(temp)

    # thermodynamics
    ρ = atmosphere.water_density                        # density of water [kg/m³]
    g = planet.gravity                                  # gravity [m/s²]
    Lᵥ = clausius_clapeyron.latent_heat_condensation    # latent heat of vaporization
    cₚ = clausius_clapeyron.heat_capacity               # heat capacity

    # use scratch arrays for temp_ref_profile, humid_ref_profile
    temp_ref_profile =  diagn.dynamics.a_grid           # temperature [K] reference profile to adjust to
    humid_ref_profile = diagn.dynamics.b_grid           # specific humidity [kg/kg] profile to adjust to
    
    # CONVECTIVE CRITERIA AND FIRST GUESS RELAXATION
    level_zero_buoyancy = pseudo_adiabat!(ij, temp_ref_profile, temp, humid, geopotential, pₛ, σ, clausius_clapeyron)
            
    for k in level_zero_buoyancy:nlayers
        qsat = saturation_humidity(temp_ref_profile[ij, k], pₛ*σ[k], clausius_clapeyron)
        humid_ref_profile[ij, k] = qsat * convection.relative_humidity
    end

    Pq::NF = 0        # precipitation due to drying
    PT::NF = 0        # precipitation due to cooling
    ΔT::NF = 0        # vertically uniform temperature profile adjustment
    Qref::NF = 0      # = ∫_pₛ^p_LZB -humid_ref_profile dp

    # skip constants compared to Frierson 2007, i.e. no /τ, /gravity, *cₚ/Lᵥ
    for k in level_zero_buoyancy:nlayers
        # Frierson's equation (1)
        # δq = -(humid[ij, k] - humid_ref_profile[ij, k])/SBM.time_scale.value
        # Pq -= δq*Δσ[k]/gravity
        #
        # δT = -(temp[ij, k] - temp_ref_profile[ij, k])/SBM.time_scale.value
        # PT += δT*Δσ[k]/gravity*cₚ/Lᵥ

        # shorter form with same sign (τ, gravity, cₚ, Lᵥ all positive) to be reused
        Pq += (humid[ij, k] - humid_ref_profile[ij, k]) * Δσ[k]
        PT -= (temp[ij, k] - temp_ref_profile[ij, k]) * Δσ[k]
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
        # changing the reference profiles for both temperature and humidity so the
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

    # Initialize rain accumulation for this grid point
    rain_convection::NF = 0

    # GET TENDENCIES FROM ADJUSTED PROFILES
    τ⁻¹ = inv(convert(NF, Second(convection.time_scale).value))
    for k in level_zero_buoyancy:nlayers
        temp_tend[ij, k] -= (temp[ij, k] - temp_ref_profile[ij, k]) * τ⁻¹
        δq = (humid[ij, k] - humid_ref_profile[ij, k]) * τ⁻¹
        humid_tend[ij, k] -= δq

        # convective precipitation (rain), integrate dq\dt [(kg/kg)/s] vertically
        rain = max(δq * Δσ[k], 0)       # only integrate excess humidity for precip (no reevaporation)
        rain_convection += rain         # integrate vertically, Formula 25, unit [m]
    end
 
    # CONVECTIVE PRECIPITATION
    # enforce no precip for shallow conv via * deep_convection flag
    pₛΔt_gρ = (pₛ * Δt / (g * ρ)) * deep_convection 
	rain_convection *= pₛΔt_gρ                  # convert to [m] of rain during Δt
    rain_convection = max(rain_convection, 0)   # ensure non-negative precipitation, rounding errors
    
    # Store precipitation in diagnostic arrays
    diagn.physics.rain_convection[ij] += rain_convection            # accumulated rain [m] for output
    rain_rate_convection = rain_convection / Δt                     # instantaneous rate [m/s] for coupling
    diagn.physics.rain_rate_convection[ij] = rain_rate_convection   # instantaneous rate [m/s] for coupling
    
    # accumulate into total rain rate including large-scale condensation [m/s]
    diagn.physics.rain_rate[ij] += rain_rate_convection             # instantaneous rate [m/s] for coupling

    # clouds reach to top of convection
    diagn.physics.cloud_top[ij] = min(diagn.physics.cloud_top[ij], level_zero_buoyancy)
    return nothing
end

function variables(::BettsMillerConvection)
    return (
        DiagnosticVariable(name=:rain_convection, dims=Grid2D(), desc="Convective precipitation (accumulated)", units="m"),
        DiagnosticVariable(name=:rain_rate_convection, dims=Grid2D(), desc="Convective precipitation rate", units="m/s"),
        DiagnosticVariable(name=:rain_rate, dims=Grid2D(), desc="Rain rate (large-scale + convection)", units="m/s"),
        DiagnosticVariable(name=:cloud_top, dims=Grid2D(), desc="Cloud top layer index", units="1"),
    )
end

"""
$(TYPEDSIGNATURES)
Calculates the moist pseudo adiabat given temperature and humidity of surface parcel.
Follows the dry adiabat till condensation and then continues on the pseudo moist-adiabat
with immediate condensation to the level of zero buoyancy. Levels above are skipped,
set to NaN instead and should be skipped in the relaxation."""
@propagate_inbounds function pseudo_adiabat!(
    ij,
    temp_ref_profile,
    temp_environment,
    humid_environment,
    geopotential,
    pres,
    σ,
    clausius_clapeyron,
)
    NF = eltype(temp_ref_profile)                       # number format
    nlayers = length(σ)                                 # number of vertical layers

    # thermodynamics
    (; R_dry, R_vapour) = clausius_clapeyron
    Lᵥ = clausius_clapeyron.latent_heat_condensation    # latent heat of vaporization
    cₚ = clausius_clapeyron.heat_capacity               # heat capacity

    κ = R_dry/cₚ                         
    ε = clausius_clapeyron.mol_ratio
    μ = (1-ε)/ε                                         # for virtual temperature

    temp_parcel::NF = temp_environment[ij, nlayers]     # start at surface with environment temperature [K]
    humid_parcel::NF = humid_environment[ij, nlayers]   # and humidity [kg/kg]

    for k in 1:nlayers
        temp_ref_profile[ij, k] = NaN           # reset profile from any previous calculation, TODO necessary?
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
                
            ΔΦ = geopotential[ij, k] - geopotential[ij, k+1]        # vertical gradient in geopotential
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
        temp_virt_parcel = virtual_temperature(temp_parcel, humid_parcel, μ)         # virtual temperature of parcel
        buoyant = temp_virt_parcel > virtual_temperature(temp_environment[ij, k], humid_environment[ij, k], μ)     
    end
    
    # if parcel isn't buoyant anymore set last temperature (with negative buoyancy) back to NaN
    temp_ref_profile[ij, k] = !buoyant ? NaN32 : temp_ref_profile[ij, k]
    
    # level of zero buoyancy is reached when the loop stops, but in case it's at the top it's still buoyant
    level_zero_buoyancy = k + (1-buoyant)
    return level_zero_buoyancy
end

export BettsMillerDryConvection

"""
The simplified Betts-Miller convection scheme from Frierson, 2007,
https://doi.org/10.1175/JAS3935.1 but with humidity set to zero.
Fields and options are
$(TYPEDFIELDS)"""
@kwdef struct BettsMillerDryConvection{NF} <: AbstractConvection
    "[OPTION] Relaxation time for profile adjustment"
    time_scale::Second = Hour(4)
end

Adapt.@adapt_structure BettsMillerDryConvection

# generator function
BettsMillerDryConvection(SG::SpectralGrid; kwargs...) = BettsMillerDryConvection{SG.NF}(; kwargs...)
initialize!(::BettsMillerDryConvection, ::PrimitiveEquation) = nothing

# function barrier
@propagate_inbounds parameterization!(ij, diagn, progn, convection_scheme::BettsMillerDryConvection, model) =
    convection!(ij, diagn, convection_scheme, model)

"""
$(TYPEDSIGNATURES)
calculates temperature tendency for the dry convection scheme following the
simplified Betts-Miller convection from Frierson 2007 but with zero humidity.
Starts with a first-guess relaxation to determine the convective criterion,
then adjusts the reference profiles
for thermodynamic consistency (e.g. in dry convection the humidity profile is non-precipitating),
and relaxes current vertical profiles to the adjusted references."""
@propagate_inbounds function convection!(ij, diagn, DBM::BettsMillerDryConvection, model)

    (; geometry, atmosphere) = model
    NF = eltype(diagn.grid.temp_grid_prev)
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
@propagate_inbounds function dry_adiabat!(
    ij,
    temp_ref_profile,
    temp_environment,
    temp_parcel,
    σ,
    atmosphere,
)
    NF = eltype(temp_ref_profile)
    κ = atmosphere.heat_capacity

    nlayers = length(σ)                     # number of vertical levels

    for k in 1:nlayers
        temp_ref_profile[ij, k] = NaN       # reset profile from any previous calculation
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
@propagate_inbounds parameterization!(ij, diagn, progn, convection_scheme::ConvectiveHeating, model) =
    convection!(ij, diagn, convection_scheme, model)

@propagate_inbounds function convection!(
    ij,
    diagn::DiagnosticVariables,
    scheme::ConvectiveHeating,
    model,
)
    pₛ = diagn.grid.pres_grid_prev
    temp_tend = diagn.tendencies.temp_tend_grid
    nlayers = size(temp_tend, 2)
    NF = eltype(temp_tend)

    # Get latitude ring index and latitude
    j = model.geometry.whichring[ij]
    latd = model.geometry.latd[j]
    σ = model.geometry.σ_levels_full
    
    # escape immediately if not in the tropics
    abs(latd) >= scheme.σθ && return nothing

    p₀ = scheme.p₀ * 100     # hPa -> Pa
    σₚ = scheme.σₚ * 100     # hPa -> Pa
    cos²θ_term = scheme.lat_mask[j]
    Qmax = inv(convert(NF, Second(scheme.time_scale).value))

    for k in 1:nlayers
        p = pₛ[ij] * σ[k]      # Pressure in Pa on layer k

        # Lee and Kim, 2003, eq. 2
        temp_tend[ij, k] += Qmax*exp(-((p-p₀)/σₚ)^2 / 2)*cos²θ_term
    end
    return nothing
end