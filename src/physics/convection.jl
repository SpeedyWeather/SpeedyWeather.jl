abstract type AbstractConvection <: AbstractParameterization end

export NoConvection
struct NoConvection <: AbstractConvection end
NoConvection(::SpectralGrid) = NoConvection()
initialize!(::NoConvection, ::PrimitiveEquation) = nothing
convection!(::ColumnVariables, ::NoConvection, ::PrimitiveEquation) = nothing

export SimplifiedBettsMiller

"""
The simplified Betts-Miller convection scheme from Frierson, 2007,
https://doi.org/10.1175/JAS3935.1. This implements the qref-formulation
in their paper. Fields and options are
$(TYPEDFIELDS)"""
@kwdef struct SimplifiedBettsMiller{NF} <: AbstractConvection
    "number of vertical layers"
    nlayers::Int

    "[OPTION] Relaxation time for profile adjustment"
    time_scale::Second = Hour(4)

    "[OPTION] Relative humidity for reference profile"
    relative_humidity::NF = 0.7
end

SimplifiedBettsMiller(SG::SpectralGrid; kwargs...) = SimplifiedBettsMiller{SG.NF}(nlayers=SG.nlayers; kwargs...)
initialize!(::SimplifiedBettsMiller, ::PrimitiveEquation) = nothing

# function barrier for all AbstractConvection
function convection!(
    column::ColumnVariables,
    model::PrimitiveEquation,
)
    convection!(column, model.convection, model)
end

# function barrier to unpack model
function convection!(
    column::ColumnVariables,
    scheme::SimplifiedBettsMiller,
    model::PrimitiveWet,
)
    convection!(column, scheme, model.clausius_clapeyron,
                    model.geometry, model.planet, model.atmosphere, model.time_stepping)
end

"""
$(TYPEDSIGNATURES)
calculates temperature and humidity tendencies for the convection scheme following the
simplified Betts-Miller convection. Starts with a first-guess relaxation to determine
the convective criteria (none, dry/shallow or deep), then adjusts reference profiles
for thermodynamic consistency (e.g. in dry convection the humidity profile is non-precipitating),
and relaxes current vertical profiles to the adjusted references."""
function convection!(
    column::ColumnVariables{NF},
    SBM::SimplifiedBettsMiller,
    clausius_clapeyron::AbstractClausiusClapeyron,
    geometry::Geometry,
    planet::AbstractPlanet,
    atmosphere::AbstractAtmosphere,
    time_stepping::AbstractTimeStepper,
) where NF

    σ = geometry.σ_levels_full
    σ_half = geometry.σ_levels_half
    Δσ = geometry.σ_levels_thick
    (; geopot, nlayers, temp, temp_virt, humid, temp_tend, humid_tend) = column
    pₛ = column.pres[end]
    (; Lᵥ, cₚ) = clausius_clapeyron
    
    # use work arrays for temp_ref_profile, humid_ref_profile
    temp_ref_profile = column.a     # temperature [K] reference profile to adjust to
    humid_ref_profile = column.b    # specific humidity [kg/kg] profile to adjust to
    
    # CONVECTIVE CRITERIA AND FIRST GUESS RELAXATION
    temp_parcel = temp[nlayers]
    humid_parcel = humid[nlayers]
    level_zero_buoyancy = pseudo_adiabat!(temp_ref_profile,
                                            temp_parcel, humid_parcel,
                                            temp_virt, geopot, pₛ, σ,
                                            clausius_clapeyron)
            
    for k in level_zero_buoyancy:nlayers
        qsat = saturation_humidity(temp_ref_profile[k], pₛ*σ[k], clausius_clapeyron)
        humid_ref_profile[k] = qsat*SBM.relative_humidity
    end

    local Pq::NF = 0        # precipitation due to drying
    local PT::NF = 0        # precipitation due to cooling
    local ΔT::NF = 0        # vertically uniform temperature profile adjustment
    local Qref::NF = 0      # = ∫_pₛ^p_LZB -humid_ref_profile dp

    # skip constants compared to Frierson 2007, i.e. no /τ, /gravity, *cₚ/Lᵥ
    for k in level_zero_buoyancy:nlayers
        # Frierson's equation (1)
        # δq = -(humid[k] - humid_ref_profile[k])/SBM.time_scale.value
        # Pq -= δq*Δσ[k]/gravity
        #
        # δT = -(temp[k] - temp_ref_profile[k])/SBM.time_scale.value
        # PT += δT*Δσ[k]/gravity*cₚ/Lᵥ

        # shorter form with same sign (τ, gravity, cₚ, Lᵥ all positive) to be reused
        Pq += (humid[k] - humid_ref_profile[k])*Δσ[k]
        PT -= (temp[k] - temp_ref_profile[k])*Δσ[k]
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
            temp_ref_profile[k] -= ΔT       # equation (6)
        end
    
    elseif shallow_convection
        
        # FRIERSON'S QREF SCHEME
        # "changing the reference profiles for both temperature and humidity so the
        # precipitation is zero.

        for k in level_zero_buoyancy:nlayers
            Qref -= humid_ref_profile[k]*Δσ[k]  # eq (11) but in σ coordinates
        end
        fq = 1 - Pq/Qref                    # = 1 - Δq/Qref in eq (12) but we reuse Pq

        ΔT = PT/Δσ_lzb                      # equation (14), reuse PT and in σ coordinates
        for k in level_zero_buoyancy:nlayers
            humid_ref_profile[k] *= fq      # update humidity profile, eq (13)
            temp_ref_profile[k] -= ΔT       # update temperature profile, eq (15)
        end
    end

    # GET TENDENCIES FROM ADJUSTED PROFILES
    for k in level_zero_buoyancy:nlayers
        temp_tend[k] -= (temp[k] - temp_ref_profile[k]) / SBM.time_scale.value
        δq = (humid[k] - humid_ref_profile[k]) / SBM.time_scale.value
        humid_tend[k] -= δq

        # convective precipiation, integrate dq\dt [(kg/kg)/s] vertically
        column.precip_convection += δq * Δσ[k]
    end

    (; gravity) = planet
    (; water_density) = atmosphere
    (; Δt_sec) = time_stepping
    pₛΔt_gρ = (pₛ * Δt_sec / gravity / water_density) * deep_convection # enfore no precip for shallow conv 
    column.precip_convection *= pₛΔt_gρ                                 # convert to [m] of rain during Δt
    column.cloud_top = min(column.cloud_top, level_zero_buoyancy)       # clouds reach to top of convection
end

"""
$(TYPEDSIGNATURES)
Calculates the moist pseudo adiabat given temperature and humidity of surface parcel.
Follows the dry adiabat till condensation and then continues on the pseudo moist-adiabat
with immediate condensation to the level of zero buoyancy. Levels above are skipped,
set to NaN instead and should be skipped in the relaxation."""
function pseudo_adiabat!(
    temp_ref_profile::AbstractVector,
    temp_parcel::NF,
    humid_parcel::Real,
    temp_virt_environment::AbstractVector,
    geopot::AbstractVector,
    pres::Real,
    σ::AbstractVector,
    clausius_clapeyron::AbstractClausiusClapeyron,
) where NF

    (; Lᵥ, R_dry, R_vapour, cₚ) = clausius_clapeyron
    R_cₚ = R_dry/cₚ
    ε = clausius_clapeyron.mol_ratio
    μ = (1-ε)/ε                             # for virtual temperature

    @boundscheck length(temp_ref_profile) == length(geopot) ==
        length(σ) == length(temp_virt_environment) || throw(BoundsError)

    nlayers = length(temp_ref_profile)      # number of vertical levels
    temp_ref_profile .= NaN                 # reset profile from any previous calculation
    temp_ref_profile[nlayers] = temp_parcel # start profile with parcel temperature

    local saturated::Bool = false           # did the parcel reach saturation yet?
    local buoyant::Bool = true              # is the parcel still buoyant?
    local k::Int = nlayers                  # layer index top to surface
    local temp_virt_parcel::NF = temp_parcel * (1 + μ*humid_parcel)

    while buoyant && k > 1                  # calculate moist adiabat while buoyant till top
        k -= 1                              # one level up
            
        if !saturated                       # if not saturated yet follow dry adiabat
            # dry adiabatic ascent and saturation humidity of that temperature 
            temp_parcel_dry = temp_parcel*(σ[k]/σ[k+1])^R_cₚ
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
                
            ΔΦ = geopot[k] - geopot[k+1]                            # vertical gradient in geopotential
            temp_parcel = temp_parcel - ΔΦ/cₚ*Γ                     # new temperature of parcel at k
                
            # at new (lower) temperature condensation occurs immediately
            # new humidity equals to that saturation humidity
            humid_parcel = saturation_humidity(temp_parcel, σ[k]*pres, clausius_clapeyron)
        else
            temp_parcel = temp_parcel_dry       # else parcel temperature following dry adiabat
        end
    
        # use dry/moist adiabatic ascent for reference profile
        temp_ref_profile[k] = temp_parcel

        # check whether parcel is still buoyant wrt to environment
        # use virtual temperature as it's equivalent to density
        temp_virt_parcel = temp_parcel*(1 + μ*humid_parcel)         # virtual temperature of parcel
        buoyant = temp_virt_parcel > temp_virt_environment[k] ? true : false      
    end
    
    # if parcel isn't buoyant anymore set last temperature (with negative buoyancy) back to NaN
    temp_ref_profile[k] = !buoyant ? NaN : temp_ref_profile[k]    
    
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
    "number of vertical layers/levels"
    nlayers::Int

    "[OPTION] Relaxation time for profile adjustment"
    time_scale::Second = Hour(4)
end

DryBettsMiller(SG::SpectralGrid; kwargs...) = DryBettsMiller{SG.NF}(nlayers=SG.nlayers; kwargs...)
initialize!(::DryBettsMiller, ::PrimitiveEquation) = nothing

# function barrier to unpack model
function convection!(
    column::ColumnVariables,
    scheme::DryBettsMiller,
    model::PrimitiveEquation,
)
    convection!(column, scheme, model.geometry, model.atmosphere)
end

"""
$(TYPEDSIGNATURES)
calculates temperature tendency for the dry convection scheme following the
simplified Betts-Miller convection from Frierson 2007 but with zero humidity.
Starts with a first-guess relaxation to determine the convective criterion,
then adjusts the reference profiles
for thermodynamic consistency (e.g. in dry convection the humidity profile is non-precipitating),
and relaxes current vertical profiles to the adjusted references."""
function convection!(
    column::ColumnVariables{NF},
    DBM::DryBettsMiller,
    geometry::Geometry,
    atmosphere::AbstractAtmosphere,
) where NF

    σ = geometry.σ_levels_full
    σ_half = geometry.σ_levels_half
    Δσ = geometry.σ_levels_thick
    (; nlayers, temp, temp_tend) = column

    # use work arrays for temp_ref_profile
    temp_ref_profile = column.a     # temperature [K] reference profile to adjust to

    # CONVECTIVE CRITERIA AND FIRST GUESS RELAXATION
    level_zero_buoyancy = dry_adiabat!(temp_ref_profile,
                                            temp, σ,
                                            atmosphere)

    local PT::NF = 0        # precipitation due to coolinga
    local ΔT::NF = 0        # vertically uniform temperature profile adjustment

    # skip constants compared to Frierson 2007, i.e. no /τ, /gravity, *cₚ/Lᵥ
    for k in level_zero_buoyancy:nlayers
        # Frierson's equation (1)
        # δT = -(temp[k] - temp_ref_profile[k])/SBM.time_scale.value
        # PT += δT*Δσ[k]/gravity*cₚ/Lᵥ

        # shorter form with same sign (τ, gravity, cₚ, Lᵥ all positive) to be reused
        PT -= (temp[k] - temp_ref_profile[k])*Δσ[k]
    end

    # ADJUST PROFILES FOLLOWING FRIERSON 2007
    convection = PT > 0
    convection || return nothing            # escape immediately for no convection

    # height of zero buoyancy level in σ coordinates
    Δσ_lzb = σ_half[nlayers+1] - σ_half[level_zero_buoyancy]   
    ΔT = PT/Δσ_lzb                          # eq (5) or (14) but reusing PT
    for k in level_zero_buoyancy:nlayers
        temp_ref_profile[k] -= ΔT           # equation (6) or equation (15)
        temp_tend[k] -= (temp[k] - temp_ref_profile[k]) / DBM.time_scale.value
    end
end

"""
$(TYPEDSIGNATURES)
Calculates the moist pseudo adiabat given temperature and humidity of surface parcel.
Follows the dry adiabat till condensation and then continues on the pseudo moist-adiabat
with immediate condensation to the level of zero buoyancy. Levels above are skipped,
set to NaN instead and should be skipped in the relaxation."""
function dry_adiabat!(
    temp_ref_profile::AbstractVector,
    temp_environment::AbstractVector,
    σ::AbstractVector,
    atmosphere::AbstractAtmosphere,
)

    cₚ = atmosphere.heat_capacity
    R_cₚ = atmosphere.R_dry/cₚ

    @boundscheck length(temp_ref_profile) ==
        length(σ) == length(temp_environment) || throw(BoundsError)

    nlayers = length(temp_ref_profile)      # number of vertical levels
    temp_parcel = temp_environment[nlayers] # parcel is at lowermost layer temperature
    temp_ref_profile .= NaN                 # reset profile from any previous calculation
    temp_ref_profile[nlayers] = temp_parcel    # start profile with parcel temperature

    local buoyant::Bool = true              # is the parcel still buoyant?
    local k::Int = nlayers                  # layer index top to surface

    while buoyant && k > 1                  # calculate moist adiabat while buoyant till top
        k -= 1                              # one level up
            
        # dry adiabatic ascent
        temp_parcel = temp_parcel*(σ[k]/σ[k+1])^R_cₚ
        temp_ref_profile[k] = temp_parcel

        # check whether parcel is still buoyant wrt to environment
        buoyant = temp_parcel > temp_environment[k] ? true : false      
    end
    
    # if parcel isn't buoyant anymore set last temperature (with negative buoyancy) back to NaN
    temp_ref_profile[k] = !buoyant ? NaN : temp_ref_profile[k]    
    
    # level of zero buoyancy is reached when the loop stops, but in case it's at the top it's still buoyant
    level_zero_buoyancy = k + (1-buoyant)
    return level_zero_buoyancy
end