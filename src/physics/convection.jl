abstract type AbstractConvection{NF} <: AbstractParameterization{NF} end

Base.@kwdef struct SpeedyConvection{NF} <: AbstractConvection{NF}

    "Number of vertical levels"
    nlev::Int

    "Minimum (normalised) surface pressure for the occurrence of convection [1]"
    pres_threshold::NF = 0.8

    "Relative humidity threshold for convection in PBL [1]"
    humid_threshold_boundary::NF = 0.9

    "Relative humidity threshold for convection in the troposphere [1]"
    humid_threshold_troposphere::NF = 0.7

    "Relaxation time for PBL humidity"
    time_scale::Second = Hour(6)

    "Maximum entrainment as a fraction of cloud-base mass flux"
    max_entrainment::NF = 0.5

    "Ratio between secondary and primary mass flux at cloud-base"
    ratio_secondary_mass_flux::NF = 0.8

    "Super saturation of relative humidity in the planteray boundar layer [1]"
    super_saturation::NF = 1.01

    "Mass flux limiter for humidity [kg/kg]"
    mass_flux_limiter_humid::NF = 5e-3

    # precomputed in initialize!
    "Reference surface pressure [Pa]"
    pres_ref::Base.RefValue{NF} = Base.Ref(zero(NF))

    "latent heat of condensation [J/kg] for consistency with specific humidity [kg/kg], also called alhc"
    latent_heat_condensation::Base.RefValue{NF} = Base.Ref(zero(NF))

    "specific heat [J/kg/K] of air"
    specific_heat::Base.RefValue{NF} = Base.Ref(zero(NF))

    "Number of vertical levels for stratosphere"
    n_stratosphere_levels::Base.RefValue{Int} = Base.Ref(0)

    "Number of vertical levels for planetary boundary layer"
    n_boundary_levels::Base.RefValue{Int} = Base.Ref(0)

    "Mass flux entrainment profile in the vertical [1]"
    entrainment_profile::Vector{NF} = zeros(NF,nlev)
end

SpeedyConvection(SG::SpectralGrid;kwargs...) = SpeedyConvection{SG.NF}(nlev=SG.nlev;kwargs...)

function initialize!(convection::SpeedyConvection,model::PrimitiveWet)

    (;σ_levels_full) = model.geometry
    (;σ_tropopause, σ_boundary_layer) = model.atmosphere
    (;entrainment_profile, nlev) = convection

    # number of stratospheric levels, for nlev very small n can be nothing, then use 0
    n = findlast(σ->σ<=σ_tropopause,σ_levels_full)
    convection.n_stratosphere_levels[] = isnothing(n) ? 0 : n

    # number of levels for the planetary boundary layer, same as above
    n = findlast(σ->σ<=σ_boundary_layer,σ_levels_full)
    convection.n_boundary_levels[] = isnothing(n) ? 0 : nlev - n
    
    # reference pressure
    convection.pres_ref[] = model.atmosphere.pres_ref*100     # [hPa] -> [Pa]
    convection.latent_heat_condensation[] = model.atmosphere.latent_heat_condensation
    convection.specific_heat[] = model.atmosphere.cₚ

    # Mass entrainment profile
    entrainment_profile[1] = 0      # no entrainment in top layer
    entrainment_profile[nlev] = 0   # no entrainment in bottom layer
    for k = 2:nlev-1                # intermediate layers with minimum at σ=0
        entrainment_profile[k] = max(0, (σ_levels_full[k] - 0.5))
    end

    # profile as fraction of cloud-base mass flux, normalise to max entrainment at nlev-1
    entrainment_profile .*= convection.max_entrainment/entrainment_profile[nlev-1]
end

# dry model doesn't have convection TODO maybe it should?
convection!(column::ColumnVariables,model::PrimitiveDry) = nothing

# function barrier for all AbstractConvection
function convection!(
    column::ColumnVariables,
    model::PrimitiveEquation,
)
    convection!(column,model.convection,model)
end

function convection!(
    column::ColumnVariables,
    scheme::SpeedyConvection,
    model::PrimitiveEquation,
)

    moist_static_energy!(column, model.clausis_clapeyron)
    vertical_interpolate!(column, model)                # Interpolate certain variables to half-levels

    # always diagnose convection
    diagnose_convection!(column, model.convection)

    # but only execute if conditions are met
    if column.conditional_instability && column.activate_convection
        convection!(column,model.convection,model.constants,model.geometry,model.time_stepping)
    end
end

"""
$(TYPEDSIGNATURES)
Check whether the convection scheme should be activated in the given atmospheric column.

1. A conditional instability exists when the saturation moist energy (MSS) decreases with
height, that is, there exists an atmospheric level k such that,

    MSS(N) > MSS(k+h)

where N is the planetary boundary layer (PBL) and k+h is the half-level at the lower
boundary of the full level k.

2. When a conditional instability exists, the convection scheme is activated when, either,

    a. the actual moist static energy (MSE) at level N-1 (directly above the PBL) is greater
       than the saturation moist static energy at some half-level k+h,

            MSE(N-1) > MSS(k+h)

    b. the humidity in both the PBL and one layer above exceeds a prescribed threshold,

            Q(N)   > RH_cnv * Qˢᵃᵗ(N)
            Q(N-1) > RH_cnv * Qˢᵃᵗ(N-1)

The top-of-convection (TCN) layer, or cloud-top, is the largest value of k for which
condition 1 is satisfied. The cloud-top layer may be subsequently adjusted upwards by the
large-scale condensation parameterization, which is executed after this one."""
function diagnose_convection!(column::ColumnVariables,convection::SpeedyConvection)

    (; pres_ref, pres_threshold, humid_threshold_boundary) = convection
    n_stratosphere_levels = convection.n_stratosphere_levels[]
    n_boundary_levels = convection.n_boundary_levels[]
    latent_heat = convection.latent_heat_condensation[]

    (; nlev ) = column
    (; humid, pres, sat_humid, moist_static_energy,
    sat_moist_static_energy, sat_moist_static_energy_half) = column

    # effectively disables convection over Himalaya/Tibet, Greenland and Antarctica
    if pres[end] > pres_threshold*pres_ref[]

        # First we pre-compute some values which we will need inside the loop
        # 1. Saturation (or super-saturated) moist static energy in the PBL
        sat_moist_static_energy_pbl =
            max(moist_static_energy[nlev], sat_moist_static_energy[nlev])

        # 2. Minimum of moist static energy in the lowest two levels
        moist_static_energy_lower_trop =
            min(moist_static_energy[nlev], moist_static_energy[nlev-1])

        # 3. Humidity threshold for convection, defined in the PBL and one level above
        humid_threshold_pbl = humid_threshold_boundary * sat_humid[nlev]
        humid_threshold_above_pbl = humid_threshold_boundary * sat_humid[nlev-1]

        # The range of this loop requires clarification, but in its current form it means
        # that the top-of-convection level may be any tropospheric level, excluding the two
        # layers directly above the PBL.
        for k = (nlev-(n_boundary_levels+1)):-1:(n_stratosphere_levels+1)
            # Condition 1: Conditional instability (MSS in PBL < MSS at this half-level)
            if sat_moist_static_energy_pbl > sat_moist_static_energy_half[k]
                column.conditional_instability = true
                column.cloud_top = k
            end

            # Condition 2a: Gradient of actual moist static energy between lower and upper troposphere
            if moist_static_energy_lower_trop > sat_moist_static_energy_half[k]
                column.activate_convection = true
                column.excess_humid = max(
                    humid[nlev] - humid_threshold_pbl,
                    (moist_static_energy[nlev] - sat_moist_static_energy_half[k]) / latent_heat,
                )
            end
        end

        if column.conditional_instability && column.activate_convection
            return nothing  # Condition for convection already satisfied
        end

        # Condition 2b: Humidity exceeds threshold in both PBL and one layer above
        if column.conditional_instability &&
           (humid[nlev] > humid_threshold_pbl) &&
           (humid[nlev-1] > humid_threshold_above_pbl)
            column.activate_convection = true
            column.excess_humid = humid[nlev] - humid_threshold_pbl
        end
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
Compute fluxes and precipitation due to convection in the given atmospheric column.

The scheme computes fluxes of mass, humidity and dry static energy. A part of the upward
moisture flux at the lower boundary of the cloud-top (TCN) layer is converted into
convective precipitation.

For full details of the scheme see: http://users.ictp.it/~kucharsk/speedy_description/km_ver41_appendixA.pdf
"""
function convection!(
    column::ColumnVariables{NF},
    convection::SpeedyConvection,
    C::DynamicsConstants,
    G::Geometry,
    T::TimeStepper,
) where NF

    (; gravity, pres_ref ) = C
    (; nlev, pres, humid, humid_half, sat_humid, sat_humid_half, cloud_top) = column
    (; dry_static_energy, dry_static_energy_half, cloud_top) = column
    (; flux_temp_downward, flux_temp_upward, flux_humid_downward, flux_humid_upward) = column
    (; temp_tend, humid_tend) = column
    (; σ_levels_thick) = G

    n_boundary_levels = convection.n_boundary_levels[]
    time_scale = convection.time_scale.value
    latent_heat_cₚ = convection.latent_heat_condensation[]/convection.specific_heat[]
    (;pres_threshold, humid_threshold_troposphere, ratio_secondary_mass_flux) = convection
    (;mass_flux_limiter_humid) = convection

    # 1. Fluxes in the planetary boundary layer (PBL)
    # Humidity at the upper boundary of the PBL
    humid_top_of_pbl = min(humid_half[nlev-n_boundary_levels], humid[nlev-n_boundary_levels+1])
    
    # Maximum specific humidity in the PBL
    max_humid_pbl = zero(NF)
    for k in nlev-n_boundary_levels+1:nlev
        max_humid_pbl = max(convection.super_saturation * humid[k], sat_humid[k])
    end

    # Cloud-base mass flux
    pₛ = pres[end]                          # surface pressure
    p_norm = pₛ/pres_ref                    # normalised surface pressure
    Δp = pres_ref*σ_levels_thick[nlev]      # Pressure difference between bottom and top of PBL
    
    # excess humidity relative to humid difference across PBL
    excess_humid = column.excess_humid / (max_humid_pbl - humid_top_of_pbl)

    # Fortran SPEEDY documentation eq. (12) and (13)
    mass_flux =
        Δp / (gravity * time_scale) *
        # Fortran SPEEDY extends this with some flux limiters:
        # (original SPEEDY formulation doesn't make sense to me given it's supposed to be for
        # "Minimum (normalised) surface pressure for the occurrence of convection")
        # p_norm*min(1,  2 * (p_norm - pres_threshold) / (1 - pres_threshold)) *
        max(0,(p_norm - pres_threshold) / (1 - pres_threshold)) *       # new version
        min(mass_flux_limiter_humid, excess_humid)
    
    column.cloud_base_mass_flux = mass_flux

    # Upward fluxes at upper boundary of PBL
    flux_up_humid = mass_flux * max_humid_pbl                                       # eq. (10)
    flux_up_static_energy = mass_flux * dry_static_energy[nlev-n_boundary_levels+1] # eq. (11)

    # Downward fluxes at upper boundary of PBL
    flux_down_humid = mass_flux * humid_top_of_pbl                                          # eq. (10) 
    flux_down_static_energy = mass_flux * dry_static_energy_half[nlev-n_boundary_levels]    # eq. (11)

    # Accumulate in fluxes for all parameterizations
    flux_temp_upward[nlev-n_boundary_levels+1]   += flux_up_static_energy
    flux_temp_downward[nlev-n_boundary_levels+1] += flux_down_static_energy
    
    flux_humid_upward[nlev-n_boundary_levels+1]   += flux_up_humid
    flux_humid_downward[nlev-n_boundary_levels+1] += flux_down_humid

    # 2. Fluxes for intermediate layers
    for k = (nlev-n_boundary_levels):-1:(cloud_top+1)

        # Mass entrainment
        mass_entrainment = p_norm * convection.entrainment_profile[k] * column.cloud_base_mass_flux
        mass_flux += mass_entrainment

        # Upward fluxes at upper boundary
        flux_up_static_energy += mass_entrainment * dry_static_energy[k]
        flux_up_humid         += mass_entrainment * humid[k]

        # Downward fluxes at upper boundary, _half vectors skip top (k=1/2)
        flux_down_static_energy = mass_flux * dry_static_energy_half[k-1]
        flux_down_humid         = mass_flux * humid_half[k-1]

        # accumulate in fluxes that are translated to tendencies in tendencies.jl
        flux_temp_upward[k]   += flux_up_static_energy
        flux_temp_downward[k] += flux_down_static_energy
        
        flux_humid_upward[k]   += flux_up_humid
        flux_humid_downward[k] += flux_down_humid

        # Secondary moisture flux representing shallower, non-precipitating convective systems
        # Occurs when RH in an intermediate layer falls below a threshold
        Δhumid = humid_threshold_troposphere * sat_humid[k] - humid[k]
        if Δhumid > 0
            Δflux_humid = ratio_secondary_mass_flux * mass_flux * Δhumid

            # a flux from bottom nlev to layer k, equivalent to net flux of Δflux_humid into k but out of nlev
            for kk in k+1:nlev+1
                flux_humid_upward[kk] += Δflux_humid
            end
        end
    end

    # 3. Convective precipitation in top-of-convection layer
    precip_convection = max(flux_up_humid - mass_flux * sat_humid_half[cloud_top], 0)   # in [kg/m²/s]
    column.precip_convection = precip_convection*T.Δt_sec/C.water_density               # convert to [m]

    # Condensation of convectice precipiation creates a humidity and temperature (heating)
    # tendency at the top of convection layer cloud_top
    humid_tend_k = -precip_convection*gravity/(pₛ*σ_levels_thick[cloud_top])
    humid_tend[cloud_top] += humid_tend_k
    temp_tend[cloud_top] -= humid_tend_k*latent_heat_cₚ

    return nothing
end

Base.@kwdef struct SimplifiedBettsMiller{NF} <: AbstractConvection{NF}
    "number of vertical layers/levels"
    nlev::Int

    "Relaxation time for profile adjustment"
    time_scale::Second = Hour(4)

    "Relative humidity for reference profile"
    relative_humidity::NF = 0.8

    "temperature [K] reference profile to adjust to"
    temp_ref_profile::Vector{NF} = zeros(NF,nlev)

    "specific humidity [kg/kg] profile to adjust to"
    humid_ref_profile::Vector{NF} = zeros(NF,nlev)
end

SimplifiedBettsMiller(SG::SpectralGrid;kwargs...) = SimplifiedBettsMiller{SG.NF}(nlev=SG.nlev;kwargs...)

initialize!(::SimplifiedBettsMiller,::PrimitiveWet) = nothing

# function barrier to unpack model
function convection!(
    column::ColumnVariables,
    scheme::SimplifiedBettsMiller,
    model::PrimitiveEquation,
)
    convection!(column, scheme, model.clausis_clapeyron, model.geometry, model.constants, model.time_stepping)
end

function convection!(
    column::ColumnVariables{NF},
    SBM::SimplifiedBettsMiller,
    clausius_clapeyron::AbstractClausiusClapeyron,
    geometry::Geometry,
    constants::DynamicsConstants,
    time_stepping::TimeStepper,
) where NF

    σ = geometry.σ_levels_full
    σ_half = geometry.σ_levels_half
    Δσ = geometry.σ_levels_thick
    (;temp_ref_profile, humid_ref_profile) = SBM
    (;geopot, nlev, temp, temp_virt, humid, temp_tend, humid_tend) = column
    pₛ = column.pres[end]
    (;Lᵥ, cₚ) = clausius_clapeyron

    # CONVECTIVE CRITERIA AND FIRST GUESS RELAXATION
    temp_parcel = temp[nlev]
    humid_parcel = humid[nlev]
    level_zero_buoyancy = pseudo_adiabat!(temp_ref_profile,
                                            temp_parcel, humid_parcel,
                                            temp_virt, geopot, pₛ, σ,
                                            clausius_clapeyron)
            
    for k in level_zero_buoyancy:nlev
        qsat = saturation_humidity(temp_ref_profile[k],pₛ*σ[k],clausius_clapeyron)
        humid_ref_profile[k] = qsat*SBM.relative_humidity
    end

    local Pq::NF = 0        # precipitation due to drying
    local PT::NF = 0        # precipitation due to cooling
    local ΔT::NF = 0        # vertically uniform temperature profile adjustment
    local Qref::NF = 0      # = ∫_pₛ^p_LZB -humid_ref_profile dp

    # skip constants compared to Frierson 2007, i.e. no /τ, /gravity, *cₚ/Lᵥ
    for k in level_zero_buoyancy:nlev
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

    # escape immediately for no convection
    no_convection = !(deep_convection || shallow_convection)
    no_convection && return nothing

    # height of zero buoyancy level in σ coordinates
    Δσ_lzb = σ_half[nlev+1] - σ_half[level_zero_buoyancy]   

    if deep_convection

        ΔT = (PT + Pq*Lᵥ/cₚ)/Δσ_lzb         # minus eq (5) but reusing PT, Pq, and /cₚ already included

        for k in level_zero_buoyancy:nlev
            temp_ref_profile[k] += ΔT       # equation (6)
        end
    
    elseif shallow_convection
        
        # FRIERSON'S QREF SCHEME
        # "changing the reference profiles for both temperature and humidity so the
        # precipitation is zero.

        for k in level_zero_buoyancy:nlev
            Qref -= humid_ref_profile[k]*Δσ[k]  # eq (11) but in σ coordinates
        end
        fq = 1 - Pq/Qref                    # = 1 - Δq/Qref in eq (12) but we reuse Pq

        ΔT = PT/Δσ_lzb                      # equation (14), reuse PT and in σ coordinates
        for k in level_zero_buoyancy:nlev
            humid_ref_profile[k] *= fq      # update humidity profile, eq (13)
            temp_ref_profile[k] -= ΔT       # update temperature profile, eq (15)
        end
    end

    # GET TENDENCIES FROM ADJUSTED PROFILES
    for k in level_zero_buoyancy:nlev
        temp_tend[k] -= (temp[k] - temp_ref_profile[k]) / SBM.time_scale.value
        δq = (humid[k] - humid_ref_profile[k]) / SBM.time_scale.value
        humid_tend[k] -= δq

        # convective precipiation, integrate dq\dt [(kg/kg)/s] vertically
        column.precip_convection += δq * Δσ[k]
    end

    (;gravity, water_density) = constants
    (;Δt_sec) = time_stepping
    pₛΔt_gρ = (pₛ * Δt_sec / gravity / water_density) * deep_convection # enfore no precip for shallow conv 
    column.precip_convection *= pₛΔt_gρ                                 # convert to [m] of rain during Δt
    
    # TODO at the moment the cloud top include shallow convection, should that be?
    column.cloud_top = min(column.cloud_top,level_zero_buoyancy)        # clouds reach to top of convection
end

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

    (;Lᵥ, R_dry, R_vapour, cₚ) = clausius_clapeyron
    R_cₚ = R_dry/cₚ
    ε = clausius_clapeyron.mol_ratio
    μ = (1-ε)/ε                             # for virtual temperature

    @boundscheck length(temp_ref_profile) == length(geopot) ==
        length(σ) == length(temp_virt_environment) || throw(BoundsError)

    nlev = length(temp_ref_profile)         # number of vertical levels
    temp_ref_profile .= NaN                 # reset profile from any previous calculation
    temp_ref_profile[nlev] = temp_parcel    # start profile with parcel temperature

    local saturated::Bool = false           # did the parcel reach saturation yet?
    local buoyant::Bool = true              # is the parcel still buoyant?
    local k::Int = nlev                     # layer index top to surface
    local temp_virt_parcel::NF = temp_parcel * (1 + μ*humid_parcel)

    while buoyant && k > 1                  # calculate moist adiabat while buoyant till top
        k -= 1                              # one level up
            
        if !saturated                       # if not saturated yet follow dry adiabat
            # dry adiabatic ascent and saturation humidity of that temperature 
            temp_parcel_dry = temp_parcel*(σ[k]/σ[k+1])^R_cₚ
            sat_humid = saturation_humidity(temp_parcel_dry,σ[k]*pres,clausius_clapeyron)
                    
            # set to saturated when the dry adiabatic ascent would reach saturation
            # then follow moist adiabat instead (see below)
            saturated = humid_parcel >= sat_humid
        end
    
        if saturated            
            # calculate moist/pseudo adiabatic lapse rate, dT/dΦ = -Γ/cp
            T, Tᵥ, q = temp_parcel, temp_virt_parcel, humid_parcel  # for brevity
            A = q*Lᵥ / ((1-q)^2 * R_dry)
            B = q*Lᵥ^2 / ((1-q)^2 * cₚ * R_vapour)
            Γ = (1 + A/Tᵥ) / (1 + B/T^2)
                
            ΔΦ = geopot[k] - geopot[k+1]                            # vertical gradient in geopotential
            temp_parcel = temp_parcel - ΔΦ/cₚ*Γ                     # new temperature of parcel at k
                
            # at new (lower) temperature condensation occurs immediately
            # new humidity equals to that saturation humidity
            humid_parcel = saturation_humidity(temp_parcel,σ[k]*pres,clausius_clapeyron)
        else
            temp_parcel = temp_parcel_dry       # else parcel temperature following dry adiabat
        end
    
        # use dry/moist adiabatic ascent for reference profile
        temp_ref_profile[k] = temp_parcel

        # check whether parcel is still buoyant wrt to environment
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