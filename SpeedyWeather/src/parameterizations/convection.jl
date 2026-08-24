abstract type AbstractConvection <: AbstractParameterization end

export BettsMillerConvection

"""The simplified Betts-Miller convection scheme from Frierson, 2007,
https://doi.org/10.1175/JAS3935.1. This implements the qref-formulation
in their paper. Fields and options are $(TYPEDFIELDS)"""
@parameterized @kwdef struct BettsMillerConvection{NF, Entrainment <: AbstractEntrainment} <: AbstractConvection
    "[OPTION] Relaxation time for profile adjustment"
    time_scale::Second = Hour(4)

    "[OPTION] Relative humidity for reference profile [1]"
    @param relative_humidity::NF = 0.7

    "[OPTION] Entrainment profile mixing environmental air into the rising parcel"
    @component entrainment::Entrainment = NoEntrainment()

    "[OPTION] Convert convective rain below freezing to snow?"
    snow::Bool = true

    "[OPTION] Temperature of the lowermost layer below which convective rain becomes snow [K]"
    @param freezing_threshold::NF = 273.15 (bounds = Positive,)
end

Adapt.@adapt_structure BettsMillerConvection

# generator function
BettsMillerConvection(SG::SpectralGrid; entrainment = NoEntrainment(), kwargs...) =
    BettsMillerConvection{SG.NF, typeof(entrainment)}(; entrainment, kwargs...)
initialize!(::BettsMillerConvection, ::PrimitiveEquation) = nothing

# function barrier
@propagate_inbounds parameterization!(ij, vars, convection_scheme::BettsMillerConvection, model) =
    convection!(ij, vars, convection_scheme, model)

"""
$(TYPEDSIGNATURES)
calculates temperature and humidity tendencies for the convection scheme following the
simplified Betts-Miller convection. Starts with a first-guess relaxation to determine
the convective criteria (none, dry/shallow or deep), then adjusts reference profiles
for thermodynamic consistency (e.g. in dry convection the humidity profile is non-precipitating),
and relaxes current vertical profiles to the adjusted references."""
@propagate_inbounds function convection!(ij, vars, convection::BettsMillerConvection, model)

    (; geometry, planet, atmosphere, time_stepping) = model
    # TODO: σ, σ_half, Δσ are used for buoyancy level detection and moisture flux
    # calculations in the Betts-Miller scheme. These are baked into sigma-coordinate
    # formulations and would need revisiting for hybrid coordinates.
    σ = geometry.σ_levels_full
    σ_half = geometry.σ_levels_half
    Δσ = geometry.σ_levels_thick
    nlayers = length(σ)
    (; Δt) = model.time_stepping                            # time step in [s]

    # use previous time step for more stable calculations
    temp = get_prognostic_step(vars.grid.temperature, time_stepping, convection)
    humid = get_prognostic_step(vars.grid.humidity, time_stepping, convection)
    geopotential = vars.dynamics.geopotential
    temp_tend = get_tendency_step(vars.tendencies.grid.temperature, time_stepping, convection)
    humid_tend = get_tendency_step(vars.tendencies.grid.humidity, time_stepping, convection)
    pₛ = vars.parameterizations.surface_pressure[ij]        # surface pressure [Pa]
    NF = eltype(temp)

    # thermodynamics
    ρ = atmosphere.water_density                # density of water [kg/m³]
    g = planet.gravity                          # gravity [m/s²]
    Lᵥ = atmosphere.latent_heat_condensation    # latent heat of vaporization
    cₚ = atmosphere.heat_capacity               # heat capacity

    # use scratch arrays for temp_ref_profile, humid_ref_profile
    temp_ref_profile = vars.scratch.grid.a      # temperature [K] reference profile to adjust to
    humid_ref_profile = vars.scratch.grid.b     # specific humidity [kg/kg] profile to adjust to

    # CONVECTIVE CRITERIA AND FIRST GUESS RELAXATION
    # pseudo_adiabat! also fills humid_ref_profile; above the level of zero buoyancy (LZB)
    # both reference profiles equal the environment, so they contribute exact zero below
    # without needing to mask the loop range to level_zero_buoyancy:nlayers
    level_zero_buoyancy = pseudo_adiabat!(
        ij, temp_ref_profile, humid_ref_profile, temp, humid,
        geopotential, pₛ, σ, atmosphere, convection.relative_humidity, convection.entrainment,
    )

    Pq::NF = 0        # precipitation due to drying
    PT::NF = 0        # precipitation due to cooling
    Qref::NF = 0      # = ∫_pₛ^p_LZB -humid_ref_profile dp

    # skip constants compared to Frierson 2007, i.e. no /τ, /gravity, *cₚ/Lᵥ
    # full-range, branchless: above LZB, temp_ref_profile == temp and humid_ref_profile == humid
    # (see pseudo_adiabat!) so Pq, PT accumulate exact zero there without a mask; Qref needs a
    # mask since humid_ref_profile above LZB is a real (environmental) humidity, not zero
    for k in 1:nlayers
        # Frierson's equation (1)
        # δq = -(humid[ij, k] - humid_ref_profile[ij, k])/SBM.time_scale.value
        # Pq -= δq*Δσ[k]/gravity
        #
        # δT = -(temp[ij, k] - temp_ref_profile[ij, k])/SBM.time_scale.value
        # PT += δT*Δσ[k]/gravity*cₚ/Lᵥ

        # shorter form with same sign (τ, gravity, cₚ, Lᵥ all positive) to be reused
        Pq += (humid[ij, k] - humid_ref_profile[ij, k]) * Δσ[k]
        PT -= (temp[ij, k] - temp_ref_profile[ij, k]) * Δσ[k]

        in_column = k >= level_zero_buoyancy
        Qref -= ifelse(in_column, humid_ref_profile[ij, k], zero(NF)) * Δσ[k]  # eq (11) but in σ coordinates
    end

    # ADJUST PROFILES FOLLOWING FRIERSON 2007
    deep_convection = Pq > 0 && PT > 0
    shallow_convection = Pq <= 0 && PT > 0

    # escape immediately for no convection
    no_convection = !(deep_convection || shallow_convection)
    no_convection && return nothing

    # height of zero buoyancy level in σ coordinates
    Δσ_lzb = σ_half[nlayers + 1] - σ_half[level_zero_buoyancy]

    # branchless selection between the deep (eq 5-6) and shallow "qref" (eq 11-15) corrections;
    # fq = 1 for deep convection, i.e. no humidity profile correction there
    ΔT = ifelse(deep_convection, (PT - Pq * Lᵥ / cₚ) / Δσ_lzb, PT / Δσ_lzb)
    fq = ifelse(deep_convection, one(NF), 1 - Pq / Qref)

    # Initialize rain accumulation for this grid point
    rain_convection::NF = 0

    # GET TENDENCIES FROM ADJUSTED PROFILES
    # the corrected reference profiles are computed inline (not written back to the scratch
    # arrays) as temp_ref_profile/humid_ref_profile only hold the *uncorrected* adiabat; the
    # correction is masked to the in-column range so levels above LZB stay exactly at the
    # environment and contribute exact zero, full-range and branchless
    τ⁻¹ = inv(convert(NF, Second(convection.time_scale).value))
    for k in 1:nlayers
        in_column = k >= level_zero_buoyancy
        T_ref = temp_ref_profile[ij, k] - ifelse(in_column, ΔT, zero(NF))    # eq (6) / (15)
        q_ref = humid_ref_profile[ij, k] * ifelse(in_column, fq, one(NF))    # eq (13)

        temp_tend[ij, k] -= (temp[ij, k] - T_ref) * τ⁻¹
        δq = (humid[ij, k] - q_ref) * τ⁻¹
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

    # CONVECTIVE SNOW: below freezing, all convective rain falls as snow instead. Convective
    # precipitation is deposited immediately (not fluxed through layers like large-scale
    # precipitation), so a single check on the lowermost layer's temperature is sufficient,
    # unlike the falling-flux melt cascade in large-scale condensation
    snow_convection::NF = 0
    freezing = convection.snow && temp[ij, nlayers] < convection.freezing_threshold
    rain_convection, snow_convection = ifelse(
        freezing,
        (snow_convection, rain_convection),
        (rain_convection, snow_convection),
    )

    # Store precipitation in diagnostic arrays
    vars.parameterizations.rain_convection[ij] += rain_convection            # accumulated rain [m] for output
    vars.parameterizations.snow_convection[ij] += snow_convection            # accumulated snow [m] for output
    rain_rate_convection = rain_convection / Δt                     # instantaneous rate [m/s] for coupling
    snow_rate_convection = snow_convection / Δt                     # instantaneous rate [m/s] for coupling
    vars.parameterizations.rain_rate_convection[ij] = rain_rate_convection   # instantaneous rate [m/s] for coupling
    vars.parameterizations.snow_rate_convection[ij] = snow_rate_convection   # instantaneous rate [m/s] for coupling

    # accumulate into total rain/snow rate including large-scale condensation [m/s]
    vars.parameterizations.rain_rate[ij] += rain_rate_convection             # instantaneous rate [m/s] for coupling
    vars.parameterizations.snow_rate[ij] += snow_rate_convection             # instantaneous rate [m/s] for coupling

    # clouds reach to top of convection
    vars.parameterizations.cloud_top[ij] = min(vars.parameterizations.cloud_top[ij], level_zero_buoyancy)
    return nothing
end

function variables(::BettsMillerConvection)
    return (
        ParameterizationVariable(:rain_convection, Grid2D(), desc = "Convective precipitation (accumulated)", units = "m"),
        ParameterizationVariable(:rain_rate_convection, Grid2D(), desc = "Convective precipitation rate", units = "m/s"),
        ParameterizationVariable(:snow_convection, Grid2D(), desc = "Convective precipitation (snow, accumulated)", units = "m"),
        ParameterizationVariable(:snow_rate_convection, Grid2D(), desc = "Convective snow rate", units = "m/s"),
        ParameterizationVariable(:rain_rate, Grid2D(), desc = "Rain rate (large-scale + convection)", units = "m/s"),
        ParameterizationVariable(:snow_rate, Grid2D(), desc = "Snow rate (large-scale + convection)", units = "m/s"),
        ParameterizationVariable(:cloud_top, Grid2D(), desc = "Cloud top layer index", units = "1"),
    )
end

"""
$(TYPEDSIGNATURES)
Calculates the moist pseudo adiabat given temperature and humidity of surface parcel, and the
associated reference humidity profile ``q_{ref} = RH \\, q^\\star(T_{ref})``. Follows the dry
adiabat till condensation and then continues on the pseudo moist-adiabat with immediate
condensation to the level of zero buoyancy (LZB). Levels above the LZB are set to the
*environmental* temperature/humidity (rather than left undefined), so that
``T - T_{ref} = 0`` and ``q - q_{ref} = 0`` there exactly: downstream loops over the full
column then need no branch to skip those levels, they contribute zero on their own."""
@propagate_inbounds function pseudo_adiabat!(
        ij,
        temp_ref_profile,
        humid_ref_profile,
        temp_environment,
        humid_environment,
        geopotential,
        pres,
        σ,
        atmosphere,
        relative_humidity,
        entrainment,
    )
    NF = eltype(temp_ref_profile)             # number format
    nlayers = length(σ)                       # number of vertical layers

    # thermodynamics
    (; κ, R_dry, R_vapor) = atmosphere
    Lᵥ = atmosphere.latent_heat_condensation    # latent heat of vaporization
    cₚ = atmosphere.heat_capacity               # heat capacity

    # levels above the LZB are never touched again below, so they stay at the environment
    for k in 1:nlayers
        temp_ref_profile[ij, k] = temp_environment[ij, k]
        humid_ref_profile[ij, k] = humid_environment[ij, k]
    end

    temp_parcel::NF = temp_environment[ij, nlayers]     # start at surface with environment temperature [K]
    humid_parcel::NF = humid_environment[ij, nlayers]   # and humidity [kg/kg]
    temp_ref_profile[ij, nlayers] = temp_parcel     # start profile at surface with parcel temperature
    sat_humid = saturation_humidity(temp_parcel, σ[nlayers] * pres, atmosphere)
    humid_ref_profile[ij, nlayers] = relative_humidity * sat_humid

    saturated::Bool = false           # did the parcel reach saturation yet?
    buoyant::Bool = true              # is the parcel still buoyant?
    k::Int = nlayers                  # layer index top to surface
    temp_virt_parcel::NF = virtual_temperature(temp_parcel, humid_parcel, atmosphere)

    while buoyant && k > 1                  # calculate moist adiabat while buoyant till top
        k -= 1                              # one level up

        if !saturated                       # if not saturated yet follow dry adiabat
            # dry adiabatic ascent and saturation humidity of that temperature
            temp_parcel_dry = temp_parcel * (σ[k] / σ[k + 1])^κ
            sat_humid = saturation_humidity(temp_parcel_dry, σ[k] * pres, atmosphere)

            # set to saturated when the dry adiabatic ascent would reach saturation
            # then follow moist adiabat instead (see below)
            saturated = humid_parcel >= sat_humid
        end

        if saturated
            # calculate moist/pseudo adiabatic lapse rate, dT/dΦ = -Γ/cp
            T, Tᵥ, q = temp_parcel, temp_virt_parcel, humid_parcel  # for brevity
            A = q * Lᵥ / ((1 - q)^2 * R_dry)
            B = q * Lᵥ^2 / ((1 - q)^2 * cₚ * R_vapor)
            Γ = (1 + A / Tᵥ) / (1 + B / T^2)

            ΔΦ = geopotential[ij, k] - geopotential[ij, k + 1]        # vertical gradient in geopotential
            temp_parcel = temp_parcel - ΔΦ / cₚ * Γ                     # new temperature of parcel at k

            # at new (lower) temperature condensation occurs immediately
            # new humidity equals to that saturation humidity
            humid_parcel = saturation_humidity(temp_parcel, σ[k] * pres, atmosphere)
            sat_humid = humid_parcel               # reused below for the reference humidity

            if entraining(entrainment)             # compiled away entirely for NoEntrainment
                # entrainment: mix the rising parcel with environmental air, diluting its
                # buoyancy and moisture; sat_humid is recomputed at the mixed temperature so
                # the reference humidity below stays consistent with the (now mixed) T_ref
                ε = entrainment(σ[k])
                temp_parcel = (1 - ε) * temp_parcel + ε * temp_environment[ij, k]
                humid_parcel = (1 - ε) * humid_parcel + ε * humid_environment[ij, k]
                sat_humid = saturation_humidity(temp_parcel, σ[k] * pres, atmosphere)
            end
        else
            temp_parcel = temp_parcel_dry       # else parcel temperature following dry adiabat
            # sat_humid already holds saturation_humidity(temp_parcel_dry, σ[k]*pres, atmosphere)
        end

        # use dry/moist adiabatic ascent for reference profile, and its saturation humidity
        # (already computed above in both branches) scaled by RH for the humidity reference
        temp_ref_profile[ij, k] = temp_parcel
        humid_ref_profile[ij, k] = relative_humidity * sat_humid

        # check whether parcel is still buoyant wrt to environment
        # use virtual temperature as it's equivalent to density
        temp_virt_parcel = virtual_temperature(temp_parcel, humid_parcel, atmosphere)         # virtual temperature of parcel
        buoyant = temp_virt_parcel > virtual_temperature(temp_environment[ij, k], humid_environment[ij, k], atmosphere)
    end

    # if parcel isn't buoyant anymore restore the environment at the level buoyancy was lost
    # (replaces the previous NaN marker)
    temp_ref_profile[ij, k] = ifelse(buoyant, temp_ref_profile[ij, k], temp_environment[ij, k])
    humid_ref_profile[ij, k] = ifelse(buoyant, humid_ref_profile[ij, k], humid_environment[ij, k])

    # level of zero buoyancy is reached when the loop stops, but in case it's at the top it's still buoyant
    level_zero_buoyancy = k + (1 - buoyant)
    return level_zero_buoyancy
end

export BettsMillerDryConvection

"""
The simplified Betts-Miller convection scheme from Frierson, 2007,
https://doi.org/10.1175/JAS3935.1 but with humidity set to zero.
Fields and options are
$(TYPEDFIELDS)"""
@parameterized @kwdef struct BettsMillerDryConvection{NF, Entrainment <: AbstractEntrainment} <: AbstractConvection
    "[OPTION] Relaxation time for profile adjustment"
    time_scale::Second = Hour(4)

    "[OPTION] Entrainment profile mixing environmental air into the rising parcel"
    @component entrainment::Entrainment = NoEntrainment()
end

Adapt.@adapt_structure BettsMillerDryConvection

# generator function
BettsMillerDryConvection(SG::SpectralGrid; entrainment = NoEntrainment(), kwargs...) =
    BettsMillerDryConvection{SG.NF, typeof(entrainment)}(; entrainment, kwargs...)
initialize!(::BettsMillerDryConvection, ::PrimitiveEquation) = nothing

# function barrier
@propagate_inbounds parameterization!(ij, vars, convection_scheme::BettsMillerDryConvection, model) =
    convection!(ij, vars, convection_scheme, model)

"""
$(TYPEDSIGNATURES)
calculates temperature tendency for the dry convection scheme following the
simplified Betts-Miller convection from Frierson 2007 but with zero humidity.
Starts with a first-guess relaxation to determine the convective criterion,
then adjusts the reference profiles
for thermodynamic consistency (e.g. in dry convection the humidity profile is non-precipitating),
and relaxes current vertical profiles to the adjusted references."""
@propagate_inbounds function convection!(ij, vars, DBM::BettsMillerDryConvection, model)

    (; geometry, atmosphere, time_stepping) = model
    # TODO: same as above — σ, σ_half, Δσ baked into Betts-Miller dry convection scheme.
    σ = geometry.σ_levels_full
    σ_half = geometry.σ_levels_half
    Δσ = geometry.σ_levels_thick
    nlayers = length(σ)

    # use previous time step for more stable calculations
    temp = get_prognostic_step(vars.grid.temperature, time_stepping, DBM)
    NF = eltype(temp)
    temp_tend = get_tendency_step(vars.tendencies.grid.temperature, time_stepping, DBM)

    # use work arrays for temp_ref_profile
    temp_ref_profile = vars.scratch.grid.a     # temperature [K] reference profile to adjust to

    # CONVECTIVE CRITERIA AND FIRST GUESS RELAXATION
    # Use surface temperature directly (simplified for now)
    temp_parcel = temp[ij, nlayers]
    level_zero_buoyancy = dry_adiabat!(
        ij, temp_ref_profile,
        temp,
        temp_parcel,
        σ,
        atmosphere,
        DBM.entrainment,
    )

    PT::NF = 0        # precipitation due to cooling

    # skip constants compared to Frierson 2007, i.e. no /τ, /gravity, *cₚ/Lᵥ
    # full-range, branchless: above LZB, temp_ref_profile == temp (see dry_adiabat!),
    # contributing exact zero without needing to mask the loop range
    for k in 1:nlayers
        # Frierson's equation (1)
        # δT = -(temp[ij, k] - temp_ref_profile[ij, k])/DBM.time_scale.value
        # PT += δT*Δσ[k]/gravity*cₚ/Lᵥ

        # shorter form with same sign (τ, gravity, cₚ, Lᵥ all positive) to be reused
        PT -= (temp[ij, k] - temp_ref_profile[ij, k]) * Δσ[k]
    end

    # ADJUST PROFILES FOLLOWING FRIERSON 2007
    convection = PT > 0
    convection || return nothing            # escape immediately for no convection

    # height of zero buoyancy level in σ coordinates
    Δσ_lzb = σ_half[nlayers + 1] - σ_half[level_zero_buoyancy]
    ΔT = PT / Δσ_lzb                          # eq (5) or (14) but reusing PT

    # corrected reference profile computed inline (not written back), masked to the
    # in-column range, full-range and branchless
    for k in 1:nlayers
        T_ref = temp_ref_profile[ij, k] - ifelse(k >= level_zero_buoyancy, ΔT, zero(NF))  # eq (6) or (15)
        temp_tend[ij, k] -= (temp[ij, k] - T_ref) / DBM.time_scale.value
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
Calculates the dry adiabat given the temperature of the surface parcel. Follows the dry
adiabat to the level of zero buoyancy (LZB). Levels above the LZB are set to the
*environmental* temperature (rather than left undefined), so that ``T - T_{ref} = 0`` there
exactly: downstream loops over the full column then need no branch to skip those levels."""
@propagate_inbounds function dry_adiabat!(
        ij,
        temp_ref_profile,
        temp_environment,
        temp_parcel,
        σ,
        atmosphere,
        entrainment,
    )
    NF = eltype(temp_ref_profile)
    (; κ) = atmosphere

    nlayers = length(σ)                     # number of vertical levels

    # levels above the LZB are never touched again below, so they stay at the environment
    for k in 1:nlayers
        temp_ref_profile[ij, k] = temp_environment[ij, k]
    end
    temp_ref_profile[ij, nlayers] = temp_parcel    # start profile at surface with parcel temperature

    buoyant::Bool = true                    # is the parcel still buoyant?
    k::Int = nlayers                        # layer index top to surface

    while buoyant && k > 1                  # calculate moist adiabat while buoyant till top
        k -= 1                              # one level up

        # dry adiabatic ascent
        temp_parcel = temp_parcel * (σ[k] / σ[k + 1])^κ

        if entraining(entrainment)          # compiled away entirely for NoEntrainment
            # entrainment: mix the rising parcel with environmental air, diluting its buoyancy
            ε = entrainment(σ[k])
            temp_parcel = (1 - ε) * temp_parcel + ε * temp_environment[ij, k]
        end

        temp_ref_profile[ij, k] = temp_parcel

        # check whether parcel is still buoyant wrt to environment
        buoyant = temp_parcel > temp_environment[ij, k]
    end

    # if parcel isn't buoyant anymore restore the environment at the level buoyancy was lost
    # (replaces the previous NaN marker)
    temp_ref_profile[ij, k] = ifelse(buoyant, temp_ref_profile[ij, k], temp_environment[ij, k])

    # level of zero buoyancy is reached when the loop stops, but in case it's at the top it's still buoyant
    level_zero_buoyancy = k + (1 - buoyant)
    return level_zero_buoyancy
end

export ConvectiveHeating

"""Convective heating as defined by Lee and Kim, 2003, JAS
implemented as convection parameterization. Fields are
$(TYPEDFIELDS)"""
@parameterized @kwdef struct ConvectiveHeating{NF, VectorType} <: AbstractConvection
    "[OPTION] Q_max heating strength as 1K/time_scale"
    time_scale::Second = Hour(12)

    "[OPTION] Pressure of maximum heating [hPa]"
    @param p₀::NF = 525 (bounds = Positive,)

    "[OPTION] Vertical extent of heating [hPa]"
    @param σₚ::NF = 200 (bounds = Positive,)

    "[OPTION] Latitude of heating [˚N]"
    @param θ₀::NF = 0 (bounds = -90 .. 90,)

    "[OPTION] Latitudinal width of heating [˚]"
    @param σθ::NF = 20 (bounds = Positive,)

    "[DERIVED] Latitudinal mask"
    lat_mask::VectorType
end

Adapt.@adapt_structure ConvectiveHeating

# generator
ConvectiveHeating(SG::SpectralGrid; kwargs...) = ConvectiveHeating{SG.NF, SG.VectorType}(lat_mask = zeros(SG.nlat); kwargs...)

# precompute latitudinal mask
function initialize!(C::ConvectiveHeating, model::PrimitiveEquation)
    θ = model.geometry.latd
    (; θ₀, σθ) = C

    # Lee and Kim, 2003, eq. 2
    return @. C.lat_mask .= cos(deg2rad((θ - θ₀) / σθ))^2
end

# function barrier
@propagate_inbounds parameterization!(ij, vars, convection_scheme::ConvectiveHeating, model) =
    convection!(ij, vars, convection_scheme, model)

@propagate_inbounds function convection!(
        ij,
        vars,
        scheme::ConvectiveHeating,
        model,
    )
    time_stepping = model.time_stepping
    pₛ = vars.parameterizations.surface_pressure[ij]          # surface pressure [Pa]
    temp_tend = get_tendency_step(vars.tendencies.grid.temperature, time_stepping, scheme)
    nlayers = size(temp_tend, 2)
    NF = eltype(temp_tend)

    # Get latitude ring index and latitude
    j = model.geometry.whichring[ij]
    latd = model.geometry.latd[j]
    σ = model.geometry.σ_levels_full

    Qmax = ifelse(abs(latd) < scheme.σθ, inv(convert(NF, Second(scheme.time_scale).value)), NF(0))
    p₀ = scheme.p₀ * 100     # hPa -> Pa
    σₚ = scheme.σₚ * 100     # hPa -> Pa
    cos²θ_term = scheme.lat_mask[j]

    for k in 1:nlayers
        p = pₛ * σ[k]       # Pressure in Pa on layer k

        # Lee and Kim, 2003, eq. 2
        temp_tend[ij, k] += Qmax * exp(-((p - p₀) / σₚ)^2 / 2) * cos²θ_term
    end
    return nothing
end
