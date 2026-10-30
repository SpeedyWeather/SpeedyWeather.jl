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
    @component entrainment::Entrainment = LinearEntrainment(DEFAULT_NF)

    "[OPTION] Convert convective rain below freezing to snow?"
    snow::Bool = true

    "[OPTION] Temperature of the lowermost layer below which convective rain becomes snow [K]"
    @param freezing_threshold::NF = 273.15 (bounds = Positive,)
end

Adapt.@adapt_structure BettsMillerConvection

# generator function
BettsMillerConvection(SG::SpectralGrid; entrainment = LinearEntrainment(SG), kwargs...) =
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
    # Coordinate-agnostic: level detection and moisture/energy flux integrals go through the
    # `pressure`/`pressure_half`/`pressure_thickness_ratio` interface, so this reduces exactly
    # to the previous σ-coordinate formulation for `SigmaCoordinates`. The vertical placement
    # of the convective adjustment (i.e. which levels count as "the column") remains a
    # nominal-σ criterion (see `pseudo_adiabat!`/known limitations), only the pressures and
    # mass weights used once that column is fixed are now exact for hybrid coordinates.
    coordinate = geometry.vertical_coordinates
    nlayers = geometry.nlayers
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
        geopotential, pₛ, coordinate, atmosphere, convection.relative_humidity, convection.entrainment,
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
        # Pq -= δq*Δp_k/(pₛ*gravity)
        #
        # δT = -(temp[ij, k] - temp_ref_profile[ij, k])/SBM.time_scale.value
        # PT += δT*Δp_k/(pₛ*gravity)*cₚ/Lᵥ

        # shorter form with same sign (τ, gravity, cₚ, Lᵥ all positive) to be reused
        # Δp_k/pₛ (mass weight of layer k) replaces the σ-coordinate Δσ[k]
        δ_k = pressure_thickness_ratio(k, pₛ, coordinate)
        Pq += (humid[ij, k] - humid_ref_profile[ij, k]) * δ_k
        PT -= (temp[ij, k] - temp_ref_profile[ij, k]) * δ_k

        in_column = k >= level_zero_buoyancy
        Qref -= ifelse(in_column, humid_ref_profile[ij, k], zero(NF)) * δ_k  # eq (11), mass-weighted
    end

    # ADJUST PROFILES FOLLOWING FRIERSON 2007
    deep_convection = Pq > 0 && PT > 0
    shallow_convection = Pq <= 0 && PT > 0

    # escape immediately for no convection
    no_convection = !(deep_convection || shallow_convection)
    no_convection && return nothing

    # mass fraction of the column from the level of zero buoyancy to the surface,
    # Σ_{k=level_zero_buoyancy}^{nlayers} δ_k = (p_half[nlayers+1] - p_half[level_zero_buoyancy])/pₛ,
    # which reduces exactly to σ_half[nlayers+1] - σ_half[level_zero_buoyancy] for SigmaCoordinates
    Δσ_lzb = pressure_ratio_half(nlayers + 1, pₛ, coordinate) - pressure_ratio_half(level_zero_buoyancy, pₛ, coordinate)

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
    τ⁻¹ = 1 / convert(NF, Second(convection.time_scale).value)    # TODO: `inv` fails in kernels on Reactant+CUDA https://github.com/EnzymeAD/Reactant.jl/issues/2560
    for k in 1:nlayers
        in_column = k >= level_zero_buoyancy
        T_ref = temp_ref_profile[ij, k] - ifelse(in_column, ΔT, zero(NF))    # eq (6) / (15)
        q_ref = humid_ref_profile[ij, k] * ifelse(in_column, fq, one(NF))    # eq (13)

        temp_tend[ij, k] -= (temp[ij, k] - T_ref) * τ⁻¹
        δq = (humid[ij, k] - q_ref) * τ⁻¹
        humid_tend[ij, k] -= δq

        # convective precipitation (rain), integrate dq\dt [(kg/kg)/s] vertically *with sign*:
        # moistening levels (δq < 0) are fed by the drying ones, so rain is the net column drying
        # otherwise water would be counted twice, as vapour and as precipitation, see #1275
        rain_convection += δq * pressure_thickness_ratio(k, pₛ, coordinate)   # integrate vertically, Formula 25, unit [m]
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

    # clouds reach to top of (precipitating, i.e. deep) convection, as for large-scale condensation
    cloud_top_convection = ifelse(deep_convection, level_zero_buoyancy, nlayers + 1)
    vars.parameterizations.cloud_top[ij] = min(vars.parameterizations.cloud_top[ij], cloud_top_convection)
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
        coordinate,
        atmosphere,
        relative_humidity,
        entrainment,
    )
    NF = eltype(temp_ref_profile)             # number format
    nlayers = get_nlayers(coordinate)         # number of vertical layers

    # thermodynamics
    (; κ, R_dry, R_vapor) = atmosphere
    Lᵥ = atmosphere.latent_heat_condensation    # latent heat of vaporization
    cₚ = atmosphere.heat_capacity               # heat capacity

    temp_parcel::NF = temp_environment[ij, nlayers]     # start at surface with environment temperature [K]
    humid_parcel::NF = humid_environment[ij, nlayers]   # and humidity [kg/kg]
    temp_ref_profile[ij, nlayers] = temp_parcel     # start profile at surface with parcel temperature
    sat_humid = saturation_humidity(temp_parcel, pressure(nlayers, pres, coordinate), atmosphere)
    humid_ref_profile[ij, nlayers] = relative_humidity * sat_humid

    saturated::Bool = false           # did the parcel reach saturation yet?
    buoyant::Bool = true              # is the parcel still buoyant?
    temp_virt_parcel::NF = virtual_temperature(temp_parcel, humid_parcel, atmosphere)
    level_zero_buoyancy::Int = 1      # in case the parcel is buoyant till the top

    # single pass over the column with a fixed trip count (no while loop, see #1193): the parcel
    # ascent is only computed while buoyant, every level is written exactly once, either with
    # the ascent or with the environment above the LZB, so that T - T_ref = q - q_ref = 0 there
    for k in (nlayers - 1):-1:1             # one level up at a time, surface to top
        temp_env = temp_environment[ij, k]
        humid_env = humid_environment[ij, k]

        if buoyant                          # calculate moist adiabat while buoyant till top
            # defined on all paths: the compiler can't prove that the `else` below is only reached
            # via `!saturated` here, a possibly undefined variable leaves an UndefVarError check
            # in the kernel
            temp_parcel_dry = temp_parcel
            if !saturated                   # if not saturated yet follow dry adiabat
                # dry adiabatic ascent and saturation humidity of that temperature. The Poisson
                # relation needs the pressure ratio p_k/p_k+1, taken via pressure_ratio (= p/pₛ) so
                # that the pₛ cancels exactly rather than up to rounding, i.e. σ[k]/σ[k+1] for
                # SigmaCoordinates. Layer k+1 sits below (closer to the surface).
                temp_parcel_dry = temp_parcel * (pressure_ratio(k, pres, coordinate) / pressure_ratio(k + 1, pres, coordinate))^κ
                sat_humid = saturation_humidity(temp_parcel_dry, pressure(k, pres, coordinate), atmosphere)

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
                humid_parcel = saturation_humidity(temp_parcel, pressure(k, pres, coordinate), atmosphere)
                sat_humid = humid_parcel               # reused below for the reference humidity

                if entraining(entrainment)             # compiled away entirely for NoEntrainment
                    # entrainment: mix the rising parcel with environmental air, diluting its
                    # buoyancy and moisture; sat_humid is recomputed at the mixed temperature so
                    # the reference humidity below stays consistent with the (now mixed) T_ref
                    ε = entrainment(pressure_ratio(k, pres, coordinate))
                    temp_parcel = (1 - ε) * temp_parcel + ε * temp_env
                    humid_parcel = (1 - ε) * humid_parcel + ε * humid_env
                    sat_humid = saturation_humidity(temp_parcel, pressure(k, pres, coordinate), atmosphere)
                end
            else
                temp_parcel = temp_parcel_dry       # else parcel temperature following dry adiabat
                # sat_humid already holds saturation_humidity(temp_parcel_dry, pressure(k, pres, coordinate), atmosphere)
            end

            # check whether parcel is still buoyant wrt to environment
            # use virtual temperature as it's equivalent to density
            temp_virt_parcel = virtual_temperature(temp_parcel, humid_parcel, atmosphere)
            buoyant = temp_virt_parcel > virtual_temperature(temp_env, humid_env, atmosphere)

            # the LZB is the last buoyant level, i.e. one below the level buoyancy was lost
            level_zero_buoyancy = ifelse(buoyant, k, k + 1)
        end

        # use dry/moist adiabatic ascent for reference profile, and its saturation humidity
        # (already computed above in both branches) scaled by RH for the humidity reference;
        # the environment at and above the level buoyancy was lost (replaces the previous NaN marker)
        temp_ref_profile[ij, k] = ifelse(buoyant, temp_parcel, temp_env)
        humid_ref_profile[ij, k] = ifelse(buoyant, relative_humidity * sat_humid, humid_env)
    end

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

# NF isn't used by any field so it can't be inferred by Adapt.@adapt_structure, pass it explicitly
function Adapt.adapt_structure(to, convection::BettsMillerDryConvection{NF}) where {NF}
    entrainment = adapt_structure(to, convection.entrainment)
    return BettsMillerDryConvection{NF, typeof(entrainment)}(convection.time_scale, entrainment)
end

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
    # Coordinate-agnostic: as in the wet Betts-Miller scheme above, the flux integral and
    # zero-buoyancy layer mass use `pressure_thickness_ratio`/`pressure_half`, and the dry
    # adiabat in `dry_adiabat!` uses actual pressure ratios, so this reduces exactly to the
    # previous σ-coordinate formulation for `SigmaCoordinates`.
    coordinate = geometry.vertical_coordinates
    nlayers = geometry.nlayers
    pₛ = vars.parameterizations.surface_pressure[ij]        # surface pressure [Pa]

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
        pₛ,
        coordinate,
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
        # PT += δT*Δp_k/(pₛ*gravity)*cₚ/Lᵥ

        # shorter form with same sign (τ, gravity, cₚ, Lᵥ all positive) to be reused
        PT -= (temp[ij, k] - temp_ref_profile[ij, k]) * pressure_thickness_ratio(k, pₛ, coordinate)
    end

    # ADJUST PROFILES FOLLOWING FRIERSON 2007
    convection = PT > 0
    convection || return nothing            # escape immediately for no convection

    # mass fraction of the column from the level of zero buoyancy to the surface (see
    # BettsMillerConvection above for the derivation); replaces the σ-coordinate
    # σ_half[nlayers+1] - σ_half[level_zero_buoyancy]
    Δσ_lzb = pressure_ratio_half(nlayers + 1, pₛ, coordinate) - pressure_ratio_half(level_zero_buoyancy, pₛ, coordinate)
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
        pres,
        coordinate,
        atmosphere,
        entrainment,
    )
    NF = eltype(temp_ref_profile)
    (; κ) = atmosphere

    nlayers = get_nlayers(coordinate)       # number of vertical levels

    temp_ref_profile[ij, nlayers] = temp_parcel    # start profile at surface with parcel temperature

    buoyant::Bool = true                    # is the parcel still buoyant?
    level_zero_buoyancy::Int = 1            # in case the parcel is buoyant till the top

    # single pass over the column with a fixed trip count (no while loop, see #1193), every level
    # is written exactly once, with the environment above the LZB so that T - T_ref = 0 there
    for k in (nlayers - 1):-1:1             # one level up at a time, surface to top
        temp_env = temp_environment[ij, k]

        if buoyant                          # calculate dry adiabat while buoyant till top
            # dry adiabatic ascent. The Poisson relation needs the pressure ratio p_k/p_k+1, taken
            # via pressure_ratio (= p/pₛ) so that the pₛ cancels exactly rather than up to rounding,
            # i.e. σ[k]/σ[k+1] for SigmaCoordinates. Layer k+1 sits below (closer to the surface).
            temp_parcel = temp_parcel * (pressure_ratio(k, pres, coordinate) / pressure_ratio(k + 1, pres, coordinate))^κ

            if entraining(entrainment)      # compiled away entirely for NoEntrainment
                # entrainment: mix the rising parcel with environmental air, diluting its buoyancy
                ε = entrainment(pressure_ratio(k, pres, coordinate))
                temp_parcel = (1 - ε) * temp_parcel + ε * temp_env
            end

            # check whether parcel is still buoyant wrt to environment
            buoyant = temp_parcel > temp_env

            # the LZB is the last buoyant level, i.e. one below the level buoyancy was lost
            level_zero_buoyancy = ifelse(buoyant, k, k + 1)
        end

        # environment at and above the level buoyancy was lost (replaces the previous NaN marker)
        temp_ref_profile[ij, k] = ifelse(buoyant, temp_parcel, temp_env)
    end

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
    coordinate = model.geometry.vertical_coordinates

    Qmax = ifelse(abs(latd) < scheme.σθ, 1 / convert(NF, Second(scheme.time_scale).value), NF(0))    # TODO: `inv` fails in kernels on Reactant+CUDA https://github.com/EnzymeAD/Reactant.jl/issues/2560
    p₀ = scheme.p₀ * 100     # hPa -> Pa
    σₚ = scheme.σₚ * 100     # hPa -> Pa
    cos²θ_term = scheme.lat_mask[j]

    for k in 1:nlayers
        p = pressure(k, pₛ, coordinate)       # Pressure in Pa on layer k

        # Lee and Kim, 2003, eq. 2
        temp_tend[ij, k] += Qmax * exp(-((p - p₀) / σₚ)^2 / 2) * cos²θ_term
    end
    return nothing
end
