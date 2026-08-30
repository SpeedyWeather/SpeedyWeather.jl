abstract type AbstractConvection <: AbstractParameterization end

export BettsMillerConvection

"""The simplified Betts-Miller convection scheme from Frierson, 2007,
https://doi.org/10.1175/JAS3935.1. This implements the qref-formulation
in their paper. Fields and options are $(TYPEDFIELDS)"""
@parameterized @kwdef struct BettsMillerConvection{NF} <: AbstractConvection
    "[OPTION] Relaxation time for profile adjustment"
    time_scale::Second = Hour(4)

    "[OPTION] Relative humidity for reference profile [1]"
    @param relative_humidity::NF = 0.7
end

Adapt.@adapt_structure BettsMillerConvection

# generator function
BettsMillerConvection(SG::SpectralGrid; kwargs...) = BettsMillerConvection{SG.NF}(; kwargs...)
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
    level_zero_buoyancy = pseudo_adiabat!(ij, temp_ref_profile, temp, humid, geopotential, pₛ, coordinate, atmosphere)

    for k in level_zero_buoyancy:nlayers
        qsat = saturation_humidity(temp_ref_profile[ij, k], pressure(k, pₛ, coordinate), atmosphere)
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
        # Pq -= δq*Δp_k/(pₛ*gravity)
        #
        # δT = -(temp[ij, k] - temp_ref_profile[ij, k])/SBM.time_scale.value
        # PT += δT*Δp_k/(pₛ*gravity)*cₚ/Lᵥ

        # shorter form with same sign (τ, gravity, cₚ, Lᵥ all positive) to be reused
        # Δp_k/pₛ (mass weight of layer k) replaces the σ-coordinate Δσ[k]
        δ_k = pressure_thickness_ratio(k, pₛ, coordinate)
        Pq += (humid[ij, k] - humid_ref_profile[ij, k]) * δ_k
        PT -= (temp[ij, k] - temp_ref_profile[ij, k]) * δ_k
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

    if deep_convection

        ΔT = (PT - Pq * Lᵥ / cₚ) / Δσ_lzb         # eq (5) but reusing PT, Pq, and /cₚ already included

        for k in level_zero_buoyancy:nlayers
            temp_ref_profile[ij, k] -= ΔT   # equation (6)
        end

    elseif shallow_convection

        # FRIERSON'S QREF SCHEME
        # changing the reference profiles for both temperature and humidity so the
        # precipitation is zero.

        for k in level_zero_buoyancy:nlayers
            Qref -= humid_ref_profile[ij, k] * pressure_thickness_ratio(k, pₛ, coordinate)  # eq (11), mass-weighted
        end
        fq = 1 - Pq / Qref                    # = 1 - Δq/Qref in eq (12) but we reuse Pq

        ΔT = PT / Δσ_lzb                      # equation (14), reuse PT, mass-weighted (δ_k rather than Δσ_k)
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
        rain = max(δq * pressure_thickness_ratio(k, pₛ, coordinate), 0)   # only integrate excess humidity for precip (no reevaporation)
        rain_convection += rain         # integrate vertically, Formula 25, unit [m]
    end

    # CONVECTIVE PRECIPITATION
    # enforce no precip for shallow conv via * deep_convection flag
    pₛΔt_gρ = (pₛ * Δt / (g * ρ)) * deep_convection
    rain_convection *= pₛΔt_gρ                  # convert to [m] of rain during Δt
    rain_convection = max(rain_convection, 0)   # ensure non-negative precipitation, rounding errors

    # Store precipitation in diagnostic arrays
    vars.parameterizations.rain_convection[ij] += rain_convection            # accumulated rain [m] for output
    rain_rate_convection = rain_convection / Δt                     # instantaneous rate [m/s] for coupling
    vars.parameterizations.rain_rate_convection[ij] = rain_rate_convection   # instantaneous rate [m/s] for coupling

    # accumulate into total rain rate including large-scale condensation [m/s]
    vars.parameterizations.rain_rate[ij] += rain_rate_convection             # instantaneous rate [m/s] for coupling

    # clouds reach to top of convection
    vars.parameterizations.cloud_top[ij] = min(vars.parameterizations.cloud_top[ij], level_zero_buoyancy)
    return nothing
end

function variables(::BettsMillerConvection)
    return (
        ParameterizationVariable(:rain_convection, Grid2D(), desc = "Convective precipitation (accumulated)", units = "m"),
        ParameterizationVariable(:rain_rate_convection, Grid2D(), desc = "Convective precipitation rate", units = "m/s"),
        ParameterizationVariable(:rain_rate, Grid2D(), desc = "Rain rate (large-scale + convection)", units = "m/s"),
        ParameterizationVariable(:cloud_top, Grid2D(), desc = "Cloud top layer index", units = "1"),
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
        coordinate,
        atmosphere
    )
    NF = eltype(temp_ref_profile)             # number format
    nlayers = get_nlayers(coordinate)         # number of vertical layers

    # thermodynamics
    (; κ, R_dry, R_vapor) = atmosphere
    Lᵥ = atmosphere.latent_heat_condensation    # latent heat of vaporization
    cₚ = atmosphere.heat_capacity               # heat capacity

    temp_parcel::NF = temp_environment[ij, nlayers]     # start at surface with environment temperature [K]
    humid_parcel::NF = humid_environment[ij, nlayers]   # and humidity [kg/kg]

    for k in 1:nlayers
        temp_ref_profile[ij, k] = NaN           # reset profile from any previous calculation, TODO necessary?
    end
    temp_ref_profile[ij, nlayers] = temp_parcel     # start profile at surface with parcel temperature

    local saturated::Bool = false           # did the parcel reach saturation yet?
    local buoyant::Bool = true              # is the parcel still buoyant?
    local k::Int = nlayers                  # layer index top to surface
    local temp_virt_parcel::NF = virtual_temperature(temp_parcel, humid_parcel, atmosphere)

    while buoyant && k > 1                  # calculate moist adiabat while buoyant till top
        k -= 1                              # one level up

        if !saturated                       # if not saturated yet follow dry adiabat
            # dry adiabatic ascent and saturation humidity of that temperature. The Poisson
            # relation needs the pressure ratio p_k/p_k+1, taken via pressure_ratio (= p/pₛ) so
            # that the pₛ cancels exactly rather than up to rounding, i.e. σ[k]/σ[k+1] for
            # SigmaCoordinates. Layer k+1 sits below (closer to the surface).
            p_k = pressure(k, pres, coordinate)
            pressure_ratio_k = pressure_ratio(k, pres, coordinate)
            pressure_ratio_below = pressure_ratio(k + 1, pres, coordinate)
            temp_parcel_dry = temp_parcel * (pressure_ratio_k / pressure_ratio_below)^κ
            sat_humid = saturation_humidity(temp_parcel_dry, p_k, atmosphere)

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
        else
            temp_parcel = temp_parcel_dry       # else parcel temperature following dry adiabat
        end

        # use dry/moist adiabatic ascent for reference profile
        temp_ref_profile[ij, k] = temp_parcel

        # check whether parcel is still buoyant wrt to environment
        # use virtual temperature as it's equivalent to density
        temp_virt_parcel = virtual_temperature(temp_parcel, humid_parcel, atmosphere)         # virtual temperature of parcel
        buoyant = temp_virt_parcel > virtual_temperature(temp_environment[ij, k], humid_environment[ij, k], atmosphere)
    end

    # if parcel isn't buoyant anymore set last temperature (with negative buoyancy) back to NaN
    temp_ref_profile[ij, k] = !buoyant ? NaN32 : temp_ref_profile[ij, k]

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
@kwdef struct BettsMillerDryConvection{NF} <: AbstractConvection
    "[OPTION] Relaxation time for profile adjustment"
    time_scale::Second = Hour(4)
end

Adapt.adapt_structure(to, bmdc::BettsMillerDryConvection{NF}) where {NF} = BettsMillerDryConvection{NF}(adapt_structure(to, bmdc.time_scale))

# generator function
BettsMillerDryConvection(SG::SpectralGrid; kwargs...) = BettsMillerDryConvection{SG.NF}(; kwargs...)
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
        atmosphere
    )

    PT::NF = 0        # precipitation due to cooling
    ΔT::NF = 0        # vertically uniform temperature profile adjustment

    # skip constants compared to Frierson 2007, i.e. no /τ, /gravity, *cₚ/Lᵥ
    for k in level_zero_buoyancy:nlayers
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
        pres,
        coordinate,
        atmosphere,
    )
    NF = eltype(temp_ref_profile)
    (; κ) = atmosphere

    nlayers = get_nlayers(coordinate)       # number of vertical levels

    for k in 1:nlayers
        temp_ref_profile[ij, k] = NaN       # reset profile from any previous calculation
    end
    temp_ref_profile[ij, nlayers] = temp_parcel    # start profile at surface with parcel temperature

    buoyant::Bool = true                    # is the parcel still buoyant?
    k::Int = nlayers                        # layer index top to surface

    while buoyant && k > 1                  # calculate moist adiabat while buoyant till top
        k -= 1                              # one level up

        # dry adiabatic ascent. The Poisson relation needs the pressure ratio p_k/p_k+1, taken
        # via pressure_ratio (= p/pₛ) so that the pₛ cancels exactly rather than up to rounding,
        # i.e. σ[k]/σ[k+1] for SigmaCoordinates. Layer k+1 sits below (closer to the surface).
        pressure_ratio_k = pressure_ratio(k, pres, coordinate)
        pressure_ratio_below = pressure_ratio(k + 1, pres, coordinate)
        temp_parcel = temp_parcel * (pressure_ratio_k / pressure_ratio_below)^κ
        temp_ref_profile[ij, k] = temp_parcel

        # check whether parcel is still buoyant wrt to environment
        buoyant = temp_parcel > temp_environment[ij, k]
    end

    # if parcel isn't buoyant anymore set last temperature (with negative buoyancy) back to NaN
    temp_ref_profile[ij, k] = !buoyant ? NF(NaN) : temp_ref_profile[ij, k]

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
    coordinate = model.geometry.vertical_coordinates

    Qmax = ifelse(abs(latd) < scheme.σθ, inv(convert(NF, Second(scheme.time_scale).value)), NF(0))
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
