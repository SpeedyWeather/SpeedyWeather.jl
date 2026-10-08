export PrognosticCloudCondensation

"""Large-scale condensation with a prognostic cloud condensate (liquid + ice in one variable,
the phase decided by temperature), following the structure of Zhao and Carr (1997) with
Sundqvist et al. (1989) autoconversion. How much vapour condenses into the condensate, or
condensate evaporates, is decided by the `closure`: the implicit relaxation of
[`ImplicitCondensation`](@ref) ([`RelaxationClosure`](@ref), default) or the Sundqvist closure
([`SundqvistClosure`](@ref)). The condensate is advected and diffused
as a prognostic variable `cloud_condensate` and converted to rain and snow by autoconversion.
Rain and snow are diagnosed in one top-down sweep per time step (no storage, no fall speed),
with melting and reevaporation as in `ImplicitCondensation`. The scheme also writes the cloud
state read by radiation: layer cloud fraction (Xu and Randall 1996 in the form of NCEP's GFS),
grid-mean cloud liquid and ice water and their effective radii. Fields are
$(TYPEDFIELDS)"""
@parameterized @kwdef struct PrognosticCloudCondensation{NF, C} <: AbstractCondensation
    "[OPTION] Condensation closure: how much vapour condenses into the condensate, or evaporates from it"
    @component closure::C

    "[OPTION] Temperature [K] below which all condensate is ice, linear ramp to liquid at freezing"
    @param ice_temperature::NF = 253.15 (bounds = Positive,)

    "[OPTION] Saturation over ice for the ice fraction of the condensate (otherwise over liquid only)?"
    mixed_phase_saturation::Bool = true

    "[OPTION] Autoconversion rate of cloud liquid to rain (Sundqvist C₀) [1/s]"
    @param autoconversion_rate::NF = 1.0e-4 (bounds = Nonnegative,)

    "[OPTION] Characteristic in-cloud liquid water of the autoconversion (Sundqvist m_r) [kg/kg]"
    @param autoconversion_water::NF = 3.0e-4 (bounds = Positive,)

    "[OPTION] Enhancement of autoconversion by collection of precipitation from above (Sundqvist c₁) [(kg/m²/s)^-1/2]"
    @param collection_enhancement::NF = 300 (bounds = Nonnegative,)

    "[OPTION] Enhancement of autoconversion in supercooled cloud, Bergeron-Findeisen (Sundqvist c₂) [K^-1/2]"
    @param bergeron_findeisen_enhancement::NF = 0.5 (bounds = Nonnegative,)

    "[OPTION] Temperature [K] below which the Bergeron-Findeisen enhancement starts"
    @param bergeron_findeisen_temperature::NF = 268 (bounds = Positive,)

    "[OPTION] Autoconversion rate of cloud ice to snow at freezing, decays with exp(0.025(T-T₀)) [1/s]"
    @param ice_autoconversion_rate::NF = 3.0e-4 (bounds = Nonnegative,)

    "[OPTION] Condensate that is not converted to precipitation, at 1000 hPa, scaled with pressure [kg/kg]"
    @param minimum_condensate::NF = 1.0e-5 (bounds = Nonnegative,)

    "[OPTION] Coefficient of the Xu-Randall cloud fraction (GFS value) [1]"
    @param cloud_fraction_coefficient::NF = 2000 (bounds = Positive,)

    "[OPTION] Cloud fractions below this are set to zero [1]"
    @param cloud_fraction_min::NF = 0.001 (bounds = 0 .. 1,)

    "[OPTION] Reevaporation efficiency of rain [1/(kg/kg)], 0 for no reevaporation"
    @param reevaporation::NF = 30 (bounds = Nonnegative,)

    "[OPTION] Temperature [K] above which snow melts"
    @param melting_threshold::NF = 278 (bounds = Positive,)

    "[OPTION] Effective radius of cloud liquid over ocean [m]"
    @param liquid_effective_radius_ocean::NF = 10.0e-6 (bounds = Positive,)

    "[OPTION] Effective radius of warm cloud liquid over land, increasing to the ocean value with the ice fraction [m]"
    @param liquid_effective_radius_land::NF = 5.0e-6 (bounds = Positive,)

    "[OPTION] Effective radius of cloud ice [m]"
    @param ice_effective_radius::NF = 50.0e-6 (bounds = Positive,)
end

Adapt.@adapt_structure PrognosticCloudCondensation
PrognosticCloudCondensation(SG::SpectralGrid; closure = RelaxationClosure(SG), kwargs...) =
    PrognosticCloudCondensation{SG.NF, typeof(closure)}(; closure, kwargs...)
initialize!(::PrognosticCloudCondensation, ::PrimitiveEquation) = nothing

function variables(scheme::PrognosticCloudCondensation, model::AbstractModel)
    return (
        cloud_condensate_variables(get_nsteps(model.time_stepping, model))...,
        cloud_state_variables()...,
        large_scale_precipitation_variables()...,
        variables(scheme.closure, model)...,
    )
end

abstract type AbstractCondensationClosure <: AbstractModelComponent end

export RelaxationClosure

"""Condensation closure of [`PrognosticCloudCondensation`](@ref) by implicit relaxation as in
[`ImplicitCondensation`](@ref): humidity above `relative_humidity_threshold` condenses into the
condensate over `time_scale` time steps, with the implicit latent heat factor of Frierson et al.
(2006); below that threshold condensate evaporates in the clear part of the cell (Xu-Randall cloud
fraction) over `evaporation_time_scale` time steps. The cloud fraction for autoconversion is the
Xu-Randall cloud fraction after condensation. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct RelaxationClosure{NF} <: AbstractCondensationClosure
    "[OPTION] Relative humidity threshold [1 = 100%] above which humidity condenses"
    @param relative_humidity_threshold::NF = 0.95 (bounds = 0 .. 1,)

    "[OPTION] Time scale of condensation in multiples of the time step Δt"
    @param time_scale::NF = 3 (bounds = Positive,)

    "[OPTION] Time scale of cloud evaporation in the clear part of a cell in multiples of Δt"
    @param evaporation_time_scale::NF = 6 (bounds = Positive,)
end

Adapt.@adapt_structure RelaxationClosure
RelaxationClosure(SG::SpectralGrid; kwargs...) = RelaxationClosure{SG.NF}(; kwargs...)
initialize!(::RelaxationClosure, ::AbstractModel) = nothing

export SundqvistClosure

"""Condensation closure of [`PrognosticCloudCondensation`](@ref) after Sundqvist et al. (1989) in the
form of NCEP's GFS (Zhao and Carr 1997): above a critical relative humidity `u` a cell is partly
cloudy with the cloud fraction `b = 1 - √((1 - f)/(1 - u))` of its relative humidity `f`. In a cloudy
cell the part `β` of the moisture supply `M` by all other processes condenses,

    C = β M / (1 + L/cₚ f ∂q*/∂T),   M = A_q - f ∂q*/∂T A_T - f ∂q*/∂p A_p,

`β` between `b` (no condensate) and 1 (much condensate), bounded so that the cell is not dried below
`u`. In a clear cell (`b` below `cloud_fraction_threshold`) condensate evaporates towards `u` over
`evaporation_time_scale` time steps. The supply `A_X` of temperature, humidity and pressure is the
change of the lagged state since the state the scheme left two calls before (same leapfrog parity),
kept as reference variables in the `clouds` namespace of the prognostic variables. Before two calls
have stored a reference the supply is zero. The cloud fraction for autoconversion is `b`.
The critical relative humidity varies with pressure as in ECHAM,
`u = u_top + (u_surface - u_top) exp(1 - (pₛ/p)^n)`. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct SundqvistClosure{NF} <: AbstractCondensationClosure
    "[OPTION] Critical relative humidity at the surface [1]"
    @param critical_relative_humidity_surface::NF = 0.9 (bounds = 0 .. 1,)

    "[OPTION] Critical relative humidity at the top of the atmosphere [1]"
    @param critical_relative_humidity_top::NF = 0.9 (bounds = 0 .. 1,)

    "[OPTION] Exponent n of the pressure dependence of the critical relative humidity [1]"
    @param critical_relative_humidity_exponent::NF = 2 (bounds = Positive,)

    "[OPTION] Cloud fraction below which a cell counts as clear and its condensate evaporates [1]"
    @param cloud_fraction_threshold::NF = 1.0e-3 (bounds = 0 .. 1,)

    "[OPTION] Time scale of condensate evaporation in clear cells in multiples of Δt"
    @param evaporation_time_scale::NF = 2 (bounds = Positive,)
end

Adapt.@adapt_structure SundqvistClosure
SundqvistClosure(SG::SpectralGrid; kwargs...) = SundqvistClosure{SG.NF}(; kwargs...)
initialize!(::SundqvistClosure, ::AbstractModel) = nothing

# reference state of the cloud scheme, two copies for the two leapfrog parities:
# step 1 = from two calls before (used for the supply), step 2 = from the last call
variables(::SundqvistClosure) = (
    PrognosticVariable(:temperature_reference, GridXYZT(2), desc = "Temperature after the cloud scheme's condensation, two calls before (1) and last call (2)", units = "K", namespace = :clouds),
    PrognosticVariable(:humidity_reference, GridXYZT(2), desc = "Humidity after the cloud scheme's condensation, two calls before (1) and last call (2)", units = "kg/kg", namespace = :clouds),
    PrognosticVariable(:surface_pressure_reference, GridXYT(2), desc = "Surface pressure at the cloud scheme's call two calls before (1) and last call (2)", units = "Pa", namespace = :clouds),
)

"""$(TYPEDSIGNATURES)
Critical relative humidity at pressure `p` [Pa] and surface pressure `pₛ` [Pa] of the
[`SundqvistClosure`](@ref): `u_top + (u_surface - u_top) exp(1 - (pₛ/p)^n)`."""
@inline function critical_relative_humidity(closure::SundqvistClosure, p, pₛ)
    u_surface = closure.critical_relative_humidity_surface
    u_top = closure.critical_relative_humidity_top
    n = closure.critical_relative_humidity_exponent
    return u_top + (u_surface - u_top) * exp(1 - (pₛ / p)^n)
end

"""$(TYPEDSIGNATURES)
Sundqvist cloud fraction `b = 1 - √((1 - f)/(1 - u))` for relative humidity `f` above the critical
relative humidity `u`, zero below. `f` is clamped below 1 for a finite derivative."""
@inline function sundqvist_cloud_fraction(f, u)
    NF = typeof(f)
    f_clamped = clamp(f, u, 1 - NF(1.0e-6))
    return clamp(1 - sqrt((1 - f_clamped) / max(1 - u, NF(1.0e-6))), zero(NF), one(NF))
end

"""$(TYPEDSIGNATURES)
Condensation [kg/kg/s] of vapour into the condensate, evaporation [kg/kg/s] of condensate and the
closure's cloud fraction for layer `k` of column `ij` with the state in `layer`, see the closures."""
condensation_closure

@inline function condensation_closure(closure::RelaxationClosure, scheme, ij, k, vars, layer)
    (; q, qc, p, q_sat, dq_sat_dT, L, cₚ, Δt, Δt_prognostic) = layer
    NF = typeof(q)
    (; relative_humidity_threshold) = closure
    implicit_factor = 1 + L / cₚ * relative_humidity_threshold * dq_sat_dT
    excess = q - relative_humidity_threshold * q_sat
    condensation = max(excess, zero(NF)) / (implicit_factor * closure.time_scale * Δt)
    cloud_fraction = xu_randall_cloud_fraction(qc, q / q_sat, q_sat, p, scheme.cloud_fraction_coefficient, scheme.cloud_fraction_min)
    evaporation = (1 - cloud_fraction) * max(-excess, zero(NF)) / (implicit_factor * closure.evaporation_time_scale * Δt)
    evaporation = min(evaporation, qc / Δt_prognostic)     # not more than available
    return condensation, evaporation, cloud_fraction
end

@inline function condensation_closure(closure::SundqvistClosure, scheme, ij, k, vars, layer)
    (; T, q_lagged, q, qc, p, pₛ, pₛ_reference, q_sat, dq_sat_dT, L, cₚ, Δt, Δt_prognostic, coord) = layer
    NF = typeof(q)
    (; temperature_reference, humidity_reference) = vars.prognostic.clouds

    u = critical_relative_humidity(closure, p, pₛ)
    f = q / q_sat
    b = sundqvist_cloud_fraction(f, u)

    # moisture supply by all other processes: change of the lagged state since the state the
    # scheme left two calls before (same leapfrog parity), zero before a reference exists
    has_reference = pₛ_reference > 0
    A_T = ifelse(has_reference, (T - temperature_reference[ij, k, 1]) / Δt_prognostic, zero(NF))
    A_q = ifelse(has_reference, (q_lagged - humidity_reference[ij, k, 1]) / Δt_prognostic, zero(NF))
    p_reference = pressure(k, ifelse(has_reference, pₛ_reference, pₛ), coord)
    A_p = (p - p_reference) / Δt_prognostic
    dq_sat_dp = -q_sat / p
    M = A_q - f * dq_sat_dT * A_T - f * dq_sat_dp * A_p

    # cloudy: the part β of the supply condenses, β = b without condensate, → 1 with much condensate
    a = b * (1 - b) * (1 - u) * q_sat
    β = (b * a + qc / 2) / max(a + qc / 2, eps(NF))
    condensation = β * M / (1 + L / cₚ * f * dq_sat_dT)
    condensation = clamp(condensation, zero(NF), max(q - u * q_sat, zero(NF)) / Δt_prognostic)

    # clear: condensate evaporates towards the critical relative humidity
    implicit_factor = 1 + L / cₚ * u * dq_sat_dT
    evaporation = max(u * q_sat - q, zero(NF)) / (implicit_factor * closure.evaporation_time_scale * Δt)
    evaporation = min(evaporation, qc / Δt_prognostic)

    cloudy = b > closure.cloud_fraction_threshold
    return ifelse(cloudy, condensation, zero(NF)), ifelse(cloudy, zero(NF), evaporation), b
end

# cloud fraction for autoconversion: Xu-Randall after condensation (relaxation), Sundqvist (closure)
@inline autoconversion_cloud_fraction(::RelaxationClosure, closure_cloud_fraction, cloud_fraction) = cloud_fraction
@inline autoconversion_cloud_fraction(::SundqvistClosure, closure_cloud_fraction, cloud_fraction) = closure_cloud_fraction

# surface pressure of the reference state (0 = no reference yet), and storing the reference state
@inline surface_pressure_reference(::RelaxationClosure, ij, vars) = zero(eltype(vars.parameterizations.surface_pressure))
@inline surface_pressure_reference(::SundqvistClosure, ij, vars) = vars.prognostic.clouds.surface_pressure_reference[ij, 1]

@inline store_reference!(::RelaxationClosure, ij, k, vars, T, q) = nothing
@inline store_surface_pressure_reference!(::RelaxationClosure, ij, vars, pₛ) = nothing

# shift the last call's reference to "two calls before", store this call's as the last
@propagate_inbounds function store_reference!(::SundqvistClosure, ij, k, vars, T, q)
    (; temperature_reference, humidity_reference) = vars.prognostic.clouds
    temperature_reference[ij, k, 1] = temperature_reference[ij, k, 2]
    temperature_reference[ij, k, 2] = T
    humidity_reference[ij, k, 1] = humidity_reference[ij, k, 2]
    humidity_reference[ij, k, 2] = q
    return nothing
end

@propagate_inbounds function store_surface_pressure_reference!(::SundqvistClosure, ij, vars, pₛ)
    (; surface_pressure_reference) = vars.prognostic.clouds
    surface_pressure_reference[ij, 1] = surface_pressure_reference[ij, 2]
    surface_pressure_reference[ij, 2] = pₛ
    return nothing
end

"""$(TYPEDSIGNATURES)
The cloud condensate as a prognostic variable fused with the atmospheric prognostic variables
(spectral state, grid copy and both tendencies), and its flux intermediates `u⋅qc`, `v⋅qc` fused
with the tendencies, so that it is transformed in the batched transforms and advected and
diffused as humidity is, see `ADVECTED_SCALARS`. `nsteps` are the step counts of the time stepper."""
function cloud_condensate_variables(nsteps)
    pg = nsteps.prognostic_grid
    ps = nsteps.prognostic_spectral
    tg = nsteps.tendency_grid
    ts = nsteps.tendency_spectral
    return (
        PrognosticVariable(:cloud_condensate, SpectralXYZT(ps), desc = "Cloud condensate (liquid + ice)", units = "kg/kg", fuse = :prognostic),
        GridVariable(:cloud_condensate, GridXYZT(pg), desc = "Cloud condensate (liquid + ice)", units = "kg/kg", fuse = :grid),
        TendencyVariable(:cloud_condensate, SpectralXYZT(ts), desc = "Tendency of cloud condensate", units = "kg/kg/s", fuse = :spectral_tendencies),
        TendencyVariable(:cloud_condensate, GridXYZT(tg), namespace = :grid, desc = "Tendency of cloud condensate on the grid", units = "kg/kg/s", fuse = :grid_tendencies),
        DynamicsVariable(:uqc, GridXYZT(tg), desc = "u*cloud condensate intermediate on grid", namespace = :grid, fuse = :grid_tendencies),
        DynamicsVariable(:vqc, GridXYZT(tg), desc = "v*cloud condensate intermediate on grid", namespace = :grid, fuse = :grid_tendencies),
        DynamicsVariable(:uqc, SpectralXYZT(ts), desc = "u*cloud condensate intermediate in spectral space", fuse = :spectral_tendencies),
        DynamicsVariable(:vqc, SpectralXYZT(ts), desc = "v*cloud condensate intermediate in spectral space", fuse = :spectral_tendencies),
    )
end

"""$(TYPEDSIGNATURES)
The cloud state that a prognostic cloud scheme writes and radiation reads: per layer the cloud
fraction, grid-mean cloud liquid and ice water and their effective radii, and the liquid and ice
water paths of the column. Declared by the cloud scheme and by the radiation schemes that read it,
so that radiation without a cloud scheme sees zeros (no clouds)."""
cloud_state_variables() = (
    ParameterizationVariable(:cloud_fraction, GridXYZ(), desc = "Cloud fraction", units = "1"),
    ParameterizationVariable(:cloud_liquid_water, GridXYZ(), desc = "Grid-mean cloud liquid water", units = "kg/kg"),
    ParameterizationVariable(:cloud_ice_water, GridXYZ(), desc = "Grid-mean cloud ice water", units = "kg/kg"),
    ParameterizationVariable(:cloud_liquid_effective_radius, GridXYZ(), desc = "Effective radius of cloud liquid", units = "m"),
    ParameterizationVariable(:cloud_ice_effective_radius, GridXYZ(), desc = "Effective radius of cloud ice", units = "m"),
    ParameterizationVariable(:liquid_water_path, Grid2D(), desc = "Liquid water path", units = "kg/m^2"),
    ParameterizationVariable(:ice_water_path, Grid2D(), desc = "Ice water path", units = "kg/m^2"),
)

"""$(TYPEDSIGNATURES)
Fraction of the condensate that is ice at temperature `T` [K]: 0 at and above `temperature_freezing`,
1 at and below `ice_temperature`, linear in between."""
@inline function ice_fraction(T, temperature_freezing, ice_temperature)
    Δ = max(temperature_freezing - ice_temperature, eps(T))
    return clamp((temperature_freezing - T) / Δ, zero(T), one(T))
end

"""$(TYPEDSIGNATURES)
Saturation specific humidity [kg/kg] over a mixture of liquid and ice with ice fraction `f_ice` at
temperature `T` [K] and pressure `p` [Pa], and its derivative with respect to temperature [1/K]
(at fixed `f_ice`). Over liquid this is `saturation_humidity` of the `atmosphere`, over ice
Clausius-Clapeyron with the latent heat of sublimation `Lᵥ + Lᵢ` with the same reference pressure at
freezing. With `mixed_phase = false` saturation is over liquid only."""
@inline function mixed_phase_saturation_humidity(T, p, f_ice, atmosphere, mixed_phase::Bool)
    Lᵥ = atmosphere.latent_heat_condensation
    Lᵢ = atmosphere.latent_heat_fusion
    Rᵥ = atmosphere.R_vapor
    T₀ = atmosphere.temperature_freezing

    q_liquid = saturation_humidity(T, p, atmosphere)
    q_ice = q_liquid * exp(Lᵢ / Rᵥ * (inv(T₀) - inv(T)))     # Lᵢ + Lᵥ in the exponent overall
    f = ifelse(mixed_phase, f_ice, zero(f_ice))
    q_sat = (1 - f) * q_liquid + f * q_ice
    dq_sat_dT = ((1 - f) * q_liquid * Lᵥ + f * q_ice * (Lᵥ + Lᵢ)) / (Rᵥ * T^2)
    return q_sat, dq_sat_dT
end

"""$(TYPEDSIGNATURES)
Cloud fraction from grid-mean condensate `q_condensate` [kg/kg], relative humidity and saturation
humidity `q_sat` [kg/kg] at pressure `p` [Pa], Xu and Randall (1996) in the form and with the
constants of NCEP's GFS (`progcld_zhao_carr`):

    C = RH^¼ [1 - exp(-α q_c / ((1 - RH) q_sat)^¼)]

with the denominator clamped to [1e-4, 1], the exponent to 50, no cloud for condensate below
`1e-6 × p/1000 hPa` and cloud fractions below `cloud_fraction_min` set to zero. No cloud
without condensate."""
@inline function xu_randall_cloud_fraction(q_condensate, relative_humidity, q_sat, p, α, cloud_fraction_min)
    NF = typeof(q_condensate)
    rh = clamp(relative_humidity, NF(1.0e-10), one(NF))         # >0 for a finite derivative of RH^¼
    q_threshold = NF(1.0e-6) * p / NF(1.0e5)
    denominator = clamp(sqrt(sqrt(max(1 - rh, NF(1.0e-10)) * q_sat)), NF(1.0e-4), one(NF))
    exponent = min(α * max(q_condensate, zero(NF)) / denominator, NF(50))
    cloud_fraction = sqrt(sqrt(rh)) * (1 - exp(-exponent))
    is_cloud = (q_condensate > q_threshold) & (cloud_fraction >= cloud_fraction_min)
    return ifelse(is_cloud, cloud_fraction, zero(NF))
end

"""$(TYPEDSIGNATURES)
Cloud liquid [kg/kg] converted to rain over the time step `Δt` [s] (Sundqvist et al. 1989 as in
NCEP's GFS): the liquid above `minimum_condensate` decays with rate
`k = C₀ F [1 - exp(-(q F/(m_r C))²)]`, enhanced by `F` from the collection of precipitation
`precipitation_above` [kg/m²/s] falling into the layer and from the Bergeron-Findeisen process
in supercooled cloud, `C` the cloud fraction. Exponential in time so that it never removes more
than available."""
@inline function liquid_autoconversion(
        q_liquid, cloud_fraction, precipitation_above, T, minimum_condensate, scheme, Δt,
    )
    NF = typeof(q_liquid)
    (; autoconversion_rate, autoconversion_water, collection_enhancement) = scheme
    (; bergeron_findeisen_enhancement, bergeron_findeisen_temperature) = scheme
    excess = max(q_liquid - minimum_condensate, zero(NF))

    # small offsets in the square roots for a finite derivative at zero
    collection = 1 + collection_enhancement * sqrt(max(precipitation_above, zero(NF)) + NF(1.0e-12))
    supercooling = clamp(bergeron_findeisen_temperature - T, zero(NF), NF(20))
    bergeron_findeisen = 1 + bergeron_findeisen_enhancement * sqrt(supercooling + NF(1.0e-6))
    F = collection * bergeron_findeisen

    x = excess * F / (autoconversion_water * max(cloud_fraction, NF(0.01)))
    rate = autoconversion_rate * F * (1 - exp(-min(x^2, NF(50))))
    return excess * (1 - exp(-rate * Δt))
end

"""$(TYPEDSIGNATURES)
Cloud ice [kg/kg] converted to snow over the time step `Δt` [s] (Lin et al. 1983 as in NCEP's
GFS): the ice above `minimum_condensate` decays with rate `c exp(0.025 (T - T₀))`, exponential in time."""
@inline function ice_autoconversion(q_ice, T, temperature_freezing, minimum_condensate, scheme, Δt)
    NF = typeof(q_ice)
    excess = max(q_ice - minimum_condensate, zero(NF))
    rate = scheme.ice_autoconversion_rate * exp(NF(0.025) * (T - temperature_freezing))
    return excess * (1 - exp(-rate * Δt))
end

# function barrier
@propagate_inbounds function parameterization!(ij, vars, scheme::PrognosticCloudCondensation, model)
    (; geometry, planet, atmosphere, land_sea_mask, time_stepping) = model
    cloud_condensation!(ij, vars, scheme, geometry, planet, atmosphere, land_sea_mask, time_stepping)
    return nothing
end

"""$(TYPEDSIGNATURES)
Column `ij` of the prognostic cloud condensation, one top-down sweep over the layers. Per layer:

1. negative condensate (from spectral transport) is filled from vapour,
2. vapour condenses into the condensate or condensate evaporates as the `closure` decides
   (with the latent heat and saturation of the mixed phase),
3. the cloud fraction follows from condensate and relative humidity (Xu-Randall),
4. liquid autoconverts to rain, ice to snow,
5. snow from above melts, rain from above reevaporates,
6. tendencies of temperature, humidity and condensate, the cloud state and precipitation.

All sinks are bounded with the step the state is advanced with (`2Δt` for leapfrog), so that
the condensate stays non-negative. The column conserves water (vapour, condensate,
precipitation) and enthalpy, with the condensate's ice fraction taken at the layer's temperature."""
@propagate_inbounds function cloud_condensation!(
        ij,
        vars,
        scheme::PrognosticCloudCondensation,
        geometry::Geometry,
        planet::AbstractPlanet,
        atmosphere::AbstractAtmosphere,
        land_sea_mask,
        time_stepping,
    )
    # previous time step for an Euler forward step of the parameterizations
    temp = get_prognostic_step(vars.grid.temperature, time_stepping, scheme)
    humid = get_prognostic_step(vars.grid.humidity, time_stepping, scheme)
    condensate = get_prognostic_step(vars.grid.cloud_condensate, time_stepping, scheme)
    temp_tend = get_tendency_step(vars.tendencies.grid.temperature, time_stepping, scheme)
    humid_tend = get_tendency_step(vars.tendencies.grid.humidity, time_stepping, scheme)
    condensate_tend = get_tendency_step(vars.tendencies.grid.cloud_condensate, time_stepping, scheme)
    (; cloud_fraction, cloud_liquid_water, cloud_ice_water) = vars.parameterizations
    (; cloud_liquid_effective_radius, cloud_ice_effective_radius) = vars.parameterizations
    NF = eltype(temp)
    nlayers = size(temp, 2)

    pₛ = vars.parameterizations.surface_pressure[ij]    # surface pressure [Pa]
    coord = geometry.vertical_coordinates
    land_fraction = land_sea_mask.land_fraction[ij]

    # time steps: Δt defines the relaxation time scales as in ImplicitCondensation,
    # the prognostic step (2Δt for leapfrog) is what the tendencies act over and bounds the sinks
    (; Δt) = time_stepping
    Δt_prognostic = default_time_step(time_stepping)

    # thermodynamics
    g = planet.gravity
    ρ = atmosphere.water_density                        # [kg/m³] to convert between [kg/m²/s] and [m/s]
    cₚ = atmosphere.heat_capacity
    Lᵥ = atmosphere.latent_heat_condensation
    Lᵢ = atmosphere.latent_heat_fusion
    T₀ = atmosphere.temperature_freezing

    (; closure, ice_temperature, mixed_phase_saturation) = scheme
    (; cloud_fraction_coefficient, cloud_fraction_min, minimum_condensate) = scheme
    r_ocean = scheme.liquid_effective_radius_ocean
    r_land = scheme.liquid_effective_radius_land

    rain_flux = zero(NF)                # downward rain and snow flux [m/s] (precipitation rate),
    snow_flux = zero(NF)                # zero at the top of the atmosphere
    liquid_water_path = zero(NF)        # [kg/m²]
    ice_water_path = zero(NF)
    cloud_top = vars.parameterizations.cloud_top[ij]
    pₛ_reference = surface_pressure_reference(closure, ij, vars)  # Sundqvist closure only

    for k in 1:nlayers
        Tₖ = temp[ij, k]
        qₖ = humid[ij, k]
        pₖ = pressure(k, pₛ, coord)
        Δpₖ = pressure_thickness(k, pₛ, coord)
        Δp_gρ = Δpₖ / (g * ρ)           # humidity rate [kg/kg/s] × Δp_gρ = precipitation rate [m/s]

        f_ice = ice_fraction(Tₖ, T₀, ice_temperature)
        L = Lᵥ + f_ice * Lᵢ             # latent heat between vapour and condensate
        q_sat, dq_sat_dT = mixed_phase_saturation_humidity(Tₖ, pₖ, f_ice, atmosphere, mixed_phase_saturation)

        # 1. FILL NEGATIVE CONDENSATE FROM VAPOUR (as far as vapour is available), NCEP GFS
        fill = min(max(-condensate[ij, k], zero(NF)), max(qₖ, zero(NF)))   # [kg/kg]
        qc₁ = max(condensate[ij, k] + fill, zero(NF))
        q₁ = qₖ - fill

        # 2. CONDENSATION of vapour into the condensate and EVAPORATION of condensate [kg/kg/s],
        # as the closure decides
        layer = (;
            T = Tₖ, q_lagged = qₖ, q = q₁, qc = qc₁, p = pₖ, pₛ, pₛ_reference, q_sat, dq_sat_dT, L, cₚ,
            Δt, Δt_prognostic, coord,
        )
        condensation, evaporation, closure_cloud_fraction = condensation_closure(closure, scheme, ij, k, vars, layer)
        net_condensation = condensation - evaporation                         # vapour to condensate [kg/kg/s]
        qc₂ = qc₁ + Δt_prognostic * net_condensation
        q₂ = q₁ - Δt_prognostic * net_condensation

        # 3. CLOUD FRACTION after condensation, for radiation and (relaxation closure) autoconversion
        cloud_fraction₂ = xu_randall_cloud_fraction(qc₂, q₂ / q_sat, q_sat, pₖ, cloud_fraction_coefficient, cloud_fraction_min)
        cloud_fraction_autoconversion = autoconversion_cloud_fraction(closure, closure_cloud_fraction, cloud_fraction₂)

        # 4. AUTOCONVERSION of liquid to rain, ice to snow [kg/kg over Δt_prognostic]
        condensate_min = minimum_condensate * pₖ / 100000     # scaled with pressure as in GFS
        precipitation_above = ρ * (rain_flux + snow_flux)     # [kg/m²/s]
        to_rain = liquid_autoconversion(
            (1 - f_ice) * qc₂, cloud_fraction_autoconversion, precipitation_above, Tₖ, condensate_min, scheme, Δt_prognostic
        )
        to_snow = ice_autoconversion(f_ice * qc₂, Tₖ, T₀, condensate_min, scheme, Δt_prognostic)
        qc₃ = qc₂ - to_rain - to_snow

        # 5. PRECIPITATION FROM ABOVE: melting of snow (energy-limited) and reevaporation of rain
        # (not beyond saturation), as in ImplicitCondensation
        melt_max = cₚ / Lᵢ * max(Tₖ - scheme.melting_threshold, zero(NF)) * Δp_gρ / Δt_prognostic  # [m/s]
        melt = min(snow_flux, melt_max)
        snow_flux -= melt
        rain_flux += melt
        subsaturation = max(q_sat - qₖ, zero(NF))
        rain_evaporated = min(
            min(scheme.reevaporation * subsaturation, one(NF)) * rain_flux,
            subsaturation * Δp_gρ / Δt_prognostic
        )
        rain_flux -= rain_evaporated

        # 6. this layer's autoconversion falls out (into the layer below)
        rain_flux += to_rain / Δt_prognostic * Δp_gρ
        snow_flux += to_snow / Δt_prognostic * Δp_gρ

        # TENDENCIES [kg/kg/s] and [K/s]
        vapour_to_condensate = fill / Δt_prognostic + net_condensation

        # state after this scheme's condensation, the reference for the next supply (Sundqvist closure)
        store_reference!(closure, ij, k, vars, Tₖ + Δt_prognostic * L / cₚ * vapour_to_condensate, qₖ - Δt_prognostic * vapour_to_condensate)
        rain_rate = to_rain / Δt_prognostic
        snow_rate = to_snow / Δt_prognostic
        melt_rate = melt / Δp_gρ
        evaporation_rate = rain_evaporated / Δp_gρ

        condensate_tend[ij, k] += vapour_to_condensate - rain_rate - snow_rate
        humid_tend[ij, k] += evaporation_rate - vapour_to_condensate

        # latent heat: vapour to condensate of ice fraction f_ice, the ice of the condensate
        # converted to snow is replaced by freezing (or the liquid to rain by melting) to keep
        # the condensate at f_ice, melting of snow, reevaporation of rain
        temp_tend[ij, k] += (
            L * vapour_to_condensate +
                Lᵢ * ((1 - f_ice) * snow_rate - f_ice * rain_rate) -
                Lᵢ * melt_rate -
                Lᵥ * evaporation_rate
        ) / cₚ

        # CLOUD STATE for radiation and output
        qc_new = max(qc₃, zero(NF))
        cloud_fraction[ij, k] = cloud_fraction₂
        cloud_liquid_water[ij, k] = (1 - f_ice) * qc_new
        cloud_ice_water[ij, k] = f_ice * qc_new
        r_liquid_land = r_land + f_ice * (r_ocean - r_land)     # GFS: 5 + 5 f_ice μm over land
        cloud_liquid_effective_radius[ij, k] = (1 - land_fraction) * r_ocean + land_fraction * r_liquid_land
        cloud_ice_effective_radius[ij, k] = scheme.ice_effective_radius
        liquid_water_path += (1 - f_ice) * qc_new * Δpₖ / g
        ice_water_path += f_ice * qc_new * Δpₖ / g

        # large-scale condensation raises the cloud top (if above the convective cloud top) as in
        # ImplicitCondensation, read by DiagnosticClouds
        cloud_top = min(cloud_top, ifelse(condensation > 0, NF(k), NF(nlayers + 1)))
    end

    vars.parameterizations.cloud_top[ij] = cloud_top
    store_surface_pressure_reference!(closure, ij, vars, pₛ)
    vars.parameterizations.liquid_water_path[ij] = liquid_water_path
    vars.parameterizations.ice_water_path[ij] = ice_water_path

    # avoid negative precipitation from rounding errors
    rain_flux = max(rain_flux, zero(NF))
    snow_flux = max(snow_flux, zero(NF))

    # accumulated precipitation [m] over the time step the clock advances, and rates [m/s]
    vars.parameterizations.rain_large_scale[ij] += Δt * rain_flux
    vars.parameterizations.snow_large_scale[ij] += Δt * snow_flux
    vars.parameterizations.rain_rate_large_scale[ij] = rain_flux
    vars.parameterizations.snow_rate_large_scale[ij] = snow_flux

    # accumulate into total rain/snow rate including convection [m/s]
    vars.parameterizations.rain_rate[ij] += rain_flux
    vars.parameterizations.snow_rate[ij] += snow_flux
    return nothing
end
