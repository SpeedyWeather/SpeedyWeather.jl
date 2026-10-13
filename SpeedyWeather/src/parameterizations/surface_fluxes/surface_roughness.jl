abstract type AbstractSurfaceRoughness <: AbstractParameterization end

export ConstantSurfaceRoughness
export LearnedSurfaceRoughness

variables(::AbstractSurfaceRoughness) = (
    # Momentum, heat and moisture roughness lengths
    ParameterizationVariable(:momentum_roughness, Grid2D(), desc = "Momentum roughness length", units = "m"),
    ParameterizationVariable(:heat_roughness, Grid2D(), desc = "Heat roughness length", units = "m"),
    ParameterizationVariable(:humidity_roughness, Grid2D(), desc = "Moisture roughness length", units = "m"),

    # Momentum, heat and moisture roughness lengths for the land
    ParameterizationVariable(:momentum_roughness, Grid2D(), desc = "Land momentum roughness length", units = "m", namespace = :land),
    ParameterizationVariable(:heat_roughness, Grid2D(), desc = "Land heat roughness length", units = "m", namespace = :land),
    ParameterizationVariable(:humidity_roughness, Grid2D(), desc = "Land moisture roughness length", units = "m", namespace = :land),

    # Momentum, heat and moisture roughness lengths for the ocean
    ParameterizationVariable(:momentum_roughness, Grid2D(), desc = "Ocean momentum roughness length", units = "m", namespace = :ocean),
    ParameterizationVariable(:heat_roughness, Grid2D(), desc = "Ocean heat roughness length", units = "m", namespace = :ocean),
    ParameterizationVariable(:humidity_roughness, Grid2D(), desc = "Ocean moisture roughness length", units = "m", namespace = :ocean),
)

"""
Surface roughness length parameterization with constant roughness length
over land and ocean. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct ConstantSurfaceRoughness{NF} <: AbstractSurfaceRoughness
    "[OPTION] constant momentum roughness length over land [m]"
    @param momentum_roughness_land::NF = 0.5 (bounds = Nonnegative,)

    "[OPTION] constant heat roughness length over land [m]"
    @param heat_roughness_land::NF = 0.5 (bounds = Nonnegative,)

    "[OPTION] constant humidity roughness length over land [m]"
    @param humidity_roughness_land::NF = 0.5 (bounds = Nonnegative,)

    "[OPTION] constant momentum roughness length over ocean [m]"
    @param momentum_roughness_ocean::NF = 1.0e-4 (bounds = Nonnegative,)

    "[OPTION] constant heat roughness length over ocean [m]"
    @param heat_roughness_ocean::NF = 1.0e-4 (bounds = Nonnegative,)

    "[OPTION] constant humidity roughness length over ocean [m]"
    @param humidity_roughness_ocean::NF = 1.55e-4 (bounds = Nonnegative,)

end

Adapt.@adapt_structure ConstantSurfaceRoughness
ConstantSurfaceRoughness(SG::SpectralGrid; kwargs...) = ConstantSurfaceRoughness{SG.NF}(; kwargs...)
initialize!(::ConstantSurfaceRoughness, ::PrimitiveEquation) = nothing

@propagate_inbounds function parameterization!(ij, vars, scheme::ConstantSurfaceRoughness, model)
    surface_roughness!(ij, vars, scheme, model.land_sea_mask)
    return nothing
end

@inline function constant_roughness(z₀_land, z₀_ocean, land_fraction)
    land = ifelse(land_fraction > 0, z₀_land, zero(z₀_land))
    ocean = ifelse(land_fraction < 1, z₀_ocean, zero(z₀_ocean))
    combined = land_fraction * land + (1 - land_fraction) * ocean
    return land, ocean, combined
end

@propagate_inbounds function surface_roughness!(ij, vars, scheme::ConstantSurfaceRoughness, land_sea_mask)
    land_fraction = land_sea_mask.land_fraction[ij]
    z₀M_land = scheme.momentum_roughness_land
    z₀M_ocean = scheme.momentum_roughness_ocean
    z₀H_land = scheme.heat_roughness_land
    z₀H_ocean = scheme.heat_roughness_ocean
    z₀Q_land = scheme.humidity_roughness_land
    z₀Q_ocean = scheme.humidity_roughness_ocean

    (; land, ocean) = vars.parameterizations

    # Momentum roughness
    land.momentum_roughness[ij], ocean.momentum_roughness[ij], vars.parameterizations.momentum_roughness[ij] = constant_roughness(z₀M_land, z₀M_ocean, land_fraction)

    # Heat roughness
    land.heat_roughness[ij], ocean.heat_roughness[ij], vars.parameterizations.heat_roughness[ij] = constant_roughness(z₀H_land, z₀H_ocean, land_fraction)

    # Humidity roughness (just heat times a constant)
    land.humidity_roughness[ij], ocean.humidity_roughness[ij], vars.parameterizations.humidity_roughness[ij] = constant_roughness(z₀Q_land, z₀Q_ocean, land_fraction)

    return nothing
end

export LearnedLandRoughness

"""
Learned surface roughness lengths over land for momentum, heat and humidity
fluxes, using symbolic regression expressions.

The momentum roughness length is

    z₀M = exp(max(m_c1, A + B) / (C + D) + m_c9)

where

    A = m_c4 * exp(S / (m_c2*C_H - m_c3))
    B = exp(C_L * (min(m_c5*L_L, m_c6*T - W) - m_c7^L_L))
    C = exp(C_H * max(L_H, m_c8))
    D = min(L_H, W)

The heat roughness length is

    z₀H = exp(min(A, B))

where

    Γ = max(
        (h_c1^h_c2 * min(C_H, L_H)) / (L_L + T)^h_c2 + L_H,
        2*C_H + h_c3
    )

    A = h_c5 * (Γ - h_c4 + C_L*min(C_L, L_L - W)) - h_c6*L_H
    B = h_c8 * h_c7^C_H

The humidity roughness length is taken to equal the heat roughness length,

    z₀Q = z₀H

Predictors:
- C_H, C_L: High and low vegetation cover fractions
- L_H, L_L: High and low vegetation leaf area indices
- S: Snow depth
- T: Soil temperature at the bottom soil layer
- W: Soil moisture at the top soil layer

Fields are $(TYPEDFIELDS)
"""
@parameterized @kwdef struct LearnedLandRoughness{NF} <: AbstractSurfaceRoughness
    "[OPTION] maximum land surface momentum roughness length [m]"
    momentum_roughness_max::NF = 2.0

    "[OPTION] maximum land surface heat roughness length [m]"
    heat_roughness_max::NF = 2.0

    "[OPTION] minimum land surface momentum roughness length [m]"
    momentum_roughness_min::NF = 1.3e-3

    "[OPTION] minimum land surface heat roughness length [m]"
    heat_roughness_min::NF = 1.3e-4

    "[OPTION] coefficient for momentum roughness closure over land"
    m_c1::NF = -7.231684
    m_c2::NF = 11.078163
    m_c3::NF = -0.14265412
    m_c4::NF = -5.916258
    m_c5::NF = 2.2420337
    m_c6::NF = 0.00613743
    m_c7::NF = 0.51263654
    m_c8::NF = 2.6413066
    m_c9::NF = 0.60989547

    "[OPTION] coefficient for heat roughness closure over land"
    h_c1::NF = 373.86444
    h_c2::NF = 5.981257
    h_c3::NF = 3.5229454
    h_c4::NF = 7.781359
    h_c5::NF = 2.1035099
    h_c6::NF = -0.41931123
    h_c7::NF = 0.00037635124
    h_c8::NF = -5.4191554
end

Adapt.@adapt_structure LearnedLandRoughness
LearnedLandRoughness(SG::SpectralGrid; kwargs...) = LearnedLandRoughness{SG.NF}(; kwargs...)
initialize!(::LearnedLandRoughness, ::PrimitiveEquation) = nothing

export Charnock

"""
Wind-dependent Charnock parameterization using a piecewise polynomial
function of neutral wind speed.

The Charnock parameter α is calculated as

    α = c1*uₙ^3 + c2*uₙ^2 + c3       for uₙ ≤ Uₜ
    α = c4*(uₙ - Uₜ) + c5           for uₙ > Uₜ

where
- uₙ: Neutral wind speed [m/s]
- Uₜ: Threshold wind speed [m/s]
- c1–c5: Fitted coefficients

The parameterization returns log(α), which is used as an input to
the learned ocean momentum roughness parameterization.

Fields are $(TYPEDFIELDS)
"""
@parameterized @kwdef struct Charnock{NF} <: AbstractParameterization
    "[OPTION] threshold wind speed [m/s] for Charnock parameterization"
    Uₜ::NF = 7.0

    "[OPTION] low wind cubic term coefficient"
    c1::NF = 1.8449e-6

    "[OPTION] low wind quadratic term coefficient"
    c2::NF = 9.7104e-5

    "[OPTION] low wind constant term coefficient"
    c3::NF = 0.0075

    "[OPTION] high wind linear term coefficient"
    c4::NF = 0.0016

    "[OPTION] high wind constant term coefficient"
    c5::NF = 0.0129

    "[OPTION] minimum Charnock parameter value"
    charnock_min::NF = 1.0e-5
end

Adapt.@adapt_structure Charnock
Charnock(SG::SpectralGrid; kwargs...) = Charnock{SG.NF}(; kwargs...)

@inline function (charnock_param::Charnock{NF})(uₙ::T) where {NF, T <: Real}
    Uₜ = T(charnock_param.Uₜ)
    c1 = T(charnock_param.c1)
    c2 = T(charnock_param.c2)
    c3 = T(charnock_param.c3)
    c4 = T(charnock_param.c4)
    c5 = T(charnock_param.c5)
    c_min = T(charnock_param.charnock_min)

    α = ifelse(
        uₙ ≤ Uₜ,
        muladd(uₙ * uₙ, muladd(c1, uₙ, c2), c3),
        muladd(c4, uₙ - Uₜ, c5) # for uₙ > Uₜ, linear increase
    )
    return max(α, c_min)  # clamp to lower bound of Charnock parameter
end

export LearnedOceanRoughness

"""
Learned surface roughness lengths over ocean for momentum, heat and humidity
fluxes, using symbolic regression expressions and sea-ice corrections.

The open-ocean momentum roughness length is

    z₀M_base = exp(A - B + C)

where

    A = u^m_c6
    B = (m_c1*u + α + m_c2) / (exp(m_c3*u - m_c4) + m_c5)
    C = m_c7*α - m_c8

Here α = log(α_charnock), obtained from the Charnock parameterization.

The open-ocean momentum roughness is capped,

    z₀M_open = min(z₀M_base, z₀M_max)

where z₀M_max is determined by the wind-dependent
maximum_momentum_roughness function.

The final momentum roughness includes sea-ice blending,

    z₀M = s*z₀M_ice + (1 - s)*z₀M_open

where

    z₀M_ice = max(
        0.93e-3*(1 - s) + 6.05e-3*exp(-17*(s - 0.5)^2),
        1.0e-3
    )

The heat roughness length is

    z₀H = s*1.0e-3 + (1 - s)*exp(A - B)

where

    A = h_c1*(u + h_c2)^h_c3
    B = exp(h_c4/u)

The humidity roughness length is proportional to the heat
roughness length,

    z₀Q = z₀H * (a_q / a_h)

Predictors:
- u: Surface wind speed [m/s]
- s: Sea-ice concentration
- α: Logarithm of the Charnock parameter, evaluated using neutral wind
- a_h, a_q: Heat and humidity scaling coefficients

Fields are $(TYPEDFIELDS)
"""
@parameterized @kwdef struct LearnedOceanRoughness{NF, CH} <: AbstractSurfaceRoughness
    # Constants from IFS documentation
    a_h::NF = 0.4 # heat related constant
    a_q::NF = 0.62 # moisture related constant

    a_ratio::NF = a_q / a_h # moisture to heat roughness ratio

    "[OPTION] Charnock parameter estimation"
    charnock_param::CH

    "[OPTION] coefficient for momentum roughness closure over ocean"
    m_c1::NF = 0.48786303
    m_c2::NF = 3.374574
    m_c3::NF = 0.5460857
    m_c4::NF = 0.45002246
    m_c5::NF = 0.9827853
    m_c6::NF = 0.043734703
    m_c7::NF = 1.0346544
    m_c8::NF = 7.1816134

    "[OPTION] minimum ocean momentum roughness length [m]"
    momentum_roughness_min::NF = 1.0e-7

    "[OPTION] coefficient for heat roughness closure over ocean"
    h_c1::NF = -8.207549
    h_c2::NF = 0.18300448
    h_c3::NF = 0.08970642
    h_c4::NF = -1.6667988

    "[OPTION] sea-ice heat roughness length"  # defined in IFS documentation
    sea_ice_heat_roughness::NF = 1.0e-3
end

Adapt.@adapt_structure LearnedOceanRoughness
function LearnedOceanRoughness(SG::SpectralGrid; charnock_param = Charnock(SG), kwargs...)
    return LearnedOceanRoughness{SG.NF, typeof(charnock_param)}(; charnock_param, kwargs...)
end
initialize!(::LearnedOceanRoughness, ::PrimitiveEquation) = nothing

export LearnedOceanRoughness

"""
Combined learned surface roughness parameterization for momentum, heat
and humidity fluxes over land and ocean.

Land and ocean roughness lengths are calculated independently using
LearnedLandRoughness and LearnedOceanRoughness, respectively.

The grid-cell mean roughness length for each quantity X ∈ {M, H, Q}
is obtained through area-weighted blending,

    z₀X = f_land*z₀X_land + (1 - f_land)*z₀X_ocean

where
- f_land: Land fraction of the grid cell
- z₀X_land: Learned land roughness length
- z₀X_ocean: Learned ocean roughness length

Fields are $(TYPEDFIELDS)
"""
@parameterized @kwdef struct LearnedSurfaceRoughness{LLR, LOR} <: AbstractSurfaceRoughness
    "[OPTION] learned land surface roughness scheme"
    @component land_surface_roughness::LLR

    "[OPTION] learned ocean surface roughness scheme"
    @component ocean_surface_roughness::LOR
end

Adapt.@adapt_structure LearnedSurfaceRoughness
function LearnedSurfaceRoughness(
        SG::SpectralGrid;
        land_surface_roughness = LearnedLandRoughness(SG),
        ocean_surface_roughness = LearnedOceanRoughness(SG),
    )
    return LearnedSurfaceRoughness(land_surface_roughness, ocean_surface_roughness)
end

function initialize!(LSR::LearnedSurfaceRoughness, model::PrimitiveEquation)
    initialize!(LSR.land_surface_roughness, model)
    initialize!(LSR.ocean_surface_roughness, model)
    return nothing
end

# First land, then ocean, then combine!
@propagate_inbounds function parameterization!(ij, vars, LSR::LearnedSurfaceRoughness, model)
    surface_roughness!(ij, vars, LSR.land_surface_roughness, model.land_sea_mask)
    surface_roughness!(ij, vars, LSR.ocean_surface_roughness, model)
    surface_roughness!(ij, vars, LSR, model.land_sea_mask)
    return nothing
end

@propagate_inbounds function surface_roughness!(ij, vars, LSR::LearnedSurfaceRoughness, land_sea_mask)
    land_fraction = land_sea_mask.land_fraction[ij]
    (; land, ocean) = vars.parameterizations

    # Momentum roughness
    vars.parameterizations.momentum_roughness[ij] = land_fraction * land.momentum_roughness[ij] +
        (1 - land_fraction) * ocean.momentum_roughness[ij]

    # Heat roughness
    vars.parameterizations.heat_roughness[ij] = land_fraction * land.heat_roughness[ij] +
        (1 - land_fraction) * ocean.heat_roughness[ij]

    # Moisture roughness
    vars.parameterizations.humidity_roughness[ij] = land_fraction * land.humidity_roughness[ij] +
        (1 - land_fraction) * ocean.humidity_roughness[ij]
    return nothing
end

@propagate_inbounds function surface_roughness!(ij, vars, scheme::LearnedLandRoughness, land_sea_mask)
    land_fraction = land_sea_mask.land_fraction[ij]
    (; land) = vars.parameterizations

    if land_fraction <= 0
        zero_val = zero(land_fraction)
        land.momentum_roughness[ij] = zero_val
        land.heat_roughness[ij] = zero_val
        land.humidity_roughness[ij] = zero_val
        return nothing
    end


    # Unpack predictors, using defaults when optional fields are absent
    (; land) = vars.parameterizations

    vₕ = haskey(land, :high_vegetation_cover) ? land.high_vegetation_cover[ij] : zero(land_fraction)
    vₗ = haskey(land, :low_vegetation_cover) ? land.low_vegetation_cover[ij] : zero(land_fraction)
    laiₕ = haskey(land, :high_vegetation_leaf_area) ? land.high_vegetation_leaf_area[ij] : zero(land_fraction)
    laiₗ = haskey(land, :low_vegetation_leaf_area) ? land.low_vegetation_leaf_area[ij] : zero(land_fraction)
    sd = haskey(vars.prognostic.land, :snow_depth) ? vars.prognostic.land.snow_depth[ij] : zero(land_fraction)
    soil_moisture = haskey(vars.prognostic.land, :soil_moisture) ? vars.prognostic.land.soil_moisture[ij, begin] : zero(land_fraction)
    soil_temperature = haskey(vars.prognostic.land, :soil_temperature) ? vars.prognostic.land.soil_temperature[ij, end] : zero(land_fraction)

    z₀M_land = land_momentum_roughness(scheme, vₕ, vₗ, laiₕ, laiₗ, sd, soil_temperature, soil_moisture)
    z₀H_land = land_heat_roughness(scheme, vₕ, vₗ, laiₕ, laiₗ, soil_temperature, soil_moisture)

    land.momentum_roughness[ij] = z₀M_land
    land.heat_roughness[ij] = z₀H_land
    land.humidity_roughness[ij] = z₀H_land # same as heat roughness over land

    return nothing
end


@propagate_inbounds function surface_roughness!(ij, vars, scheme::LearnedOceanRoughness, model)
    (; ocean) = vars.parameterizations
    (; land_sea_mask) = model
    land_fraction = land_sea_mask.land_fraction[ij]

    if land_fraction >= 1
        zero_val = zero(land_fraction)
        ocean.momentum_roughness[ij] = zero_val
        ocean.heat_roughness[ij] = zero_val
        ocean.humidity_roughness[ij] = zero_val
        return nothing
    end

    (; surface_wind_speed, neutral_wind_speed) = vars.parameterizations
    u = surface_wind_speed[ij]
    progn_ocean = vars.prognostic.ocean

    ℵ = if haskey(progn_ocean, :sea_ice_concentration)
        ice = get_prognostic_step(
            progn_ocean.sea_ice_concentration,
            time_stepper(model.time_stepping, :ocean),
            scheme,
        )
        ice[ij]
    else
        zero(u)
    end
    log_αᵪ = log(scheme.charnock_param(neutral_wind_speed[ij]))

    z₀M_ocean = ocean_momentum_roughness(scheme, surface_wind_speed[ij], log_αᵪ, ℵ)
    z₀H_ocean = ocean_heat_roughness(scheme, surface_wind_speed[ij], ℵ)

    ocean.momentum_roughness[ij] = z₀M_ocean
    ocean.heat_roughness[ij] = z₀H_ocean
    ocean.humidity_roughness[ij] = z₀H_ocean * scheme.a_ratio # moisture roughness is heat times a constant

    return nothing
end

@inline function sea_ice_momentum_roughness(ℵ::NF) where {NF <: Real}
    # Following IFS documentation
    base = NF(0.93e-3) * (1 - ℵ) + NF(6.05e-3) * exp(NF(-17) * (ℵ - NF(0.5))^2)
    return max(base, NF(1.0e-3))
end

@inline function maximum_momentum_roughness(uₙ::NF) where {NF <: Real}
    # Momentum roughness cap following Willett et al. (2026)
    z_cap = NF(6.74e-3)
    z_high = NF(1.3e-3)
    U_cap = NF(33.0)
    U_high = NF(55.0)

    m = (z_high - z_cap) / (U_high - U_cap)
    z_linear = muladd(m, uₙ - U_cap, z_cap)

    return clamp(z_linear, z_high, z_cap)
end

# Ocean momentum roughness
@inline function ocean_momentum_roughness(scheme, uₙ::NF, log_αᵪ::NF, ℵ::NF) where {NF <: AbstractFloat}

    num = muladd(scheme.m_c3, uₙ, log_αᵪ + exp(scheme.m_c2))
    den = exp(muladd(scheme.m_c4, uₙ, - scheme.m_c5)) + scheme.m_c6

    z₀ = exp(
        uₙ^scheme.m_c1 - num / den +
            muladd(scheme.m_c7, log_αᵪ, -scheme.m_c8)
    )

    z_max = maximum_momentum_roughness(u)
    z₀ = clamp(z₀, NF(scheme.momentum_roughness_min), z_max)

    z_ice = sea_ice_momentum_roughness(ℵ)

    return muladd(ℵ, z_ice - z₀, z₀)
end


# Ocean heat roughness
@inline function ocean_heat_roughness(scheme, uₙ::NF, ℵ::NF) where {NF <: AbstractFloat}
    z₀ = exp(scheme.h_c1 * (uₙ + scheme.h_c2)^scheme.h_c3 - exp(scheme.h_c4 / uₙ))

    return muladd(ℵ, scheme.sea_ice_heat_roughness - z₀, z₀) # weight by sea-ice concentration
end


# Land momentum roughness
@inline function land_momentum_roughness(
        scheme, high_vegetation_cover::NF, low_vegetation_cover::NF,
        high_vegetation_leaf_area::NF, low_vegetation_leaf_area::NF,
        snow_depth::NF, soil_temperature_bottom::NF, soil_moisture_bottom::NF
    ) where {NF <: AbstractFloat}
    # From Munday et al. 2026, in prep.
    C_H = high_vegetation_cover
    C_L = low_vegetation_cover
    L_H = high_vegetation_leaf_area
    L_L = low_vegetation_leaf_area
    S = snow_depth
    T = soil_temperature_bottom
    W = soil_moisture_bottom

    # Numerator
    den1 = muladd(scheme.m_c2, C_H, -scheme.m_c3)
    term1 = scheme.m_c4 * exp(S / den1)

    min_inner = min(
        scheme.m_c5 * L_L,
        muladd(scheme.m_c6, T, -W)
    )

    term2 = exp(
        C_L * (min_inner - scheme.m_c7^L_L)
    )

    numerator = max(scheme.m_c1, term1 + term2)

    # Denominator
    term3 = exp(C_H * max(L_H, scheme.m_c8))
    term4 = min(L_H, W)
    denominator = term3 + term4

    z₀ = exp(numerator / denominator + scheme.m_c9)

    return clamp(z₀, scheme.momentum_roughness_min, scheme.momentum_roughness_max)  # Physical limits from IFS documentation
end


# Land heat roughness
@inline function land_heat_roughness(
        scheme, high_vegetation_cover::NF, low_vegetation_cover::NF,
        high_vegetation_leaf_area::NF, low_vegetation_leaf_area::NF,
        soil_temperature_bottom::NF, soil_moisture_bottom::NF
    ) where {NF <: AbstractFloat}
    # From Munday et al. 2026, in prep.
    C_H = high_vegetation_cover
    C_L = low_vegetation_cover
    L_H = high_vegetation_leaf_area
    L_L = low_vegetation_leaf_area
    T = soil_temperature_bottom
    W = soil_moisture_bottom

    # Gamma
    gamma_term1_num = (scheme.h_c1^scheme.h_c2) * min(C_H, L_H)
    gamma_term1_den = (L_L + T)^scheme.h_c2
    gamma_term1 = gamma_term1_num / gamma_term1_den + L_H

    gamma_term2 = muladd(NF(2), C_H, scheme.h_c3)
    Γ = max(gamma_term1, gamma_term2)

    # Learned heat roughness
    inner_min = min(C_L, L_L - W)

    arg1 = muladd(
        scheme.h_c5,
        Γ - scheme.h_c4 + C_L * inner_min,
        - (scheme.h_c6 * L_H)
    )

    arg2 = scheme.h_c8 * (scheme.h_c7^C_H)

    z₀ = exp(min(arg1, arg2))

    return clamp(z₀, scheme.heat_roughness_min, scheme.heat_roughness_max)  # Physical limits from IFS documentation
end
