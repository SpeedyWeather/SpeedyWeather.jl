abstract type AbstractSurfaceRoughness <: AbstractParameterization end

export ConstantSurfaceRoughness
export LearnedSurfaceRoughness

variables(::AbstractSurfaceRoughness) = (
    # Momentum, heat and moisture roughness lengths
    ParameterizationVariable(:momentum_roughness, Grid2D(), desc = "Momentum roughness length", units = "m"),
    ParameterizationVariable(:heat_roughness, Grid2D(), desc = "Heat roughness length", units = "m"),
    ParameterizationVariable(:moisture_roughness, Grid2D(), desc = "Moisture roughness length", units = "m"),

    # Momentum, heat and moisture roughness lengths for the land
    ParameterizationVariable(:momentum_roughness, Grid2D(), desc = "Land momentum roughness length", units = "m", namespace = :land),
    ParameterizationVariable(:heat_roughness, Grid2D(), desc = "Land heat roughness length", units = "m", namespace = :land),
    ParameterizationVariable(:moisture_roughness, Grid2D(), desc = "Land moisture roughness length", units = "m", namespace = :land),

    # Momentum, heat and moisture roughness lengths for the ocean
    ParameterizationVariable(:momentum_roughness, Grid2D(), desc = "Ocean momentum roughness length", units = "m", namespace = :ocean),
    ParameterizationVariable(:heat_roughness, Grid2D(), desc = "Ocean heat roughness length", units = "m", namespace = :ocean),
    ParameterizationVariable(:moisture_roughness, Grid2D(), desc = "Ocean moisture roughness length", units = "m", namespace = :ocean),
)

"""
Surface roughness length parameterization with constant roughness length
over land and ocean. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct ConstantSurfaceRoughness{NF} <: AbstractSurfaceRoughness
    "[OPTION] constant momentum roughness length over land [m]"
    @param momentum_roughness_land::NF = 0.5 (bounds = Nonnegative,)

    "[OPTION] constant heat roughness length over land [m]"
    @param heat_roughness_land::NF = 0.5 (bounds = Nonnegative,)

    "[OPTION] constant momentum roughness length over ocean [m]"
    @param momentum_roughness_ocean::NF = 1.0e-4 (bounds = Nonnegative,)

    "[OPTION] constant heat roughness length over ocean [m]"
    @param heat_roughness_ocean::NF = 1.0e-4 (bounds = Nonnegative,)

end

Adapt.@adapt_structure ConstantSurfaceRoughness
ConstantSurfaceRoughness(SG::SpectralGrid; kwargs...) = ConstantSurfaceRoughness{SG.NF}(; kwargs...)
initialize!(::ConstantSurfaceRoughness, ::PrimitiveEquation) = nothing


@propagate_inbounds function surface_roughness!(ij, vars, scheme::ConstantSurfaceRoughness, land_sea_mask)
    land_fraction = land_sea_mask.land_fraction[ij]
    z₀M_land = scheme.momentum_roughness_land
    z₀M_ocean = scheme.momentum_roughness_ocean
    z₀H_land = scheme.heat_roughness_land
    z₀H_ocean = scheme.heat_roughness_ocean

    (; land, ocean) = vars.parameterizations

    # Momentum roughness
    land.momentum_roughness[ij] = ifelse(land_fraction > 0, z₀M_land, zero(z₀M_land))
    ocean.momentum_roughness[ij] = ifelse(land_fraction < 1, z₀M_ocean, zero(z₀M_ocean))
    vars.parameterizations.momentum_roughness[ij] = land_fraction * land.momentum_roughness[ij] +
        (1 - land_fraction) * ocean.momentum_roughness[ij]

    # Heat roughness
    land.heat_roughness[ij] = ifelse(land_fraction > 0, z₀H_land, zero(z₀H_land))
    ocean.heat_roughness[ij] = ifelse(land_fraction < 1, z₀H_ocean, zero(z₀H_ocean))
    vars.parameterizations.heat_roughness[ij] = land_fraction * land.heat_roughness[ij] +
        (1 - land_fraction) * ocean.heat_roughness[ij]

    # Moisture roughness (just heat times a constant)
    land.moisture_roughness[ij] = ifelse(land_fraction > 0, z₀H_land * (0.62 / 0.4), zero(z₀H_land))
    ocean.moisture_roughness[ij] = ifelse(land_fraction < 1, z₀H_ocean * (0.62 / 0.4), zero(z₀H_ocean))
    vars.parameterizations.moisture_roughness[ij] = land_fraction * land.moisture_roughness[ij] +
        (1 - land_fraction) * ocean.moisture_roughness[ij]

    return nothing
end

@kwdef struct LearnedLandRoughness{NF} <: AbstractSurfaceRoughness
    "[OPTION] constant momentum roughness length over land [m]"
    momentum_roughness_land::NF = 0.5

    "[OPTION] constant heat roughness length over land [m]"
    heat_roughness_land::NF = 0.5
end

Adapt.@adapt_structure LearnedLandRoughness
LearnedLandRoughness(SG::SpectralGrid; kwargs...) = LearnedLandRoughness{SG.NF}(; kwargs...)
initialize!(::LearnedLandRoughness, ::PrimitiveEquation) = nothing

@kwdef struct LearnedOceanRoughness{NF} <: AbstractSurfaceRoughness
    a_h::NF = NF(0.4) # heat related constant
    a_q::NF = NF(0.62) # moisture related constant
end

Adapt.@adapt_structure LearnedOceanRoughness
LearnedOceanRoughness(SG::SpectralGrid; kwargs...) = LearnedOceanRoughness{SG.NF}(; kwargs...)
initialize!(::LearnedOceanRoughness, ::PrimitiveEquation) = nothing

"""
Learned ocean surface roughness for momentum, heat and moisture fluxes."""
@kwdef struct LearnedSurfaceRoughness{LLR, LOR} <: AbstractSurfaceRoughness
    "[OPTION] learned land surface roughness scheme"
    land_surface_roughness::LLR

    "[OPTION] learned ocean surface roughness scheme"
    ocean_surface_roughness::LOR

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
    surface_roughness!(ij, vars, LSR.ocean_surface_roughness, model.land_sea_mask)
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
    vars.parameterizations.moisture_roughness[ij] = land_fraction * land.moisture_roughness[ij] +
        (1 - land_fraction) * ocean.moisture_roughness[ij]
    return nothing
end

@propagate_inbounds function surface_roughness!(ij, vars, scheme::LearnedLandRoughness, land_sea_mask)
    land_fraction = land_sea_mask.land_fraction[ij]
    (; land) = vars.parameterizations

    if land_fraction <= 0
        zero_val = zero(land_fraction)
        land.momentum_roughness[ij] = zero_val
        land.heat_roughness[ij] = zero_val
        land.moisture_roughness[ij] = zero_val
        return nothing
    end

    # Unpack predictors from the land parameterization
    vₕ = vars.parameterizations.land.vegetation_high[ij]
    vₗ = vars.parameterizations.land.vegetation_low[ij]
    laiₕ = vars.parameterizations.land.lai_vegetation_high[ij]
    laiₗ = vars.parameterizations.land.lai_vegetation_low[ij]
    sd = vars.prognostic.land.snow_depth[ij]
    soil_moisture = vars.prognostic.land.soil_moisture[ij, begin]  # currently top layer
    soil_temperature = vars.prognostic.land.soil_temperature[ij, end]  # currently bottom layer

    z₀M_land = land_momentum_roughness(vₕ, vₗ, laiₕ, laiₗ, sd, soil_temperature, soil_moisture)
    z₀H_land = land_heat_roughness(vₕ, vₗ, laiₕ, laiₗ, soil_temperature, soil_moisture)

    land.momentum_roughness[ij] = z₀M_land
    land.heat_roughness[ij] = z₀H_land
    land.moisture_roughness[ij] = z₀H_land # same as heat roughness over land

    return nothing
end

@propagate_inbounds function surface_roughness!(ij, vars, scheme::LearnedOceanRoughness, land_sea_mask)
    land_fraction = land_sea_mask.land_fraction[ij]
    (; ocean) = vars.parameterizations

    if land_fraction >= 1
        zero_val = zero(land_fraction)
        ocean.momentum_roughness[ij] = zero_val
        ocean.heat_roughness[ij] = zero_val
        ocean.moisture_roughness[ij] = zero_val
        return nothing
    end

    siconc = vars.prognostic.ocean.sea_ice_concentration[ij]
    (; surface_wind_speed) = vars.parameterizations
    αᵪ = ocean.charnock_parameter[ij] # in log form

    z₀M_ocean = ocean_momentum_roughness(surface_wind_speed[ij], αᵪ, siconc)
    z₀H_ocean = ocean_heat_roughness(surface_wind_speed[ij], siconc)

    ocean.momentum_roughness[ij] = z₀M_ocean
    ocean.heat_roughness[ij] = z₀H_ocean
    ocean.moisture_roughness[ij] = z₀H_ocean * (scheme.a_q / scheme.a_h) # moisture roughness is heat times a constant

    return nothing
end

@inline function sea_ice_momentum_roughness(siconc::NF) where {NF <: Real}
    # Following IFS documentation
    base = NF(0.93e-3) * (1 - siconc) + NF(6.05e-3) * exp(NF(-17) * (siconc - 0.5)^2)
    return max(base, NF(1.0e-3))
end

@inline function maximum_momentum_roughness(uₙ::NF) where {NF <: Real}
    z_cap = NF(6.74e-3)
    z_high = NF(1.3e-3)
    U_cap = NF(33.0)
    U_high = NF(55.0)

    # The compiler folds this constant arithmetic at compile-time (m ≈ -2.4727e-4)
    m = (z_high - z_cap) / (U_high - U_cap)

    # Evaluate linear regime using FMA: z_cap + m*(uₙ - U_cap)
    z_linear = muladd(m, uₙ - U_cap, z_cap)

    # Branchless clamp replaces max(z_high, min(z_linear, z_cap))
    return clamp(z_linear, z_high, z_cap)
end

@inline function ocean_momentum_roughness(uₙ::NF, αᵪ::NF, siconc::NF) where {NF <: Real}
    c1 = NF(0.48786303)
    c2 = NF(3.374574)
    c3 = NF(0.5460857)
    c4 = NF(0.45002246)
    c5 = NF(0.9827853)
    c6 = NF(0.043734703)
    c7 = NF(1.0346544)
    c8 = NF(7.1816134)

    num = muladd(c3, uₙ, αᵪ + c2)
    den = exp(muladd(c4, uₙ, -c5)) + c6

    log_roughness = (uₙ^c1) - (num / den) + muladd(c7, αᵪ, -c8)
    z_base = exp(log_roughness)

    # Calculate the max roughness
    z_max = maximum_momentum_roughness(uₙ)

    z₀M = min(z_base, z_max)

    # Weight by sea ice concentration
    return sea_ice_momentum_roughness(siconc) * siconc + (1 - siconc) * z₀M
end

@inline function ocean_heat_roughness(uₙ::NF, siconc::NF) where {NF <: Real}
    c1 = NF(-8.207549)
    c2 = NF(0.18300448)
    c3 = NF(0.08970642)
    c4 = NF(-1.6667988)
    log_roughness = (c1 * ((uₙ + c2)^c3)) - exp(c4 / uₙ)
    return NF(1.0e-3) * siconc + (1 - siconc) * exp(log_roughness) # following IFS documentation, ice roughness is 1e-3 m
end

@inline function land_momentum_roughness(C_H::NF, C_L::NF, L_H::NF, L_L::NF, S::NF, T::NF, W::NF) where {NF <: Real}
    # Coefficients for z_0M
    c1 = NF(-7.231684)
    c2 = NF(11.078163)
    c3 = NF(-0.14265412)
    c4 = NF(-5.916258)
    c5 = NF(2.2420337)
    c6 = NF(0.00613743)
    c7 = NF(0.51263654)
    c8 = NF(2.6413066)
    c9 = NF(0.60989547)

    # Numerator terms
    exp_arg1 = S / muladd(c2, C_H, -c3) # c2*C_H - c3
    term1 = c4 * exp(exp_arg1)

    min_inner = min(c5 * L_L, muladd(c6, T, -W)) # min(c5*L_L, c6*T - W)
    exp_arg2 = C_L * (min_inner - (c7^L_L))
    term2 = exp(exp_arg2)

    numerator = max(c1, term1 + term2)

    # Denominator terms
    term3 = exp(C_H * max(L_H, c8))
    term4 = min(L_H, W)

    denominator = term3 + term4

    log_roughness = (numerator / denominator) + c9
    return exp(log_roughness)
end

@inline function land_heat_roughness(C_H::NF, C_L::NF, L_H::NF, L_L::NF, T::NF, W::NF) where {NF <: Real}
    # Coefficients for z_0H
    c1 = NF(373.86444)
    c2 = NF(5.981257)
    c3 = NF(3.5229454)
    c4 = NF(7.781359)
    c5 = NF(2.1035099)
    c6 = NF(-0.41931123)
    c7 = NF(0.00037635124)
    c8 = NF(-5.4191554)

    # Calculate intermediate Gamma term (eq. land_z0H line 1)
    # Note: c1^c2 is folded as a constant if c1 and c2 are hardcoded, but evaluated safely here
    gamma_term1 = ((c1^c2) * min(C_H, L_H)) / ((L_L + T)^c2) + L_H
    gamma_term2 = muladd(NF(2), C_H, c3) # 2*C_H + c3
    Γ = max(gamma_term1, gamma_term2)

    # Calculate z_0H (eq. land_z0H line 2)
    inner_min = min(C_L, L_L - W)
    arg1 = muladd(c5, Γ - c4 + C_L * inner_min, -(c6 * L_H)) # c5*(...) - c6*L_H
    arg2 = c8 * (c7^C_H)

    log_roughness = min(arg1, arg2)
    return exp(log_roughness)
end
