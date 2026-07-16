abstract type AbstractSurfaceRoughness <: AbstractParameterization end

export ConstantSurfaceRoughness
export LearnedSurfaceRoughness

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

"""
Learned ocean surface roughness for momentum, heat and moisture fluxes."""
@kwdef struct LearnedSurfaceRoughness{NF} <: AbstractSurfaceRoughness
    a_h::NF = NF(0.4) # heat related constant
    a_q::NF = NF(0.62) # moisture related constant

    # temporary
    momentum_roughness_land::NF = NF(0.5)
    heat_roughness_land::NF = NF(0.5)
end

Adapt.@adapt_structure LearnedSurfaceRoughness
LearnedSurfaceRoughness(SG::SpectralGrid; kwargs...) = LearnedSurfaceRoughness{SG.NF}(; kwargs...)
initialize!(::LearnedSurfaceRoughness, ::PrimitiveEquation) = nothing

@propagate_inbounds parameterization!(ij, vars, scheme::AbstractSurfaceRoughness, model) =
    surface_roughness!(ij, vars, scheme, model.land_sea_mask)

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

@propagate_inbounds function surface_roughness!(ij, vars, scheme::LearnedSurfaceRoughness, land_sea_mask)
    land_fraction = land_sea_mask.land_fraction[ij]
    siconc = vars.prognostic.ocean.sea_ice_concentration[ij] # sea ice concentration [0-1]
    (; land, ocean, surface_wind_speed) = vars.parameterizations
    αᵪ = ocean.charnock_parameter[ij] # in log form

    z₀M_land = scheme.momentum_roughness_land
    z₀M_ocean = ocean_momentum_roughness(surface_wind_speed[ij], αᵪ, siconc)
    z₀H_land = scheme.heat_roughness_land
    z₀H_ocean = ocean_heat_roughness(surface_wind_speed[ij], siconc)

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

    # Moisture roughness
    land.moisture_roughness[ij] = ifelse(land_fraction > 0, z₀H_land, zero(z₀H_land)) # same as heat roughness for land
    ocean.moisture_roughness[ij] = ifelse(land_fraction < 1, z₀H_ocean * (scheme.a_q / scheme.a_h), zero(z₀H_ocean)) # (just heat times a constant)
    vars.parameterizations.moisture_roughness[ij] = land_fraction * land.moisture_roughness[ij] +
        (1 - land_fraction) * ocean.moisture_roughness[ij]

    return nothing
end

@inline function sea_ice_momentum_roughness(siconc::NF) where {NF <: Real}
    # Following IFS documentation
    base = NF(0.93e-3) * (1 - siconc) + NF(6.05e-3) * exp(NF(-17) * (siconc - 0.5)^2)
    return max(base, NF(1.0e-3))
end

"""
    maximum_momentum_roughness(uₙ::NF) where {NF<:Real}

Calculates the high-wind speed aerodynamic ceiling z_max(uₙ) for momentum roughness.
Uses branchless clamping for GPU efficiency.
"""
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

"""
    ocean_momentum_roughness(uₙ::NF, αᵪ::NF) where {NF<:Real}

Evaluates the symbolic regression closure for ocean momentum roughness z₀M.
"""
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

"""
    ocean_heat_roughness(uₙ::NF) where {NF <: Real}

Evaluates the symbolic regression closure for ocean heat roughness z₀H.
"""
@inline function ocean_heat_roughness(uₙ::NF, siconc::NF) where {NF <: Real}
    c1 = NF(-8.207549)
    c2 = NF(0.18300448)
    c3 = NF(0.08970642)
    c4 = NF(-1.6667988)
    log_roughness = (c1 * ((uₙ + c2)^c3)) - exp(c4 / uₙ)
    return 1.0e-3 * siconc + (1 - siconc) * exp(log_roughness) # following IFS documentation, ice roughness is 1e-3 m
end

@inline function land_momentum_roughness()

end

@inline function land_heat_roughness()

end