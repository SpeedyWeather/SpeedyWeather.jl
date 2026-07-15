abstract type AbstractSurfaceRoughness <: AbstractParameterization end

export ConstantSurfaceRoughness
export LearnedSurfaceRoughness

"""
Surface roughness length parameterization with constant roughness length
over land and ocean. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct ConstantSurfaceRoughness{NF} <: AbstractSurfaceRoughness
    "[OPTION] constant roughness length over land [m]"
    @param roughness_length_land::NF = 0.5 (bounds = Nonnegative,)

    "[OPTION] constant roughness length over ocean [m]"
    @param roughness_length_ocean::NF = 1.0e-4 (bounds = Nonnegative,)
end

Adapt.@adapt_structure ConstantSurfaceRoughness
ConstantSurfaceRoughness(SG::SpectralGrid; kwargs...) = ConstantSurfaceRoughness{SG.NF}(; kwargs...)
initialize!(::ConstantSurfaceRoughness, ::PrimitiveEquation) = nothing

"""
Learned ocean surface roughness for momentum, heat and moisture fluxes."""
@parameterized @kwdef struct LearnedSurfaceRoughness{NF} <: AbstractSurfaceRoughness
    "[OPTION] Learned roughness length over ocean [m]"
    @param roughness_length_ocean::NF = 1.0e-4 (bounds = Nonnegative,)
    @param roughness_length_land::NF = 0.5 (bounds = Nonnegative,)
end

Adapt.@adapt_structure LearnedSurfaceRoughness
LearnedSurfaceRoughness(SG::SpectralGrid; kwargs...) = LearnedSurfaceRoughness{SG.NF}(; kwargs...)
initialize!(::LearnedSurfaceRoughness, ::PrimitiveEquation) = nothing

@propagate_inbounds parameterization!(ij, vars, scheme::AbstractSurfaceRoughness, model) =
    surface_roughness!(ij, vars, scheme, model.land_sea_mask)

variables(::AbstractSurfaceRoughness) = (
    ParameterizationVariable(:surface_roughness, Grid2D(), desc = "Surface roughness length", units = "m"),
    ParameterizationVariable(:surface_roughness, Grid2D(), desc = "Land surface roughness length", units = "m", namespace = :land),
    ParameterizationVariable(:surface_roughness, Grid2D(), desc = "Ocean surface roughness length", units = "m", namespace = :ocean),

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
    z₀_land = scheme.roughness_length_land
    z₀_ocean = scheme.roughness_length_ocean
    (; land, ocean) = vars.parameterizations

    land.surface_roughness[ij] = ifelse(land_fraction > 0, z₀_land, zero(z₀_land))
    ocean.surface_roughness[ij] = ifelse(land_fraction < 1, z₀_ocean, zero(z₀_ocean))

    vars.parameterizations.surface_roughness[ij] = land_fraction * land.surface_roughness[ij] +
        (1 - land_fraction) * ocean.surface_roughness[ij]
    return nothing
end

@propagate_inbounds function surface_roughness!(ij, vars, scheme::LearnedSurfaceRoughness, land_sea_mask)
    land_fraction = land_sea_mask.land_fraction[ij]
    (; land, ocean, surface_wind_speed) = vars.parameterizations
    αᵪ = ocean.charnock_parameter[ij] # in log form

    z₀_land = scheme.roughness_length_land
    z₀_ocean = ocean_momentum_roughness(surface_wind_speed[ij], αᵪ)

    land.surface_roughness[ij] = ifelse(land_fraction > 0, z₀_land, zero(z₀_land))
    ocean.surface_roughness[ij] = ifelse(land_fraction < 1, z₀_ocean, zero(z₀_ocean))

    vars.parameterizations.surface_roughness[ij] = land_fraction * land.surface_roughness[ij] +
        (1 - land_fraction) * ocean.surface_roughness[ij]
    return nothing
end

"""
    maximum_momentum_roughness(uₙ::T) where {T<:Real}

Calculates the high-wind speed aerodynamic ceiling z_max(uₙ) for momentum roughness.
Uses branchless clamping for GPU efficiency.
"""
@inline function maximum_momentum_roughness(uₙ::T) where {T<:Real}
    z_cap  = T(6.74e-3)
    z_high = T(1.30e-3)
    U_cap  = T(33.0)
    U_high = T(55.0)

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
@inline function ocean_momentum_roughness(uₙ::NF, αᵪ::NF) where {NF <: Real}
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

    return min(z_base, z_max)
end
