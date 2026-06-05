abstract type AbstractSurfaceRoughness <: AbstractParameterization end

export ConstantSurfaceRoughness

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

@propagate_inbounds parameterization!(ij, vars, scheme::AbstractSurfaceRoughness, model) = 
    surface_roughness!(ij, vars, scheme, model.land_sea_mask)

variables(::AbstractSurfaceRoughness) = (
    ParameterizationVariable(:surface_roughness, Grid2D(), desc = "Surface roughness length", units = "m"),
    ParameterizationVariable(:surface_roughness, Grid2D(), desc = "Land surface roughness length", units = "m", namespace = :land),
    ParameterizationVariable(:surface_roughness, Grid2D(), desc = "Ocean surface roughness length", units = "m", namespace = :ocean),
)

@propagate_inbounds function surface_roughness!(ij, vars, scheme::ConstantSurfaceRoughness, land_sea_mask)
    land_fraction = land_sea_mask.mask[ij]
    z₀_land = scheme.roughness_length_land
    z₀_ocean = scheme.roughness_length_ocean
    (; land, ocean) = vars.parameterizations
    
    land.surface_roughness[ij] = ifelse(land_fraction > 0, z₀_land, zero(z₀_land))      
    ocean.surface_roughness[ij] = ifelse(land_fraction < 1, z₀_ocean, zero(z₀_ocean))  

    vars.parameterizations.surface_roughness[ij] = land_fraction * land.surface_roughness[ij] +
                                                (1 - land_fraction) * ocean.surface_roughness[ij]
    return nothing
end