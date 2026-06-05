abstract type AbstractBoundaryLayer <: AbstractParameterization end

variables(::AbstractBoundaryLayer) = (
    ParameterizationVariable(:boundary_layer_drag, Grid2D(), desc = "Boundary layer drag coefficient", units = "1"),
)

export ConstantDrag
"""Constant boundary layer drag coefficient. Fields are $(TYPEDFIELDS)"""
@kwdef struct ConstantDrag{NF} <: AbstractBoundaryLayer
    "[OPTION] Constant drag coefficient [1]"
    drag::NF = 1.0e-3
end

Adapt.@adapt_structure ConstantDrag
ConstantDrag(SG::SpectralGrid; kwargs...) = ConstantDrag{SG.NF}(; kwargs...)
initialize!(::ConstantDrag, ::PrimitiveEquation) = nothing
@propagate_inbounds parameterization!(ij, diagn, progn, drag::ConstantDrag, model) =
    boundary_layer_drag!(ij, diagn, drag)

@propagate_inbounds function boundary_layer_drag!(ij, diagn, scheme::ConstantDrag)
    return diagn.physics.boundary_layer_drag[ij] = scheme.drag
end

abstract type AbstractSurfaceRoughness <: AbstractParameterization end

export ConstantSurfaceRoughness

"""
Surface roughness length parameterization with constant roughness length
over land and ocean. Fields are $(TYPEDFIELDS)"""
@kwdef struct ConstantSurfaceRoughness{NF} <: AbstractSurfaceRoughness
    "[OPTION] constant roughness length over land [m]"
    roughness_length_land::NF = 0.5

    "[OPTION] constant roughness length over ocean [m]"
    roughness_length_ocean::NF = 1.0e-4
end

Adapt.@adapt_structure ConstantSurfaceRoughness
ConstantSurfaceRoughness(SG::SpectralGrid, kwargs...) = ConstantSurfaceRoughness{SG.NF}(; kwargs...)
initialize!(::ConstantSurfaceRoughness, ::PrimitiveEquation) = nothing

@propagate_inbounds parameterization!(ij, vars, scheme::AbstractSurfaceRoughness, model) = 
    surface_roughness!(ij, vars, scheme, model.land_sea_mask)

variables(::AbstractSurfaceRoughness) = (
    ParameterizationVariable(name = :surface_roughness, dims = Grid2D(), desc = "Surface roughness length", units = "m"),
    ParameterizationVariable(name = :surface_roughness, dims = Grid2D(), desc = "Land surface roughness length", units = "m", namespace = :land),
    ParameterizationVariable(name = :surface_roughness, dims = Grid2D(), desc = "Ocean surface roughness length", units = "m", namespace = :ocean),
)

@propagate_inbounds function surface_roughness!(ij, vars, scheme::ConstantSurfaceRoughness, land_sea_mask)
    land_fraction = land_sea_mask.mask[ij]

    vars.parameterizations.land.surface_roughness[ij] = ifelse(land_fraction > 0, scheme.roughness_length_land, zero(land_fraction))    
    vars.parameterizations.ocean.surface_roughness[ij] = ifelse(land_fraction < 1, scheme.roughness_length_ocean, zero(land_fraction))  

    vars.parameterizations.surface_roughness[ij] = land_fraction * vars.parameterizations.land.surface_roughness[ij] + (1 - land_fraction) * vars.parameterizations.ocean.surface_roughness[ij]
    return nothing
end

function Base.show(io::IO, scheme::ConstantSurfaceRoughness{NF}) where NF
    print(io, "ConstantSurfaceRoughness{$NF}")
    println(io)
    println(io, "├ Ocean: $(scheme.roughness_length_ocean) m")
    return println(io, "└ Land: $(scheme.roughness_length_land) m")
end

export BulkRichardsonDrag

"""Boundary layer drag coefficient from the bulk Richardson number,
following Frierson, 2006, Journal of the Atmospheric Sciences.
$(TYPEDFIELDS)"""
@kwdef struct BulkRichardsonDrag{NF} <: AbstractBoundaryLayer
    "[OPTION] von Kármán constant [1]"
    von_Karman::NF = 0.4

    "[OPTION] Critical Richardson number for stable mixing cutoff [1]"
    critical_Richardson::NF = 10

    "[OPTION] Drag minimum to avoid zero surface fluxes in stable conditions [1]"
    drag_min::NF = 1.0e-5
end

Adapt.@adapt_structure BulkRichardsonDrag
BulkRichardsonDrag(SG::SpectralGrid, kwargs...) = BulkRichardsonDrag{SG.NF}(; kwargs...)
initialize!(::BulkRichardsonDrag, ::PrimitiveEquation) = nothing

# function barrier
@propagate_inbounds parameterization!(ij, vars, drag::BulkRichardsonDrag, model) =
    boundary_layer_drag!(ij, vars, drag, model.atmosphere, model.planet, model.orography)

@propagate_inbounds function boundary_layer_drag!(
        ij,
        vars,
        drag::BulkRichardsonDrag,
        atmosphere,
        planet,
        orography,
    )

    # Height z [m] of lowermost layer above ground
    surface = size(vars.grid.geopotential, 2)    # surface index = nlayers
    (; gravity) = planet
    z = vars.grid.geopotential[ij, surface] / gravity - orography.orography[ij]

    # maximum drag Cmax from that height, stable conditions would decrease Cmax towards 0
    # Frierson 2006, eq (12)
    κ = drag.von_Karman

    # Get surface roughness length (computed by the surface_roughness parameterization)
    z₀ = vars.parameterizations.surface_roughness[ij]

    # should be z > z₀, z=z₀ means an infinitely high drag
    # 0 < z < z₀ doesn't make sense so cap here
    z = max(z, z₀)
    drag_max = (κ / log(z / z₀))^2

    # bulk Richardson number at lowermost layer from Frierson, 2006, eq. (15)
    # they call it Ri_a = Ri here
    ΔΦ₀ = gravity * z     # geopotential high relative to surface
    Ri = bulk_richardson_surface(ij, ΔΦ₀, vars, atmosphere)
    Ri_c = drag.critical_Richardson
    (; drag_min) = drag

    # clamp to get the cases, eq (12-14)
    # if Ri > Ri_c then C = 0
    # if Ri_c > Ri > 0 then = κ^2/log(z/z₀)^2 * (1-Ri/Ri_c)^2
    # if Ri_c < 0 then κ^2/log(z/z₀)^2
    Ri = clamp(Ri, 0, Ri_c)
    vars.parameterizations.boundary_layer_drag[ij] = max(drag_min, drag_max * (1 - Ri / Ri_c)^2)
    return nothing
end

"""
$(TYPEDSIGNATURES)
Calculate the bulk Richardson number following Frierson, 2006.
For vertical stability in the boundary layer."""
@propagate_inbounds function bulk_richardson_surface(ij, ΔΦ₀, vars, atmosphere)
    cₚ = atmosphere.heat_capacity
    NF = eltype(vars.grid.temperature_prev)
    surface = size(vars.grid.temperature_prev, 2)     # surface index = nlayers

    Vₛ = vars.parameterizations.surface_wind_speed[ij]
    T = vars.grid.temperature_prev[ij, surface]
    q = haskey(vars.grid, :humidity_prev) ? vars.grid.humidity_prev[ij, surface] : zero(NF)
    Tᵥ = virtual_temperature(T, q, atmosphere)

    # bulk Richardson number at lowermost layer N from Frierson, 2006, eq. (15)
    Θ₀ = cₚ * Tᵥ          # virtual dry static energy at surface (z=0)
    Θ₁ = Θ₀ + ΔΦ₀       # virtual dry static energy at first model level (z=z)
    bulk_richardson = ΔΦ₀ * (Θ₁ - Θ₀) / (Θ₀ * Vₛ^2)
    return bulk_richardson
end
