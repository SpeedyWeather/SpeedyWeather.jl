abstract type AbstractBoundaryLayer <: AbstractParameterization end

variables(::AbstractBoundaryLayer) = (
    ParameterizationVariable(:boundary_layer_drag, Grid2D(), desc = "Boundary layer drag coefficient", units = "1"),
)

export ConstantDrag
"""Constant boundary layer drag coefficient. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct ConstantDrag{NF} <: AbstractBoundaryLayer
    "[OPTION] Constant drag coefficient [1]"
    @param drag::NF = 1.0e-3
end

Adapt.@adapt_structure ConstantDrag
ConstantDrag(SG::SpectralGrid; kwargs...) = ConstantDrag{SG.NF}(; kwargs...)
initialize!(::ConstantDrag, ::PrimitiveEquation) = nothing
@propagate_inbounds parameterization!(ij, vars, drag::ConstantDrag, model) =
    boundary_layer_drag!(ij, vars, drag)

@propagate_inbounds function boundary_layer_drag!(ij, vars, scheme::ConstantDrag)
    vars.parameterizations.boundary_layer_drag[ij] = scheme.drag
end

export BoundaryLayer
"""Composite type, containing surface roughness computation
and drag coefficient computation. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct BoundaryLayer{SR, D} <: AbstractBoundaryLayer
    @component surface_roughness::SR
    @component drag::D
end

Adapt.@adapt_structure BoundaryLayer
function BoundaryLayer(SG::SpectralGrid;
        surface_roughness = ConstantSurfaceRoughness(SG),
        drag = BulkRichardsonDrag(SG),
    )
    return BoundaryLayer(surface_roughness, drag)
end

function initialize!(BL::BoundaryLayer)
    initialize!(BL.surface_roughness)
    initialize!(BL.drag)
end

# variables of boundary layer are the union of surface roughness and drag variables
variables(BL::BoundaryLayer) = (
    variables(BL.surface_roughness)...,
    variables(BL.drag)...,
)

# just call the sub-paramterizations one after another
@propagate_inbounds function parameterization!(ij, vars, BL::BoundaryLayer, model)
    parameterization!(ij, vars, BL.surface_roughness, model)
    parameterization!(ij, vars, BL.drag, model)
end

export BulkRichardsonDrag

"""Boundary layer drag coefficient from the bulk Richardson number,
following Frierson, 2006, Journal of the Atmospheric Sciences.
$(TYPEDFIELDS)"""
@parameterized @kwdef struct BulkRichardsonDrag{NF} <: AbstractBoundaryLayer
    "[OPTION] von Kármán constant [1]"
    @param von_Karman::NF = 0.4 (bounds = 0 .. 1,)

    "[OPTION] Critical Richardson number for stable mixing cutoff [1]"
    @param critical_Richardson::NF = 10 (bounds = Positive,)

    "[OPTION] Drag minimum to avoid zero surface fluxes in stable conditions [1]"
    @param drag_min::NF = 1.0e-5 (bounds = Nonnegative,)
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
