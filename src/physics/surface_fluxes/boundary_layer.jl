abstract type AbstractBoundaryLayer <: AbstractParameterization end

variables(::AbstractBoundaryLayer) = (
    DiagnosticVariable(name=:boundary_layer_drag, dims=Grid2D(), desc="Boundary layer drag coefficient", units="1"),
    )

export ConstantDrag
@kwdef struct ConstantDrag{NF} <: AbstractBoundaryLayer
    drag::NF = 1e-3
end

Adapt.@adapt_structure ConstantDrag

ConstantDrag(SG::SpectralGrid; kwargs...) = ConstantDrag{SG.NF}(; kwargs...)
initialize!(::ConstantDrag, ::PrimitiveEquation) = nothing
parameterization!(ij, diagn, drag::ConstantDrag, model) = boundary_layer_drag!(ij, diagn, drag)
function boundary_layer_drag!(ij, diagn, scheme::ConstantDrag)
    diagn.physics.boundary_layer_drag[ij] = scheme.drag
end

export BulkRichardsonDrag

"""Boundary layer drag coefficient from the bulk Richardson number,
following Frierson, 2006, Journal of the Atmospheric Sciences.
$(TYPEDFIELDS)"""
@kwdef struct BulkRichardsonDrag{NF} <: AbstractBoundaryLayer
    "von Kármán constant [1]"
    κ::NF = 0.4

    "roughness length [m]"
    z₀::NF = 3.21e-5

    "Critical Richardson number for stable mixing cutoff [1]"
    Ri_c::NF = 10

    "Maximum drag coefficient, κ²/log(zₐ/z₀) for zₐ from reference temperature"
    drag_max::Base.RefValue{NF} = Ref(zero(NF))
end

BulkRichardsonDrag(SG::SpectralGrid, kwargs...) = BulkRichardsonDrag{SG.NF}(; kwargs...)

function initialize!(scheme::BulkRichardsonDrag, model::PrimitiveEquation)

    # Typical height Z of lowermost layer from geopotential of reference surface temperature
    # minus surface geopotential (orography * gravity)
    (; temp_ref) = model.atmosphere
    (; gravity) = model.planet
    (; Δp_geopot_full) = model.geopotential
    GPUArrays.@allowscalar Z = temp_ref * Δp_geopot_full[end] / gravity

    # maximum drag Cmax from that height, stable conditions would decrease Cmax towards 0
    # Frierson 2006, eq (12)
    (; κ, z₀) = scheme
    scheme.drag_max[] = (κ/log(Z/z₀))^2
end

# function barrier
function parameterization!(ij, diagn, drag::BulkRichardsonDrag, model)
    boundary_layer_drag!(ij, diagn, drag, model.atmosphere)
end

function boundary_layer_drag!(
    ij,
    diagn,
    drag::BulkRichardsonDrag,
    atmosphere::AbstractAtmosphere,
)
    (; Ri_c) = drag
    drag_max = drag.drag_max[]

    # bulk Richardson number at lowermost layer N from Frierson, 2006, eq. (15)
    # they call it Ri_a = Ri_N here
    Ri_N = bulk_richardson_surface(ij, diagn, atmosphere)

    # clamp to get the cases, eq (12-14)
    # if Ri_N > Ri_c then C = 0
    # if Ri_c > Ri_N > 0 then = κ^2/log(z_N/z₀)^2 * (1-Ri_N/Ri_c)^2
    # if Ri_c < 0 then κ^2/log(z_N/z₀)^2
    # but Z ≈ z_N given that the lowermost layer height is proportional to temperature
    # which doesn't change much with instantaneous temperature variations but with
    # vertical resolution, hence κ^2/log(Z/z₀)^2 is precalculated in initialize!
    Ri_N = clamp(Ri_N, 0, Ri_c)
    diagn.physics.boundary_layer_drag[ij] = drag_max*(1-Ri_N/Ri_c)^2
end

"""
$(TYPEDSIGNATURES)
Calculate the bulk richardson number following Frierson, 2007.
For vertical stability in the boundary layer."""
function bulk_richardson_surface(ij, diagn, atmosphere)
    cₚ = atmosphere.heat_capacity
    surface = diagn.grid.nlayers    # surface index = nlayers

    u = diagn.grid.u_grid_prev[ij, surface]
    v = diagn.grid.v_grid_prev[ij, surface]
    geopot = diagn.grid.geopotential[ij, surface]
    temp_virt = diagn.grid.temp_virt_grid[ij, surface]

    V² = u^2 + v^2
    Θ₀ = cₚ*temp_virt
    Θ₁ = Θ₀ + geopot
    bulk_richardson = geopot*(Θ₁ - Θ₀) / (Θ₀*V²)
    return bulk_richardson
end