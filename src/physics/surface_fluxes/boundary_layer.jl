abstract type AbstractBoundaryLayer <: AbstractParameterization end

variables(::AbstractBoundaryLayer) = (
    DiagnosticVariable(name=:boundary_layer_drag, dims=Grid2D(), desc="Boundary layer drag coefficient", units="1"),
    )

export ConstantDrag
"""Constant boundary layer drag coefficient. Fields are $(TYPEDFIELDS)"""
@kwdef struct ConstantDrag{NF} <: AbstractBoundaryLayer
    "[OPTION] Constant drag coefficient [1]"
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
@kwdef struct BulkRichardsonDrag{NF, RefNF <: Ref{NF}} <: AbstractBoundaryLayer
    "[OPTION] von Kármán constant [1]"
    von_Karman::NF = 0.4

    "[OPTION] roughness length [m]"
    roughness_length::NF = 3.21e-5

    "[OPTION] Critical Richardson number for stable mixing cutoff [1]"
    critical_Richardson::NF = 10
end

Adapt.@adapt_structure BulkRichardsonDrag
BulkRichardsonDrag(SG::SpectralGrid, kwargs...) = BulkRichardsonDrag{SG.NF}(; kwargs...)
initialize!(::BulkRichardsonDrag, ::PrimitiveEquation) = nothing

# function barrier
function parameterization!(ij, diagn, progn, drag::BulkRichardsonDrag, model)
    boundary_layer_drag!(ij, diagn, drag, model.atmosphere)
end

function boundary_layer_drag!(
    ij,
    diagn,
    drag::BulkRichardsonDrag,
    atmosphere::AbstractAtmosphere,
    planet::AbstractPlanet,
    geopotential::AbstractGeopotential
)
    # Typical height Z of lowermost layer from geopotential of reference surface temperature
    # minus surface geopotential (orography * gravity)
    (; temp_ref) = atmosphere
    (; gravity) = planet
    (; Δp_geopot_full) = geopotential
    Z = temp_ref * Δp_geopot_full[end] / gravity

    # maximum drag Cmax from that height, stable conditions would decrease Cmax towards 0
    # Frierson 2006, eq (12)
    κ = drag.von_Karman
    z₀ = drag.roughness_length
    drag_max = (κ/log(Z/z₀))^2

    # bulk Richardson number at lowermost layer N from Frierson, 2006, eq. (15)
    # they call it Ri_a = Ri_N here
    Ri_N = bulk_richardson_surface(ij, diagn, atmosphere)
    Ri_c = drag.critical_Richardson

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
@inline function bulk_richardson_surface(ij, diagn, atmosphere)
    cₚ = atmosphere.heat_capacity
    surface = diagn.grid.nlayers    # surface index = nlayers

    @inbounds begin
        u = diagn.grid.u_grid_prev[ij, surface]
        v = diagn.grid.v_grid_prev[ij, surface]
        Φ = diagn.grid.geopotential[ij, surface]
        T = diagn.grid.temp_grid_prev[ij, surface]
        q = diagn.grid.humid_grid_prev[ij, surface]
        Tᵥ = virtual_temperature(T, q, atmosphere)
    end

    V² = u^2 + v^2
    Θ₀ = cₚ*Tᵥ
    Θ₁ = Θ₀ + Φ
    bulk_richardson = Φ*(Θ₁ - Θ₀) / (Θ₀*V²)
    return bulk_richardson
end