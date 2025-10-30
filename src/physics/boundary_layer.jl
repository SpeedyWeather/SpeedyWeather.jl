abstract type AbstractBoundaryLayer <: AbstractParameterization end

# function barrier for all boundary layer drags
function boundary_layer_drag!(  column::ColumnVariables,
                                model::PrimitiveEquation)
    boundary_layer_drag!(column, model.boundary_layer_drag, model)
end

# no boundary layer drag
boundary_layer_drag!(::ColumnVariables, ::Nothing, ::PrimitiveEquation) = nothing

export ConstantDrag
@kwdef struct ConstantDrag{NF} <: AbstractBoundaryLayer
    drag::NF = 1e-3
end

ConstantDrag(SG::SpectralGrid; kwargs...) = ConstantDrag{SG.NF}(; kwargs...)
initialize!(::ConstantDrag, ::PrimitiveEquation) = nothing
function boundary_layer_drag!(  column::ColumnVariables,
                                scheme::ConstantDrag,
                                model::PrimitiveEquation)
    column.boundary_layer_drag = scheme.drag
end

export BulkRichardsonDrag

"""Boundary layer drag coefficient from the bulk Richardson number,
following Frierson, 2006, Journal of the Atmospheric Sciences.
$(TYPEDFIELDS)"""
@kwdef struct BulkRichardsonDrag{NF} <: AbstractBoundaryLayer
    "[OPTION] von Kármán constant [1]"
    κ::NF = 0.4

    "[OPTION] roughness length [m]"
    z₀::NF = 0.01

    "[OPTION] Critical Richardson number for stable mixing cutoff [1]"
    Ri_c::NF = 3

    "[OPTION] Drag minimum to avoid zero surface fluxes in stable conditions [1]"
    drag_min::NF = 1e-5

    "[OPTION] Gust speed to be added to wind speed for drag calculation [m/s]"
    gust_speed::NF = 5

    "[DERIVED] Maximum drag coefficient, κ²/log(zₐ/z₀) for zₐ from reference temperature"
    drag_max::Base.RefValue{NF} = Ref(zero(NF))
end

BulkRichardsonDrag(SG::SpectralGrid; kwargs...) = BulkRichardsonDrag{SG.NF}(; kwargs...)

function initialize!(scheme::BulkRichardsonDrag, model::PrimitiveEquation)

    # Typical height Z of lowermost layer from geopotential of reference surface temperature
    # minus surface geopotential (orography * gravity)
    (; temp_ref) = model.atmosphere
    (; gravity) = model.planet
    (; Δp_geopot_full) = model.geopotential
    Z = temp_ref * Δp_geopot_full[end] / gravity

    # maximum drag Cmax from that height, stable conditions would decrease Cmax towards 0
    # Frierson 2006, eq (12)
    (; κ, z₀) = scheme
    scheme.drag_max[] = (κ/log(Z/z₀))^2
end

function boundary_layer_drag!(  column::ColumnVariables,
                                scheme::BulkRichardsonDrag,
                                model::PrimitiveEquation)
    boundary_layer_drag!(column, scheme, model.atmosphere)
end

"""$(TYPEDSIGNATURES)
Calculate a boundary_layer_drag from the bulk richardson number
following Frierson, 2006."""
function boundary_layer_drag!(
    column::ColumnVariables,
    drag::BulkRichardsonDrag,
    atmosphere::AbstractAtmosphere,
)
    
    (; drag_min, Ri_c) = drag
    drag_max = drag.drag_max[]

    # bulk Richardson number at lowermost layer N from Frierson, 2006, eq. (15)
    # they call it Ri_a = Ri_N here
    Ri_N = bulk_richardson_surface(column, drag, atmosphere)

    # clamp to get the cases, eq (12-14)
    # if Ri_N > Ri_c then C = 0
    # if Ri_c > Ri_N > 0 then = κ^2/log(z_N/z₀)^2 * (1-Ri_N/Ri_c)^2
    # if Ri_c < 0 then κ^2/log(z_N/z₀)^2
    # but Z ≈ z_N given that the lowermost layer height is proportional to temperature
    # which doesn't change much with instantaneous temperature variations but with
    # vertical resolution, hence κ^2/log(Z/z₀)^2 is precalculated in initialize!
    Ri_N = clamp(Ri_N, 0, Ri_c)

    # extension over Frierson, 2006: minimum drag to avoid zero surface fluxes
    column.boundary_layer_drag = max(drag_min, drag_max*(1-Ri_N/Ri_c)^2)
end

"""
$(TYPEDSIGNATURES)
Calculate the bulk richardson number following Frierson, 2007.
For vertical stability in the boundary layer."""
function bulk_richardson_surface(
    column::ColumnVariables,
    drag::BulkRichardsonDrag,
    atmosphere::AbstractAtmosphere,
)
    cₚ = atmosphere.heat_capacity
    (; u, v, geopot, temp_virt, surface_geopotential) = column
    (; gust_speed) = drag
    surface = column.nlayers    # surface index = nlayers

    V² = u[surface]^2 + v[surface]^2 + gust_speed^2
    gz = geopot[surface] - surface_geopotential
    Θ₀ = cₚ*temp_virt[surface]
    Θ₁ = Θ₀ + gz
    bulk_richardson = gz*(Θ₁ - Θ₀) / (Θ₀*V²)
    return bulk_richardson
end

