abstract type AbstractBoundaryLayer <: AbstractParameterization end

# function barrier for all boundary layer drags
function boundary_layer_drag!(  column::ColumnVariables,
                                model::PrimitiveEquation)
    boundary_layer_drag!(column, model.boundary_layer_drag, model)
end

# no boundary layer drag
boundary_layer_drag!(::ColumnVariables, ::Nothing, ::PrimitiveEquation) = nothing

export LinearDrag

"""Linear boundary layer drag following Held and Suarez, 1996 BAMS
$(TYPEDFIELDS)"""
@kwdef struct LinearDrag{NF<:AbstractFloat} <: AbstractBoundaryLayer
    
    # DIMENSIONS
    nlayers::Int
    
    # PARAMETERS
    σb::NF = 0.7                    # sigma coordinate below which linear drag is applied
    time_scale::Second = Hour(24)   # time scale for linear drag coefficient at σ=1 (=1/kf in HS96)

    # PRECOMPUTED CONSTANTS
    drag_coefs::Vector{NF} = zeros(NF, nlayers)
end

"""
$(TYPEDSIGNATURES)
Generator function using `nlayers` from `SG::SpectralGrid`"""
LinearDrag(SG::SpectralGrid; kwargs...) = LinearDrag{SG.NF}(nlayers=SG.nlayers; kwargs...)

"""
$(TYPEDSIGNATURES)
Precomputes the drag coefficients for this `BoundaryLayerDrag` scheme."""
function initialize!(   scheme::LinearDrag,
                        model::PrimitiveEquation)

    (; σ_levels_full) = model.geometry
    (; σb, time_scale, drag_coefs) = scheme

    kf = 1/time_scale.value

    for (k, σ) in enumerate(σ_levels_full)
        drag_coefs[k] = -kf*max(0, (σ-σb)/(1-σb))    # drag only below σb, lin increasing to kf at σ=1
    end
end 

# function barrier
function boundary_layer_drag!(  column::ColumnVariables,
                                scheme::LinearDrag,
                                model::PrimitiveEquation)
    boundary_layer_drag!(column, scheme)
end

"""
$(TYPEDSIGNATURES)
Compute tendency for boundary layer drag of a `column` and add to its tendencies fields"""
function boundary_layer_drag!(  column::ColumnVariables,
                                scheme::LinearDrag)

    (; u, v, u_tend, v_tend) = column
    (; drag_coefs) = scheme

    @inbounds for k in eachlayer(column)
        kᵥ = drag_coefs[k]
        if kᵥ > 0
            u_tend[k] += kᵥ*u[k]    # Held and Suarez 1996, equation 1
            v_tend[k] += kᵥ*v[k]
        end
    end
end

export ConstantDrag
Base.@kwdef struct ConstantDrag{NF} <: AbstractBoundaryLayer
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

function boundary_layer_drag!(
    column::ColumnVariables,
    scheme::BulkRichardsonDrag,
    atmopshere::AbstractAtmosphere,
)
    
    (; Ri_c) = scheme
    drag_max = scheme.drag_max[]

    # bulk Richardson number at lowermost layer N from Frierson, 2006, eq. (15)
    # they call it Ri_a = Ri_N here
    Ri_N = bulk_richardson_surface(column, atmopshere)

    # clamp to get the cases, eq (12-14)
    # if Ri_N > Ri_c then C = 0
    # if Ri_c > Ri_N > 0 then = κ^2/log(z_N/z₀)^2 * (1-Ri_N/Ri_c)^2
    # if Ri_c < 0 then κ^2/log(z_N/z₀)^2
    # but Z ≈ z_N given that the lowermost layer height is proportional to temperature
    # which doesn't change much with instantaneous temperature variations but with
    # vertical resolution, hence κ^2/log(Z/z₀)^2 is precalculated in initialize!
    Ri_N = clamp(Ri_N, 0, Ri_c)
    column.boundary_layer_drag = drag_max*(1-Ri_N/Ri_c)^2
end

"""
$(TYPEDSIGNATURES)
Calculate the bulk richardson number following Frierson, 2007.
For vertical stability in the boundary layer."""
function bulk_richardson_surface(
    column::ColumnVariables,
    atmosphere::AbstractAtmosphere,
)
    cₚ = atmosphere.heat_capacity
    (; u, v, geopot, temp_virt) = column
    surface = column.nlayers    # surface index = nlayers

    V² = u[surface]^2 + v[surface]^2
    Θ₀ = cₚ*temp_virt[surface]
    Θ₁ = Θ₀ + geopot[surface]
    bulk_richardson = geopot[surface]*(Θ₁ - Θ₀) / (Θ₀*V²)
    return bulk_richardson
end

