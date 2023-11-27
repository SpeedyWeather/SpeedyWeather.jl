"""Concrete type that disables the boundary layer drag scheme."""
struct NoBoundaryLayerDrag{NF} <: BoundaryLayerDrag{NF} end

NoBoundaryLayerDrag(SG::SpectralGrid) = NoBoundaryLayerDrag{SG.NF}()

"""NoBoundaryLayer scheme does not need any initialisation."""
function initialize!(   scheme::NoBoundaryLayerDrag,
                        model::PrimitiveEquation)
    return nothing
end 

# function barrier
function boundary_layer_drag!(  column::ColumnVariables,
                                model::PrimitiveEquation)
    boundary_layer_drag!(column,model.boundary_layer_drag)
end

"""NoBoundaryLayer scheme just passes."""
function boundary_layer_drag!(  column::ColumnVariables,
                                scheme::NoBoundaryLayerDrag)
    return nothing
end

"""Linear boundary layer drag Following Held and Suarez, 1996 BAMS
$(TYPEDFIELDS)"""
Base.@kwdef struct LinearDrag{NF<:AbstractFloat} <: BoundaryLayerDrag{NF}
    # PARAMETERS
    σb::Float64 = 0.7           # sigma coordinate below which linear drag is applied
    time_scale::Float64 = 24    # [hours] time scale for linear drag coefficient at σ=1 (=1/kf in HS96)

    # PRECOMPUTED CONSTANTS
    nlev::Int = 0
    drag_coefs::Vector{NF} = zeros(NF,nlev)
end

"""
$(TYPEDSIGNATURES)
Generator function using `nlev` from `SG::SpectralGrid`"""
LinearDrag(SG::SpectralGrid;kwargs...) = LinearDrag{SG.NF}(nlev=SG.nlev;kwargs...)

"""
$(TYPEDSIGNATURES)
Precomputes the drag coefficients for this `BoundaryLayerDrag` scheme."""
function initialize!(   scheme::LinearDrag,
                        model::PrimitiveEquation)

    (;σ_levels_full) = model.geometry
    (;σb,time_scale,drag_coefs) = scheme

    kf = 1/(time_scale*3600)                        # scale with radius as ∂ₜu is; hrs -> sec

    for (k,σ) in enumerate(σ_levels_full)
        drag_coefs[k] = -kf*max(0,(σ-σb)/(1-σb))    # drag only below σb, lin increasing to kf at σ=1
    end
end 

"""
$(TYPEDSIGNATURES)
Compute tendency for boundary layer drag of a `column` and add to its tendencies fields"""
function boundary_layer_drag!(  column::ColumnVariables,
                                scheme::LinearDrag)

    (;u,v,u_tend,v_tend) = column
    (;drag_coefs) = scheme

    @inbounds for k in eachlayer(column)
        kᵥ = drag_coefs[k]
        if kᵥ > 0
            u_tend[k] += kᵥ*u[k]    # Held and Suarez 1996, equation 1
            v_tend[k] += kᵥ*v[k]
        end
    end
end