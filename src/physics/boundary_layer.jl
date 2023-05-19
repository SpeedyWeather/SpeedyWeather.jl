"""Concrete type that disables the boundary layer scheme."""
struct NoBoundaryLayer{NF} <: BoundaryLayer{NF} end

"""NoBoundaryLayer scheme just passes."""
function boundary_layer!(   column::ColumnVariables,
                            scheme::NoBoundaryLayer,
                            model::PrimitiveEquation)
    return nothing
end

"""NoBoundaryLayer scheme does not need any initialisation."""
function initialize!(   scheme::NoBoundaryLayer,
                        model::PrimitiveEquation)
    return nothing
end 

"""Following Held and Suarez, 1996 BAMS"""
@kwdef struct LinearDrag{NF<:Real} <: BoundaryLayer{NF}
    # PARAMETERS
    σb::Float64 = 0.7           # sigma coordinate below which linear drag is applied
    time_scale::Float64 = 24    # [hours] time scale for linear drag coefficient at σ=1 (=1/kf in HS96)

    # PRECOMPUTED CONSTANTS
    nlev::Int = 0
    drag_coefs::Vector{NF} = zeros(NF,nlev)
end

LinearDrag(geospectral::Geospectral{NF},kwargs...) where NF = LinearDrag{NF}(nlev=geospectral.nlev;kwargs...)

function initialize!(   scheme::LinearDrag,
                        model::PrimitiveEquation)

    (;σ_levels_full,radius) = model.geometry
    (;σb,time_scale,drag_coefs) = scheme

    kf = radius/(time_scale*3600)               # scale with radius as ∂ₜu is; hrs -> sec

    for (k,σ) in enumerate(σ_levels_full)
        drag_coefs[k] = kf*max(0,(σ-σb)/(1-σb)) # drag only below σb, lin increasing to kf at σ=1
    end
end 
        
function boundary_layer!(   column::ColumnVariables,
                            scheme::LinearDrag,
                            model::PrimitiveEquation)

    (;u,v,u_tend,v_tend) = column
    (;drag_coefs) = scheme

    @inbounds for k in eachlayer(column)
        kᵥ = drag_coefs[k]
        if kᵥ > 0
            u_tend[k] -= kᵥ*u[k]    # Held and Suarez 1996, equation 1
            v_tend[k] -= kᵥ*v[k]
        end
    end
end