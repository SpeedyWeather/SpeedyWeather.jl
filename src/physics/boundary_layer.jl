"""Concrete type that disables the boundary layer scheme."""
struct NoBoundaryLayer{NF} <: BoundaryLayer{NF} end

"""NoBoundaryLayer scheme just passes."""
function boundary_layer!(   column::ColumnVariables,
                            scheme::NoBoundaryLayer,
                            model::PrimitiveEquation)
    return nothing
end

"""NoBoundaryLayer scheme does not need any initialisation."""
function initialize_boundary_layer!(K::ParameterizationConstants,
                                    scheme::NoBoundaryLayer,
                                    P::Parameters,
                                    G::Geometry)
    return nothing
end 

"""Following Held and Suarez, 1996 BAMS"""
Base.@kwdef struct LinearDrag{NF<:Real} <: BoundaryLayer{NF}
    σb::NF = 0.7            # sigma coordinate below which linear drag is applied
    drag_time::NF = 24.0    # [hours] time scale for linear drag coefficient at σ=1 (=1/kf in HS96)
end

# generator so that LinearDrag(drag_time=1::Int) is still possible → Float64
LinearDrag(;kwargs...) = LinearDrag{Float64}(;kwargs...)

function boundary_layer!(   column::ColumnVariables,
                            scheme::LinearDrag,
                            model::PrimitiveEquation)
    (;u,v,u_tend,v_tend) = column
    (;drag_coefs) = model.parameterization_constants

    @inbounds for k in eachlayer(column)
        kᵥ = drag_coefs[k]
        if kᵥ > 0
            u_tend[k] -= kᵥ*u[k]    # Held and Suarez 1996, equation 1
            v_tend[k] -= kᵥ*v[k]
        end
    end
end

function initialize_boundary_layer!(K::ParameterizationConstants,
                                    scheme::LinearDrag,
                                    P::Parameters,
                                    G::Geometry)

    (;σ_levels_full,radius) = G
    (;σb,drag_time) = scheme
    (;drag_coefs) = K

    kf = radius/(drag_time*3600)          # scale with radius as ∂ₜu is; hrs -> sec

    for (k,σ) in enumerate(σ_levels_full)
        drag_coefs[k] = kf*max(0,(σ-σb)/(1-σb)) # drag only below σb, lin increasing to kf at σ=1
    end
end 
        