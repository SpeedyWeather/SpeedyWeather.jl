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
        