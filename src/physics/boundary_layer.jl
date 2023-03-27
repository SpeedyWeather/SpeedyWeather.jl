function boundary_layer!(   column::ColumnVariables,
                            scheme::NoBoundaryLayer,
                            model::PrimitiveEquation)
    return nothing
end

function boundary_layer!(   column::ColumnVariables{NF},
                            scheme::LinearDrag,
                            model::PrimitiveEquation) where NF

    (;σ_levels_full,radius_earth) = model.geometry
    (;σb,drag_time) = scheme
    kf = radius_earth/(drag_time*3600)
    (;u,v,u_tend,v_tend) = column

    @inbounds for k in eachlayer(column)

        σ = σ_levels_full[k]            # σ level of layer k
        kᵥ = kf*max(0,(σ-σb)/(1-σb))    # drag only below σb, linearly increasing to kf at σ=1
        
        u_tend[k] -= kᵥ*u[k]
        v_tend[k] -= kᵥ*v[k]
    end
end