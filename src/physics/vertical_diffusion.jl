struct NoVerticalDiffusion{NF} <: VerticalDiffusion{NF} end
NoVerticalDiffusion() = NoVerticalDiffusion{DEFAULT_NF}()

function vertical_diffusion!(   column::ColumnVariables,
                                scheme::NoVerticalDiffusion,
                                model::PrimitiveEquation)
    return nothing
end

function initialize_vertical_diffusion!(K::ParameterizationConstants,
                                        scheme::NoVerticalDiffusion,
                                        P::Parameters,
                                        G::Geometry)
    return nothing
end 

Base.@kwdef struct VerticalLaplacian{NF<:Real} <: VerticalDiffusion{NF}
    time_scale::NF = 6.0        # [hours] time scale to control the strength of vertical diffusion
    height_scale::NF = 100.0    # [m] scales for Δσ so that time_scale is sensible 
    umax::NF = 80.0             # velocity scale [m/s] to scale the diffusion with |u|/umax
    ΔTmax::NF = 50.0            # temperature gradient scale [K] (divided by Δσ)

    scale_with_temperature_gradient::Bool = true
    scale_with_speed::Bool = false
    pressure_coordinates::Bool = true
end

# generator so that arguments are converted to Float64
VerticalLaplacian(;kwargs...) = VerticalLaplacian{Float64}(;kwargs...)

function vertical_diffusion!(   column::ColumnVariables{NF},
                                scheme::VerticalLaplacian,
                                model::PrimitiveEquation) where NF

    (;nlev,u_tend,v_tend,temp_tend) = column
    (;u,v,temp) = column
    ∇²_below = model.parameterization_constants.vert_diff_∇²_below
    ∇²_above = model.parameterization_constants.vert_diff_∇²_above
    Δσ = model.parameterization_constants.vert_diff_Δσ
    ΔTmax = convert(NF,scheme.ΔTmax)

    # GET DIFFUSION COEFFICIENT as a function of u,v,temp and surface pressure
    # *3600 for [hrs] → [s], *1e3 for [km] → [m]
    ν0 = model.geometry.radius*inv(scheme.time_scale*3600)

    # include a height scale, technically not needed, but so that the dimensionless
    # 1/Δσ² gets a resonable scale in meters such that the time scale is not
    # counterintuitively in seconds or years
    ν0 /= scheme.height_scale^2
    
    # scale with surface pressure, on mountains pₛ can be 500hPa, which would increase ν by *2
    scheme.pressure_coordinates && (ν0 *= (model.parameters.pres_ref*100/column.pres[end])^2)
    ν0 = convert(NF,ν0)

    # temperature gradient-based scaling ∈ [0,1] or just pass ν = ν0
    get_ν(ν0,ΔT,ΔTmax,Δσ) = scheme.scale_with_temperature_gradient ? ν0*max(ΔT/Δσ/ΔTmax + 1,0)/2 : ν0

    # DO DIFFUSION
    @inbounds begin
        
        # top layer with no flux boundary conditions at k=1/2
        ΔT = temp[2] - temp[1]          # top temperature gradient
        ν = get_ν(ν0,ΔT,ΔTmax,Δσ[1])    # calculate diffusion coefficient ν

        ν∇²_below = ν*∇²_below[1]       # diffusion operator
        u_tend[1] += ν∇²_below*(u[2] - u[1])    # diffusion of u
        v_tend[1] += ν∇²_below*(v[2] - v[1])    # diffusion of v
        temp_tend[1] += ν∇²_below*ΔT            # diffusion of temperature

        # full Laplacian in other layers
        ν_above = ν                             # reuse diffusion from above (i.e. k-1/2)
        for k in 2:nlev-1

            ΔT = temp[k+1] - temp[k]            # based on temperature gradient
            ν = get_ν(ν0,ΔT,ΔTmax,Δσ[k])        # get diffusion coefficient at

            # diffusion coefficient ν times 1/Δσ²-like operator
            ν∇²_above = ν_above*∇²_above[k]
            ν∇²_below = ν*∇²_below[k]
            ν∇²_at_k = ν∇²_above + ν∇²_below

            # discrete Laplacian, like the (1, -2, 1)-stencil but for variable Δσ
            u_tend[k] += ν∇²_below*u[k+1] - ν∇²_at_k*u[k] + ν∇²_above*u[k-1]
            v_tend[k] += ν∇²_below*v[k+1] - ν∇²_at_k*v[k] + ν∇²_above*v[k-1]
            temp_tend[k] += ν∇²_below*temp[k+1] - ν∇²_at_k*temp[k] + ν∇²_above*temp[k-1]

            ν_above = ν     # reuse diffusion from k+1/2 for next layer's k-1/2
        end

        # bottom layer with no flux boundary conditions at k=nlev+1/2
        ν∇²_above = ν_above*∇²_above[end]
        u_tend[end] += ν∇²_above*(u[end-1] - u[end])
        v_tend[end] += ν∇²_above*(v[end-1] - v[end])
        temp_tend[end] += ν∇²_above*(temp[end-1] - temp[end])
    end

    return nothing
end

function initialize_vertical_diffusion!(K::ParameterizationConstants,
                                        scheme::VerticalLaplacian,
                                        P::Parameters,
                                        G::Geometry)

    (;vert_diff_∇²_above, vert_diff_∇²_below, vert_diff_Δσ) = K
    Δσ = G.σ_levels_thick
    
    # thickness Δσ of half levels
    @. vert_diff_Δσ = 1/2*(Δσ[2:end] + Δσ[1:end-1])

    # 1/Δσ² but for variable Δσ on half levels
    # = 1/(1/2*Δσₖ(Δσ_k-1 + Δσₖ))
    @. vert_diff_∇²_above = inv(Δσ[2:end]*vert_diff_Δσ)
    
    # = 1/(1/2*Δσₖ(Δσ_k+1 + Δσₖ))
    @. vert_diff_∇²_below = inv(Δσ[1:end-1]*vert_diff_Δσ)

    return nothing
end 