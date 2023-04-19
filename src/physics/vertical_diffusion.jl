struct NoVerticalDiffusion{NF} <: VerticalDiffusion{NF} end

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
    time_scale::NF = 1  # [hours] time scale for slow global relaxation
    umax::NF = 80       # velocity scale [m/s] to scale the diffusion with |u|/umax
    ΔTmax::NF = 50      # temperature gradient scale [K] (divided by Δσ)

    scale_with_temperature_gradient::Bool = true
    scale_with_speed::Bool = false
    pressure_coordinates::Bool = true
end

# generator so that arguments are converted to Float64
VerticalLaplacian(;kwargs...) = VerticalLaplacian{Float64}(;kwargs...)

function vertical_diffusion!(   column::ColumnVariables,
                                scheme::VerticalLaplacian,
                                model::PrimitiveEquation)

    (;nlev,u_tend,v_tend,temp_tend) = column
    (;u,v,temp) = column
    ∇²_below = model.parameterization_constants.vert_diff_∇²_below
    ∇²_above = model.parameterization_constants.vert_diff_∇²_above
    ν = model.parameterization_constants.vert_diff_ν
    ∂σ = model.parameterization_constants.vert_diff_∂σ

    # GET DIFFUSION COEFFICIENT as a function of u,v,temp and surface pressure
    ν .= model.geometry.radius*inv(scheme.time_scale*3600)    # [hrs] → [s]
    scheme.pressure_coordinates && ν .*= (model.parameters.pres_ref*100/column.pres[end])^2
    if scheme.scale_with_temperature_gradient
        @inbounds for k in 1:nlev-1
            ΔT = (temp[k+1] - temp[k])/∂σ[k]
            ν[k] *= max(ΔT/scheme.ΔTmax + 1,0)/2
        end
    end

    # DO DIFFUSION
    # top layer with no flux boundary conditions at k=1/2
    @inbounds begin
    ν∇²_below = ν[1]*∇²_below[1]
    u_tend[1] += ν∇²_below*(u[2] - u[1])
    v_tend[1] += ν∇²_below*(v[2] - v[1])
    temp_tend[1] += ν∇²_below*(temp[2] - temp[1])

    # full Laplacian in other layers
    for k in 2:nlev-1
        ν∇²_above = ν[k-1]*∇²_above[k]
        ν∇²_below = ν[k]*∇²_below[k]
        ν∇²_at_k = ν∇²_above + ν∇²_below
        u_tend[k] += ν∇²_below*u[k+1] - ν∇²_at_k*u[k] + ν∇²_above*u[k-1]
        v_tend[k] += ν∇²_below*v[k+1] - ν∇²_at_k*v[k] + ν∇²_above*v[k-1]
        temp_tend[k] += ν∇²_below*temp[k+1] - ν∇²_at_k*temp[k] + ν∇²_above*temp[k-1]
    end

    # bottom layer with no flux boundary conditions at k=nlev+1/2
    ν∇²_above = ν[end]*∇²_above[end]
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
    (;vert_diff_∇²_above, vert_diff_∇²_below, vert_diff_∂σ) = K
    Δσ = G.σ_levels_thick
    
    @. vert_diff_∂σ = 2(Δσ[2:end] + Δσ[1:end-1])
    @. vert_diff_∇²_above = Δσ[2:end]*vert_diff_∂σ
    @. vert_diff_∇²_below = Δσ[1:end-1]*vert_diff_∂σ

    return nothing
end 