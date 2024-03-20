abstract type AbstractVerticalDiffusion <: AbstractParameterization end

## FUNCTION BARRIERS for all AbstractVerticalDiffusion
function vertical_diffusion!(   column::ColumnVariables,
                                model::PrimitiveEquation)
    vertical_diffusion!(column, model.vertical_diffusion, model)
end

## NO VERTICAL DIFFUSION
export NoVerticalDiffusion
struct NoVerticalDiffusion <: AbstractVerticalDiffusion end
NoVerticalDiffusion(SG::SpectralGrid) = NoVerticalDiffusion()

# define dummy functions
initialize!(::NoVerticalDiffusion,::PrimitiveEquation) = nothing
vertical_diffusion!(::ColumnVariables, ::NoVerticalDiffusion, ::PrimitiveEquation) = nothing


Base.@kwdef struct BulkRichardsonDiffusion{NF} <: AbstractVerticalDiffusion
    nlev::Int

    "[OPTION] von Kármán constant [1]"
    κ::NF = 0.4

    "[OPTION] roughness length [m]"
    z₀::NF = 3.21e-5

    "[OPTION] Critical Richardson number for stable mixing cutoff [1]"
    Ri_c::NF = 1

    "[OPTION] Surface layer fraction"
    fb::NF = 0.1

    "[OPTION] diffuse static energy?"
    diffuse_static_energy::Bool = true

    "[OPTION] diffuse momentum?"
    diffuse_momentum::Bool = true

    "[OPTION] diffuse humidity? Ignored for PrimitiveDryModels"
    diffuse_humidity::Bool = true

    # precomputed operators
    ∇²_above::Vector{NF} = zeros(NF, nlev-1)
    ∇²_below::Vector{NF} = zeros(NF, nlev-1)
    Δσ::Vector{NF} = zeros(NF, nlev-1)
end

BulkRichardsonDiffusion(SG::SpectralGrid; kwargs...) = BulkRichardsonDiffusion{SG.NF}(; nlev=SG.nlev, kwargs...)

function initialize!(scheme::BulkRichardsonDiffusion, model::PrimitiveEquation)
    Δσ = model.geometry.σ_levels_thick
    
    # 1. thickness Δσ of half levels
    # 2. 1/Δσ² but for variable Δσ on half levels = 1/(1/2*Δσₖ(Δσ_k-1 + Δσₖ))
    # 3. as 2. but for half levels below = 1/(1/2*Δσₖ(Δσ_k+1 + Δσₖ))
    @. scheme.Δσ = 1/2*(Δσ[2:end] + Δσ[1:end-1])
    @. scheme.∇²_above = inv(Δσ[2:end]*scheme.Δσ)
    @. scheme.∇²_below = inv(Δσ[1:end-1]*scheme.Δσ)

    return nothing
end

function vertical_diffusion!(
    column::ColumnVariables,
    scheme::BulkRichardsonDiffusion,
    model::PrimitiveEquation)

    (; diffuse_momentum, diffuse_static_energy, diffuse_humidity) = scheme
    any(( diffuse_momentum, diffuse_static_energy, diffuse_humidity)) || return nothing
    
    K = get_diffusion_coefficients!(column, scheme)

    (; u, u_tend, v, v_tend, dry_static_energy, temp_tend, humid, humid_tend) = column
    diffuse_momentum                            && vertical_diffusion!(u_tend, u, K, scheme)
    diffuse_momentum                            && vertical_diffusion!(v_tend, v, K, scheme)
    diffuse_static_energy                       && vertical_diffusion!(temp_tend, dry_static_energy, K, scheme)
    model isa PrimitiveWet && diffuse_humidity  && vertical_diffusion!(humid_tend, humid, K, scheme)
    return nothing
end

function get_diffusion_coefficients!(column, scheme)
    Ri = column.bulk_richardson
    K = column.a        # reuse work array for diffusion coefficients
    (; Ri_c) = scheme
    (; nlev) = column

    # calculate here
    column.boundary_layer_depth = nlev+1
    fill!(K,0)

    return K
end

function vertical_diffusion!(   
    tend::AbstractVector,
    var::AbstractVector,
    K::AbstractVector,
    scheme::BulkRichardsonDiffusion,
)
    (; Δσ, ∇²_above, ∇²_below) = scheme
    nlev = length(tend)

    # @inbounds begin
        # top layer with no flux boundary conditions at k=1/2
        tend[1] += K[1] * ∇²_below[1] * (var[2] - var[1])

        # full Laplacian in other layers
        for k in 2:nlev-1
            # discrete Laplacian, like the (1, -2, 1)-stencil but for variable Δσ
            ∇²_at_k = ∇²_above[k-1] + ∇²_below[k]
            tend[k] += K[k]*(∇²_below[k]*var[k+1] - ∇²_at_k*var[k] + ∇²_above[k-1]*var[k-1])
        end

        # bottom layer with no flux boundary conditions at k=nlev+1/2
        tend[end] += K[end] * ∇²_above[end] * (var[end-1] - var[end])
    # end
end