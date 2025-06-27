abstract type AbstractVerticalAdvection end
abstract type VerticalAdvection{NF, B} <: AbstractVerticalAdvection end

# Dispersive and diffusive advection schemes `NF` is the type, `B` the half-stencil size
abstract type DiffusiveVerticalAdvection{NF, B}  <: VerticalAdvection{NF, B} end
abstract type DispersiveVerticalAdvection{NF, B} <: VerticalAdvection{NF, B} end

export UpwindVerticalAdvection, WENOVerticalAdvection, CenteredVerticalAdvection
struct UpwindVerticalAdvection{NF, B}   <: DiffusiveVerticalAdvection{NF, B} end
struct WENOVerticalAdvection{NF}        <: DiffusiveVerticalAdvection{NF, 3} end
struct CenteredVerticalAdvection{NF, B} <: DispersiveVerticalAdvection{NF, B} end

CenteredVerticalAdvection(spectral_grid; order = 2) = CenteredVerticalAdvection{spectral_grid.NF, order       ÷ 2}()
UpwindVerticalAdvection(spectral_grid; order = 5)   =   UpwindVerticalAdvection{spectral_grid.NF, (order + 1) ÷ 2}()
WENOVerticalAdvection(spectral_grid)                =     WENOVerticalAdvection{spectral_grid.NF}()

@inline retrieve_previous_time_step(variables, var) = getproperty(variables, Symbol(var, :_grid_prev))
@inline  retrieve_current_time_step(variables, var) = getproperty(variables, Symbol(var, :_grid))

@inline retrieve_time_step(::DiffusiveVerticalAdvection,  variables, var) = retrieve_previous_time_step(variables, var)
@inline retrieve_time_step(::DispersiveVerticalAdvection, variables, var) =  retrieve_current_time_step(variables, var)

@inline function retrieve_stencil(k, nlayers, ::VerticalAdvection{NF, B}) where {NF, B}
    # creates allocation-free tuples for k-B:k+B but clamped into (1, nlayers)
    # e.g. (1, 1, 2), (1, 2, 3), (2, 3, 4) ... (for k=1, 2, 3; B=1)
    return ntuple(i -> clamp(i+k-B-1, 1, nlayers), 2B+1)
end

function vertical_advection!(
    diagn::DiagnosticVariables,
    model::PrimitiveEquation,
)
    Δσ = model.geometry.σ_levels_thick
    advection_scheme = model.vertical_advection
    (; σ_tend) = diagn.dynamics
    
    for var in (:u, :v, :temp)
        ξ_tend = getproperty(diagn.tendencies, Symbol(var, :_tend_grid))
        ξ      = retrieve_time_step(advection_scheme, diagn.grid, var)
        _vertical_advection!(ξ_tend, σ_tend, ξ, Δσ, advection_scheme)
    end

    if model isa PrimitiveWet   # advect humidity only with primitive wet core
        ξ_tend = getproperty(diagn.tendencies, :humid_tend_grid)
        ξ      = retrieve_time_step(advection_scheme, diagn.grid, :humid)
        _vertical_advection!(ξ_tend, σ_tend, ξ, Δσ, advection_scheme)
    end

    for tracer in values(model.tracers)
        if tracer.active
            ξ_tend = diagn.tendencies.tracers_tend_grid[tracer.name]
            ξ      = retrieve_time_step(advection_scheme, diagn.grid, :tracers)[tracer.name]
            _vertical_advection!(ξ_tend, σ_tend, ξ, Δσ, advection_scheme)
        end
    end
end

function _vertical_advection!(
    ξ_tend::AbstractField,      # tendency of quantity ξ
    σ_tend::AbstractField,      # vertical velocity at k+1/2
    ξ::AbstractField,           # ξ
    Δσ,                         # layer thickness on σ levels
    adv::VerticalAdvection      # vertical advection scheme of order B
)
    grids_match(ξ_tend, σ_tend, ξ) || throw(DimensionMismatch(ξ_tend, σ_tend, ξ))

    nlayers = size(ξ, 2)
    @inbounds for k in 1:nlayers
        Δσₖ⁻¹ = inv(Δσ[k])          # inverse layer thickness, compute inv only once

        # for k=1 "above" term (at k-1/2) is 0, for k==nlayers "below" term (at k+1/2) is zero
        # avoid out-of-bounds indexing with k⁻, k⁺
        k⁻ = max(1, k-1)    # TODO check that this actually zeros velocity at k=1/2
        k⁺ = k

        k_stencil = retrieve_stencil(k, nlayers, adv)

        for ij in eachgridpoint(ξ_tend)
            σ̇⁻ = σ_tend[ij, k⁻]       # velocity into layer k from above
            σ̇⁺ = σ_tend[ij, k⁺]       # velocity out of layer k to below
            
            ξᶠ⁺ = reconstructed_at_face(ξ, ij, k_stencil[2:end],   σ̇⁺, adv)
            ξᶠ⁻ = reconstructed_at_face(ξ, ij, k_stencil[1:end-1], σ̇⁻, adv)
            
            # -= as the tendencies already contain the parameterizations
            ξ_tend[ij, k] -=  Δσₖ⁻¹ * (σ̇⁺ * ξᶠ⁺ - σ̇⁻ * ξᶠ⁻ - ξ[ij, k] * (σ̇⁺ - σ̇⁻))
        end
    end
end

# 1st order upwind
@inline reconstructed_at_face(ξ, ij, k, u, ::UpwindVerticalAdvection{NF, 1}) where NF =
    ifelse(u > 0,   ξ[ij, k[1]],
                    ξ[ij, k[2]])

# 3rd order upwind
@inline reconstructed_at_face(ξ, ij, k, u, ::UpwindVerticalAdvection{NF, 2}) where NF =
    ifelse(u > 0,   (2ξ[ij, k[1]] + 5ξ[ij, k[2]] - ξ[ij, k[3]]) / 6,
                    (2ξ[ij, k[4]] + 5ξ[ij, k[3]] - ξ[ij, k[2]]) / 6)

# 5th order upwind
@inline reconstructed_at_face(ξ, ij, k, u, ::UpwindVerticalAdvection{NF, 3}) where NF =
    ifelse(u > 0,   (2ξ[ij, k[1]] - 13ξ[ij, k[2]] + 47ξ[ij, k[3]] + 27ξ[ij, k[4]] - 3ξ[ij, k[5]]) / 60,
                    (2ξ[ij, k[6]] - 13ξ[ij, k[5]] + 47ξ[ij, k[4]] + 27ξ[ij, k[3]] - 3ξ[ij, k[2]]) / 60)

# 2nd order centered
@inline reconstructed_at_face(ξ, ij, k, u, ::CenteredVerticalAdvection{NF, 1}) where NF =
    (ξ[ij, k[1]] +  ξ[ij, k[2]]) / 2

# 4th order centered
@inline reconstructed_at_face(ξ, ij, k, u, ::CenteredVerticalAdvection{NF, 2}) where NF =
    (-ξ[ij, k[1]] + 7ξ[ij, k[2]] + 7ξ[ij, k[3]] - ξ[ij, k[4]]) / 12

const ε  = 1e-6
const d₀ = 3/10
const d₁ = 3/5
const d₂ = 1/10

@inline weight_β₀(S, NF) = NF(13/12) * (S[1] - 2S[2] + S[3])^2 + NF(1/4) * (3S[1] - 4S[2] +  S[3])^2
@inline weight_β₁(S, NF) = NF(13/12) * (S[1] - 2S[2] + S[3])^2 + NF(1/4) * ( S[1]         -  S[3])^2
@inline weight_β₂(S, NF) = NF(13/12) * (S[1] - 2S[2] + S[3])^2 + NF(1/4) * ( S[1] - 4S[2] + 3S[3])^2

@inline p₀(S) = (2S[1] + 5S[2] -   S[3]) / 6 # downind stencil
@inline p₁(S) = (-S[1] + 5S[2] +  2S[3]) / 6 # upwind stencil
@inline p₂(S) = (2S[1] - 7S[2] + 11S[3]) / 6 # extrapolating stencil

@inline τ₅(β₀, β₁, β₂) = abs(β₂ - β₀)

@inline function weno_reconstruction(S₀, S₁, S₂, NF)
    β₀ = weight_β₀(S₀, NF)
    β₁ = weight_β₁(S₁, NF)
    β₂ = weight_β₂(S₂, NF)

    w₀ = NF(d₀) * (1 + (τ₅(β₀, β₁, β₂) / (β₀ + NF(ε)))^2)
    w₁ = NF(d₁) * (1 + (τ₅(β₀, β₁, β₂) / (β₁ + NF(ε)))^2)
    w₂ = NF(d₂) * (1 + (τ₅(β₀, β₁, β₂) / (β₂ + NF(ε)))^2)

    w₀, w₁, w₂ = (w₀, w₁, w₂) ./ (w₀ + w₁ + w₂) 

    return p₀(S₀) * w₀ + p₁(S₁) * w₁ + p₂(S₂) * w₂
end

@inline function reconstructed_at_face(ξ, ij, k, u, ::WENOVerticalAdvection{NF}) where NF
    if u > 0
        S₀ = (ξ[ij, k[3]], ξ[ij, k[4]], ξ[ij, k[5]])
        S₁ = (ξ[ij, k[2]], ξ[ij, k[3]], ξ[ij, k[4]])
        S₂ = (ξ[ij, k[1]], ξ[ij, k[2]], ξ[ij, k[3]])
    else
        S₀ = (ξ[ij, k[4]], ξ[ij, k[3]], ξ[ij, k[2]])
        S₁ = (ξ[ij, k[5]], ξ[ij, k[4]], ξ[ij, k[3]])
        S₂ = (ξ[ij, k[6]], ξ[ij, k[5]], ξ[ij, k[4]])
    end
    return weno_reconstruction(S₀, S₁, S₂, NF)
end
