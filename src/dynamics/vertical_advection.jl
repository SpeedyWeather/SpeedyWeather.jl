abstract type VerticalAdvection{NF,B} end

# Dispersive and diffusive advection schemes `NF` is the type, `B` the half-stencil size
abstract type DiffusiveVerticalAdvection{NF, B}  <: VerticalAdvection{NF, B} end
abstract type DispersiveVerticalAdvection{NF, B} <: VerticalAdvection{NF, B} end

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

@inline function retrieve_current_stencil(k, layers, var, nlev, ::VerticalAdvection{NF, B}) where {NF, B}
    k_stencil = max.(min.(nlev, k-B:k+B), 1)
    ξ_stencil = Tuple(retrieve_current_time_step(layers[k].grid_variables, var) for k in k_stencil)
    return ξ_stencil
end

@inline function retrieve_previous_stencil(k, layers, var, nlev, ::VerticalAdvection{NF, B}) where {NF, B}
    k_stencil = max.(min.(nlev, k-B:k+B), 1)
    ξ_stencil = Tuple(retrieve_previous_time_step(layers[k].grid_variables, var) for k in k_stencil)
    return ξ_stencil
end

@inline retrieve_stencil(k, layers, var, nlev, scheme::DiffusiveVerticalAdvection)  = retrieve_previous_stencil(k, layers, var, nlev, scheme)
@inline retrieve_stencil(k, layers, var, nlev, scheme::DispersiveVerticalAdvection) =  retrieve_current_stencil(k, layers, var, nlev, scheme)

function vertical_advection!(   layer::DiagnosticVariablesLayer,
                                diagn::DiagnosticVariables,
                                model::PrimitiveEquation)
            
    (; k ) = layer       # which layer are we on?
    
    wet_core = model isa PrimitiveWet
    (; σ_levels_thick, nlev ) = model.geometry
     
    scheme = model.vertical_advection

    # for k==1 "above" term is 0, for k==nlev "below" term is zero
    # avoid out-of-bounds indexing with k_above, k_below as follows
    k⁻ = max(1,k-1) # just saturate, because M_1/2 = 0 (which zeros that term)
        
    # mass fluxes, M_1/2 = M_nlev+1/2 = 0, but k=1/2 isn't explicitly stored
    σ_tend_above = diagn.layers[k⁻].dynamics_variables.σ_tend
    σ_tend_below = layer.dynamics_variables.σ_tend
    
    # layer thickness Δσ on level k
    Δσₖ = σ_levels_thick[k]
    
    for var in (:u, :v, :temp)
        ξ_tend = getproperty(layer.tendencies, Symbol(var, :_tend_grid))
        ξ_sten = retrieve_stencil(k, diagn.layers, var, nlev, scheme)
        ξ      = retrieve_time_step(scheme, layer.grid_variables, var)

        _vertical_advection!(ξ_tend,σ_tend_above,σ_tend_below,ξ_sten,ξ,Δσₖ,scheme)
    end

    if wet_core
        ξ_tend = getproperty(layer.tendencies, :humid_tend_grid)
        ξ_sten = retrieve_current_stencil(k, diagn.layers, :humid, nlev, scheme)
        ξ      = retrieve_current_time_step(layer.grid_variables, :humid)

        _vertical_advection!(ξ_tend,σ_tend_above,σ_tend_below,ξ_sten,ξ,Δσₖ,scheme)
    end
end

# MULTI THREADED VERSION only writes into layer k
function _vertical_advection!(  ξ_tend::Grid,                  # tendency of quantity ξ at k
                                σ_tend_above::Grid,            # vertical velocity at k-1/2
                                σ_tend_below::Grid,            # vertical velocity at k+1/2
                                ξ_sten,                        # ξ stencil for vertical advection (from k-B to k+B)
                                ξ::Grid,                       # ξ at level k
                                Δσₖ::NF,                       # layer thickness on σ levels
                                adv::VerticalAdvection{NF, B}  # vertical advection scheme
                                ) where {NF<:AbstractFloat,Grid<:AbstractGrid{NF}, B}
    Δσₖ⁻¹ = 1/Δσₖ                                      # precompute

    # += as the tendencies already contain the parameterizations
    for ij in eachgridpoint(ξ_tend)
        σ̇⁻ = σ_tend_above[ij]       # velocity into layer k from above
        σ̇⁺ = σ_tend_below[ij]       # velocity out of layer k to below

        ξᶠ⁺ = reconstructed_at_face(ij, adv, σ̇⁺, ξ_sten[2:end])
        ξᶠ⁻ = reconstructed_at_face(ij, adv, σ̇⁻, ξ_sten[1:end-1])

        ξ_tend[ij] -=  Δσₖ⁻¹ * (σ̇⁺ * ξᶠ⁺ - σ̇⁻ * ξᶠ⁻ - ξ[ij] * (σ̇⁺ - σ̇⁻))
    end
end

@inline reconstructed_at_face(ij, ::UpwindVerticalAdvection{NF, 1}, u, ξ) where NF = ifelse(u > 0, ξ[1][ij], ξ[2][ij])
@inline reconstructed_at_face(ij, ::UpwindVerticalAdvection{NF, 2}, u, ξ) where NF = ifelse(u > 0, (2ξ[1][ij] + 5ξ[2][ij] - ξ[3][ij]) / 6,
                                                                                                   (2ξ[4][ij] + 5ξ[3][ij] - ξ[2][ij]) / 6)
@inline reconstructed_at_face(ij, ::UpwindVerticalAdvection{NF, 3}, u, ξ) where NF = ifelse(u > 0, (2ξ[1][ij] - 13ξ[2][ij] + 47ξ[3][ij] + 27ξ[4][ij] - 3ξ[5][ij]) / 60,
                                                                                                   (2ξ[6][ij] - 13ξ[5][ij] + 47ξ[4][ij] + 27ξ[3][ij] - 3ξ[2][ij]) / 60)

@inline reconstructed_at_face(ij, ::CenteredVerticalAdvection{NF, 1}, u, ξ) where NF = ( ξ[1][ij] +  ξ[2][ij]) / 2
@inline reconstructed_at_face(ij, ::CenteredVerticalAdvection{NF, 2}, u, ξ) where NF = (-ξ[1][ij] + 7ξ[2][ij] + 7ξ[3][ij] - ξ[4][ij]) / 12

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

@inline function reconstructed_at_face(ij, ::WENOVerticalAdvection{NF}, u, ξ) where NF
    if u > 0
        S₀ = (ξ[3][ij], ξ[4][ij], ξ[5][ij])
        S₁ = (ξ[2][ij], ξ[3][ij], ξ[4][ij])
        S₂ = (ξ[1][ij], ξ[2][ij], ξ[3][ij])
    else
        S₀ = (ξ[4][ij], ξ[3][ij], ξ[2][ij])
        S₁ = (ξ[5][ij], ξ[4][ij], ξ[3][ij])
        S₂ = (ξ[6][ij], ξ[5][ij], ξ[4][ij])
    end
    return weno_reconstruction(S₀, S₁, S₂, NF)
end
