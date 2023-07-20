# Dispersive and diffusive advection schemes `NF` is the type, `B` the half-stencil size
abstract type DiffusiveVerticalAdvection{NF, B}  <: VerticalAdvection{NF} end
abstract type DispersiveVerticalAdvection{NF, B} <: VerticalAdvection{NF} end

struct FirstOrderUpwind{NF}    <: DiffusiveVerticalAdvection{NF, 1} end
struct ThirdOrderUpwind{NF}    <: DiffusiveVerticalAdvection{NF, 2} end
struct WENO3{NF}               <: DiffusiveVerticalAdvection{NF, 2} end

struct SecondOrderCentered{NF} <: DispersiveVerticalAdvection{NF, 1} end
struct FourthOrderCentered{NF} <: DispersiveVerticalAdvection{NF, 2} end

for T in (:FirstOrderUpwind, :ThirdOrderUpwind, :WENO3, :SecondOrderCentered, :FourthOrderCentered)
    @eval $T(spectral_grid::SpectralGrid) = $T{spectral_grid.NF}()
end

@inline retrieve_time_step(::DiffusiveVerticalAdvection,  variables, var) = getproperty(variables, Symbol(var, :_grid_prev))
@inline retrieve_time_step(::DispersiveVerticalAdvection, variables, var) = getproperty(variables, Symbol(var, :_grid_prev))

function vertical_advection!(   layer::DiagnosticVariablesLayer,
                                diagn::DiagnosticVariables,
                                model::PrimitiveEquation)
            
    (; k ) = layer       # which layer are we on?
    
    wet_core = model isa PrimitiveWet
    (; σ_levels_thick, nlev ) = model.geometry
     
    vertical_advection = model.vertical_advection

    # for k==1 "above" term is 0, for k==nlev "below" term is zero
    # avoid out-of-bounds indexing with k_above, k_below as follows
    k_above = max(1,k-1)        # just saturate, because M_1/2 = 0 (which zeros that term)
    k_below = min(k+1,nlev)     # just saturate, because M_nlev+1/2 = 0 (which zeros that term)
        
    # mass fluxes, M_1/2 = M_nlev+1/2 = 0, but k=1/2 isn't explicitly stored
    σ_tend_above = diagn.layers[k_above].dynamics_variables.σ_tend
    σ_tend_below = layer.dynamics_variables.σ_tend
    
    # layer thickness Δσ on level k
    Δσₖ = σ_levels_thick[k]
    
    for var in (:u, :v, :temp)
        ξ_tend  = getproperty(layer.tendencies, Symbol(var, :_tend_grid))
        ξ_above = retrieve_time_step(vertical_advection, diagn.layers[k_above].grid_variables, var)
        ξ       = retrieve_time_step(vertical_advection, layer.grid_variables, var)
        ξ_below = retrieve_time_step(vertical_advection, diagn.layers[k_below].grid_variables, var)

        _vertical_advection!(ξ_tend,σ_tend_above,σ_tend_below,ξ_above,ξ,ξ_below,Δσₖ,vertical_advection)
    end

    if wet_core
        ξ_tend  = getproperty(layer.tendencies, :humid_tend_grid)
        ξ_above = retrieve_time_step(vertical_advection, diagn.layers[k_above].grid_variables, :humid)
        ξ       = retrieve_time_step(vertical_advection, layer.grid_variables, :humid)
        ξ_below = retrieve_time_step(vertical_advection, diagn.layers[k_below].grid_variables, :humid)
    
        _vertical_advection!(q_tend,σ_tend_above,σ_tend_below,q_above,q,q_below,Δσₖ,vertical_advection)
    end
end

@inline reconstructed_at_face(::FirstOrderUpwind,    u, c⁻, c⁺) = ifelse(u > 0, c⁻, c⁺)
@inline reconstructed_at_face(::SecondOrderCentered, u, c⁻, c⁺) = (c⁻ + c⁺) / 2

# MULTI THREADED VERSION only writes into layer k
function _vertical_advection!(  ξ_tend::Grid,           # tendency of quantity ξ at k
                                σ_tend_above::Grid,     # vertical velocity at k-1/2
                                σ_tend_below::Grid,     # vertical velocity at k+1/2
                                ξ_above::Grid,          # quantity ξ at k-1
                                ξ::Grid,                # quantity ξ at k
                                ξ_below::Grid,          # quantity ξ at k+1
                                Δσₖ::NF,                # layer thickness on σ levels
                                vertical_advection      # vertical advection scheme
                                ) where {NF<:AbstractFloat,Grid<:AbstractGrid{NF}}
    Δσₖ⁻¹ = 1/Δσₖ                                      # precompute

    # += as the tendencies already contain the parameterizations
    for ij in eachgridpoint(ξ_tend,σ_tend_above,σ_tend_below,ξ_above,ξ,ξ_below)
        σ̇⁻ = σ_tend_above[ij]       # velocity into layer k from above
        σ̇⁺ = σ_tend_below[ij]       # velocity out of layer k to below

        ξₖ₋₁ = ξ_above[ij]
        ξₖ   = ξ[ij]
        ξₖ₊₁ = ξ_below[ij]

        ξ⁺ = reconstructed_at_face(vertical_advection, σ̇⁺, ξₖ, ξₖ₊₁)
        ξ⁻ = reconstructed_at_face(vertical_advection, σ̇⁻, ξₖ₋₁, ξₖ)

        ξ_tend[ij] -=  Δσₖ⁻¹ * (σ̇⁺ * ξ⁺ - σ̇⁻ * ξ⁻ - ξₖ * (σ̇⁺ - σ̇⁻))
    end
end
