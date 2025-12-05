abstract type AbstractAdiabaticConversion <: AbstractModelComponent end

export AdiabaticConversion

"""Terms used to compute the adiabatic conversion term in the spectral model.
Fields are: $(TYPEDFIELDS)"""
@kwdef struct AdiabaticConversion{VectorType} <: AbstractAdiabaticConversion
    "σ-related factor A needed for adiabatic conversion term"
    σ_lnp_A::VectorType

    "σ-related factor B needed for adiabatic conversion term"
    σ_lnp_B::VectorType
end

Adapt.@adapt_structure AdiabaticConversion

# generator function
AdiabaticConversion(SG::SpectralGrid) = AdiabaticConversion(
    on_architecture(SG.architecture, zeros(SG.NF, SG.nlayers)),
    on_architecture(SG.architecture, zeros(SG.NF, SG.nlayers)))

function initialize!(
    adiabatic::AdiabaticConversion,
    model::PrimitiveEquation,
)
    (; σ_lnp_A, σ_lnp_B) = adiabatic

    # Transfer arrays to CPU for computation
    σ_lnp_A_cpu = on_architecture(CPU(), σ_lnp_A)
    σ_lnp_B_cpu = on_architecture(CPU(), σ_lnp_B)

    # ADIABATIC TERM, Simmons and Burridge, 1981, eq. 3.12
    (; σ_levels_half, σ_levels_full, σ_levels_thick) = model.geometry
    σ_levels_half_cpu = on_architecture(CPU(), σ_levels_half)
    σ_levels_thick_cpu = on_architecture(CPU(), σ_levels_thick)
    
    # precompute ln(σ_k+1/2) - ln(σ_k-1/2) but swap sign, include 1/Δσₖ
    σ_lnp_A_cpu .= log.(σ_levels_half_cpu[1:end-1]./σ_levels_half_cpu[2:end]) ./ σ_levels_thick_cpu
    σ_lnp_A_cpu[1] = 0  # the corresponding sum is 1:k-1 so 0 to replace log(0) from above

    # precompute the αₖ = 1 - p_k-1/2/Δpₖ*log(p_k+1/2/p_k-1/2) term in σ coordinates
    σ_lnp_B_cpu .= 1 .- σ_levels_half_cpu[1:end-1]./σ_levels_thick_cpu .*
                    log.(σ_levels_half_cpu[2:end]./σ_levels_half_cpu[1:end-1])
    σ_lnp_B_cpu[1] = σ_levels_half_cpu[1] <= 0 ? log(2) : σ_lnp_B_cpu[1]    # set α₁ = log(2), eq. 3.19
    σ_lnp_B_cpu .*= -1  # absorb sign from -1/Δσₖ only, eq. 3.12
    
    # Transfer results back to device
    σ_lnp_A .= on_architecture(architecture(σ_lnp_A), σ_lnp_A_cpu)
    σ_lnp_B .= on_architecture(architecture(σ_lnp_B), σ_lnp_B_cpu)
    
    return nothing
end