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
    on_architecture(SG.architecture, zeros(SG.NF, SG.nlayers))
)

function initialize!(
        adiabatic::AdiabaticConversion,
        model::PrimitiveEquation,
    )
    (; σ_lnp_A, σ_lnp_B) = adiabatic

    # ADIABATIC TERM, Simmons and Burridge, 1981, eq. 3.12
    (; σ_levels_half, σ_levels_thick) = model.geometry
    nlayers = length(σ_levels_thick)

    σ_half_lower = @view σ_levels_half[1:nlayers]       # σ_k-1/2
    σ_half_upper = @view σ_levels_half[2:(nlayers + 1)]  # σ_k+1/2
    ks = 1:nlayers

    # clamp σ_half_lower to avoid log(0) at k=1 where σ_levels_half[1] = 0
    NF = eltype(σ_lnp_A)
    σ_half_lower_safe = max.(σ_half_lower, eps(NF))

    # precompute ln(σ_k+1/2) - ln(σ_k-1/2) but swap sign, include 1/Δσₖ
    # set k=1 to 0 (the corresponding sum is 1:k-1 so 0 to replace log(0))
    σ_lnp_A .= log.(σ_half_lower_safe ./ σ_half_upper) ./ σ_levels_thick .* (ks .> 1)

    # precompute the αₖ = 1 - p_k-1/2/Δpₖ*log(p_k+1/2/p_k-1/2) term in σ coordinates
    # set α₁ = log(2) when σ_levels_half[1] <= 0, eq. 3.19
    # absorb sign from -1/Δσₖ only, eq. 3.12
    σ_lnp_B .= .-(1 .- σ_half_lower_safe ./ σ_levels_thick .* log.(σ_half_upper ./ σ_half_lower_safe)) .* (ks .> 1) .-
        log(NF(2)) .* (ks .== 1)

    return nothing
end
