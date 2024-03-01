abstract type AbstractAdiabaticConversion <: AbstractModelComponent end

export AdiabaticConversion
Base.@kwdef struct AdiabaticConversion{NF} <: AbstractAdiabaticConversion
    nlev::Int

    "σ-related factor A needed for adiabatic conversion term"
    σ_lnp_A::Vector{NF} = zeros(NF,nlev)
    
    "σ-related factor B needed for adiabatic conversion term"
    σ_lnp_B::Vector{NF} = zeros(NF,nlev)
end

AdiabaticConversion(SG::SpectralGrid;kwargs...) = AdiabaticConversion{SG.NF}(;nlev=SG.nlev,kwargs...)

function initialize!(
    adiabatic::AdiabaticConversion,
    model::PrimitiveEquation,
)
    (;σ_lnp_A, σ_lnp_B) = adiabatic

    # ADIABATIC TERM, Simmons and Burridge, 1981, eq. 3.12
    (;σ_levels_half, σ_levels_full, σ_levels_thick) = model.geometry
    # precompute ln(σ_k+1/2) - ln(σ_k-1/2) but swap sign, include 1/Δσₖ
    σ_lnp_A .= log.(σ_levels_half[1:end-1]./σ_levels_half[2:end]) ./ σ_levels_thick
    σ_lnp_A[1] = 0  # the corresponding sum is 1:k-1 so 0 to replace log(0) from above
    
    # precompute the αₖ = 1 - p_k-1/2/Δpₖ*log(p_k+1/2/p_k-1/2) term in σ coordinates
    σ_lnp_B .= 1 .- σ_levels_half[1:end-1]./σ_levels_thick .*
                    log.(σ_levels_half[2:end]./σ_levels_half[1:end-1])
    σ_lnp_B[1] = σ_levels_half[1] <= 0 ? log(2) : σ_lnp_B[1]    # set α₁ = log(2), eq. 3.19
    σ_lnp_B .*= -1  # absorb sign from -1/Δσₖ only, eq. 3.12
    return nothing
end