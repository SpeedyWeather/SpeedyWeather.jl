abstract type AbstractLongwaveTransmissivity <: AbstractLongwave end

export TransparentLongwaveTransmissivity
struct TransparentLongwaveTransmissivity <: AbstractLongwaveTransmissivity end
TransparentLongwaveTransmissivity(SG::SpectralGrid) = TransparentLongwaveTransmissivity()
initialize!(::TransparentLongwaveTransmissivity, ::AbstractModel) = nothing
@propagate_inbounds function transmissivity!(ij, diagn, progn, transmissivity::TransparentLongwaveTransmissivity, model)
    # transmissivity is 1 everywhere (no absorption)
    t = diagn.dynamics.a_grid   # use scratch array
    nlayers = size(t, 2)
    for k in 1:nlayers
        t[ij, k] = one(eltype(t))
    end
    return t
end

export FriersonLongwaveTransmissivity
@kwdef mutable struct FriersonLongwaveTransmissivity{NF} <: AbstractLongwaveTransmissivity
    "[OPTION] Optical depth at the equator"
    τ₀_equator::NF = 6

    "[OPTION] Optical depth at the poles"
    τ₀_pole::NF = 1.5

    "[OPTION] Fraction to mix linear and quadratic profile"
    fₗ::NF = 0.1
end

FriersonLongwaveTransmissivity(SG::SpectralGrid; kwargs...) = FriersonLongwaveTransmissivity{SG.NF}(; kwargs...)
initialize!(::FriersonLongwaveTransmissivity, ::AbstractModel) = nothing

@propagate_inbounds function transmissivity!(ij, diagn, progn, transmissivity::FriersonLongwaveTransmissivity, model)

    # use scratch array to compute transmissivity t
    t = diagn.dynamics.a_grid
    nlayers = size(t, 2)
    NF = eltype(t)

    # but the longwave optical depth follows some idealised humidity lat-vert distribution
    (; τ₀_equator, τ₀_pole, fₗ) = transmissivity

    # coordinates 
    σ = model.geometry.σ_levels_half
    θ = model.geometry.latds[ij]

    # Frierson 2006, eq. (4), (5) but in a differential form, computing dτ between half levels below and above
    # --- τ(k=1/2)                  # half level above
    # dt = τ(k=1+1/2) - τ(k=1/2)    # differential optical depth on layer k
    # --- τ(k=1+1/2)                # half level below

    τ_above::NF = 0
    τ₀ = τ₀_equator + (τ₀_pole - τ₀_equator)*sind(θ)^2
    for k in 2:nlayers+1        # loop over half levels below
        τ_below = τ₀*(fₗ*σ[k] + (1 - fₗ)*σ[k]^4)
        t[ij, k-1] = exp(-(τ_below - τ_above))
        τ_above = τ_below
    end

    # return so the radiative_trasfer uses the right scratch array
    return t
end