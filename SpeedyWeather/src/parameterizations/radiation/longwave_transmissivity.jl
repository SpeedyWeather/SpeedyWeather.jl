abstract type AbstractLongwaveTransmissivity <: AbstractLongwave end

export TransparentLongwaveTransmissivity
TransparentLongwaveTransmissivity(SG::SpectralGrid) = ConstantLongwaveTransmissivity(SG, transmissivity = 1)

export ConstantLongwaveTransmissivity
@parameterized @kwdef struct ConstantLongwaveTransmissivity{NF} <: AbstractLongwaveTransmissivity
    @param transmissivity::NF = 0.6 (bounds = 0 .. 1,)
end
Adapt.@adapt_structure ConstantLongwaveTransmissivity
ConstantLongwaveTransmissivity(SG::SpectralGrid; kwargs...) = ConstantLongwaveTransmissivity{SG.NF}(; kwargs...)
initialize!(::ConstantLongwaveTransmissivity, ::AbstractModel) = nothing
@propagate_inbounds function transmissivity!(ij, vars, CLT::ConstantLongwaveTransmissivity, model)
    t = vars.scratch.grid.a
    nlayers = size(t, 2)

    τ = -log(CLT.transmissivity)            # total optical depth of the atmosphere
    dσ = model.geometry.σ_levels_thick      # divide optical depth wrt to pressure thickness of each layer
    for k in 1:nlayers
        t[ij, k] = exp(-τ * dσ[k])          # transmissivity through layer k
    end
    return t
end

export FriersonLongwaveTransmissivity
@parameterized @kwdef struct FriersonLongwaveTransmissivity{NF} <: AbstractLongwaveTransmissivity
    "[OPTION] Optical depth at the equator"
    @param τ₀_equator::NF = 6 (bounds = Nonnegative,)

    "[OPTION] Optical depth at the poles"
    @param τ₀_pole::NF = 1.5 (bounds = Nonnegative,)

    "[OPTION] Fraction to mix linear and quadratic profile"
    @param fₗ::NF = 0.1 (bounds = 0 .. 1,)
end

Adapt.@adapt_structure FriersonLongwaveTransmissivity
FriersonLongwaveTransmissivity(SG::SpectralGrid; kwargs...) = FriersonLongwaveTransmissivity{SG.NF}(; kwargs...)
initialize!(::FriersonLongwaveTransmissivity, ::AbstractModel) = nothing

@propagate_inbounds function transmissivity!(ij, vars, transmissivity::FriersonLongwaveTransmissivity, model)

    # use scratch array to compute transmissivity t
    t = vars.scratch.grid.a
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
    τ₀ = τ₀_equator + (τ₀_pole - τ₀_equator) * sind(θ)^2
    for k in 2:(nlayers + 1)        # loop over half levels below
        τ_below = τ₀ * (fₗ * σ[k] + (1 - fₗ) * σ[k]^4)
        t[ij, k - 1] = exp(-(τ_below - τ_above))
        τ_above = τ_below
    end

    # return so the radiative_trasfer uses the right scratch array
    return t
end
