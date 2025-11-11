abstract type AbstractLongwaveTransmittance <: AbstractLongwave end

export FriersonTransmittance
@kwdef mutable struct FriersonTransmittance{NF} <: AbstractLongwaveTransmittance
    "[OPTION] Spectral band to use"
    band::Int = 1

    "[OPTION] Optical depth at the equator"
    τ₀_equator::NF = 6

    "[OPTION] Optical depth at the poles"
    τ₀_pole::NF = 1.5

    "[OPTION] Fraction to mix linear and quadratic profile"
    fₗ::NF = 0.1
end

FriersonTransmittance(SG::SpectralGrid; kwargs...) = FriersonTransmittance{SG.NF}(; kwargs...)
initialize!(od::FriersonTransmittance, model::AbstractModel) = nothing

function transmittance!(
    column::ColumnVariables{NF},
    od::FriersonTransmittance,
    model::AbstractModel,
) where NF

    # escape immediately if fewer bands defined in longwave radiation scheme
    od.band > column.nbands_longwave && return nothing

    # Frierson et al. 2006 uses a transparent atmosphere for shortwave radiation
    column.transmittance_shortwave .= 1

    # but the longwave optical depth follows some idealised humidity lat-vert distribution
    transmittance = column.transmittance_longwave
    (; τ₀_equator, τ₀_pole, fₗ, band) = od

    # coordinates 
    σ = model.geometry.σ_levels_half
    θ = column.latd
    (; nlayers) = column

    # Frierson 2006, eq. (4), (5) but in a differential form, computing dτ between half levels below and above
    # --- τ(k=1/2)                  # half level above
    # dt = τ(k=1+1/2) - τ(k=1/2)    # differential optical depth on layer k
    # --- τ(k=1+1/2)                # half level below

    local τ_above::NF = 0
    τ₀ = τ₀_equator + (τ₀_pole - τ₀_equator)*sind(θ)^2
    for k in 2:nlayers+1     # loop over half levels below
        τ_below = τ₀*(fₗ*σ[k] + (1 - fₗ)*σ[k]^4)
        transmittance[k-1, band] = exp(-(τ_below - τ_above))
        τ_above = τ_below
    end
end