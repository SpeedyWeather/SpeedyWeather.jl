abstract type AbstractOpticalDepth <: AbstractParameterization end

# function barrier to dispatch to type of model.optical_depth
function optical_depth!(column::ColumnVariables, model::AbstractModel)
    optical_depth!(column, model.optical_depth, model)
end

export ZeroOpticalDepth
struct ZeroOpticalDepth{NF} <: AbstractOpticalDepth end
ZeroOpticalDepth(SG::SpectralGrid) = ZeroOpticalDepth{SG.NF}()
initialize!(od::ZeroOpticalDepth, ::AbstractModel) = nothing
function optical_depth!(column::ColumnVariables, od::ZeroOpticalDepth, model::AbstractModel)
    column.optical_depth_longwave .= 0
    column.optical_depth_shortwave .= 0
end

export FriersonOpticalDepth
@kwdef mutable struct FriersonOpticalDepth{NF} <: AbstractOpticalDepth
    
    "[OPTION] Spectral band to use"
    band::Int = 1

    "[OPTION] Optical depth at the equator"
    τ₀_equator::NF = 6

    "[OPTION] Optical depth at the poles"
    τ₀_pole::NF = 1.5

    "[OPTION] Fraction to mix linear and quadratic profile"
    fₗ::NF = 0.1
end

FriersonOpticalDepth(SG::SpectralGrid; kwargs...) = FriersonOpticalDepth{SG.NF}(; kwargs...)
initialize!(od::FriersonOpticalDepth, model::AbstractModel) = nothing

function optical_depth!(column::ColumnVariables, od::FriersonOpticalDepth, model::AbstractModel)

    # escape immediately if fewer bands defined in longwave radiation scheme
    od.band > column.nbands_longwave && return nothing

    # Frierson et al. 2006 uses a transparent atmosphere for shortwave radiation
    column.optical_depth_shortwave .= 0

    # but the longwave optical depth follows some idealised humidity lat-vert distribution
    optical_depth = column.optical_depth_longwave
    (; τ₀_equator, τ₀_pole, fₗ, band) = od

    # coordinates 
    σ = model.geometry.σ_levels_half
    θ = column.latd

    # Frierson 2006, eq. (4), (5)
    τ₀ = τ₀_equator + (τ₀_pole - τ₀_equator)*sind(θ)^2
    @inbounds for k in eachindex(σ)     # loop over half levels
        optical_depth[k, band] = τ₀*(fₗ*σ[k] + (1 - fₗ)*σ[k]^4)
    end
end