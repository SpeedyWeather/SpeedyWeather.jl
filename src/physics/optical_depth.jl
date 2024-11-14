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
    column.optical_depth .= 0
end

export FriersonOpticalDepth
@kwdef mutable struct FriersonOpticalDepth{NF} <: AbstractOpticalDepth
    τ₀_equator::NF = 6
    τ₀_pole::NF = 1.5
    fₗ::NF = 0.1
end

FriersonOpticalDepth(SG::SpectralGrid; kwargs...) = FriersonOpticalDepth{SG.NF}(; kwargs...)
initialize!(od::FriersonOpticalDepth, model::AbstractModel) = nothing

function optical_depth!(column::ColumnVariables, od::FriersonOpticalDepth, model::AbstractModel)

    (; optical_depth) = column
    (; τ₀_equator, τ₀_pole, fₗ) = od

    # coordinates 
    σ = model.geometry.σ_levels_half
    θ = column.latd

    # Frierson 2006, eq. (4), (5)
    τ₀ = τ₀_equator + (τ₀_pole - τ₀_equator)*sind(θ)^2
    @. optical_depth = τ₀*(fₗ*σ + (1 - fₗ)*σ^4)
    return nothing
end