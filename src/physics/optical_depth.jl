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
@parameterized @kwdef mutable struct FriersonOpticalDepth{NF} <: AbstractOpticalDepth
    
    "[OPTION] Spectral band to use"
    band::Int = 1

    "[OPTION] Optical depth at the equator"
    @param τ₀_equator::NF = 6 (bounds=Positive,)

    "[OPTION] Optical depth at the poles"
    @param τ₀_pole::NF = 1.5 (bounds=Positive,)

    "[OPTION] Fraction to mix linear and quadratic profile"
    @param fₗ::NF = 0.1 (bounds=UnitInterval,)
end

FriersonOpticalDepth(SG::SpectralGrid; kwargs...) = FriersonOpticalDepth{SG.NF}(; kwargs...)
initialize!(od::FriersonOpticalDepth, model::AbstractModel) = nothing

function optical_depth!(
    column::ColumnVariables{NF},
    od::FriersonOpticalDepth,
    model::AbstractModel,
) where NF

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
    (; nlayers) = column

    # Frierson 2006, eq. (4), (5) but in a differential form, computing dτ between half levels below and above
    # --- τ(k=1/2)                  # half level above
    # dt = τ(k=1+1/2) - τ(k=1/2)    # differential optical depth on layer k
    # --- τ(k=1+1/2)                # half level below

    local τ_above::NF = 0
    τ₀ = τ₀_equator + (τ₀_pole - τ₀_equator)*sind(θ)^2
    for k in 2:nlayers+1     # loop over half levels below
        τ_below = τ₀*(fₗ*σ[k] + (1 - fₗ)*σ[k]^4)
        optical_depth[k-1, band] = τ_below - τ_above
        τ_above = τ_below
    end
end