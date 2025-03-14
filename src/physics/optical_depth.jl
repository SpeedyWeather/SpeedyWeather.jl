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
    # transmittance = exp(-optical_depth), so set to 1
    column.transmittance_shortwave .= 1
    column.transmittance_longwave .= 1
end

export FriersonOpticalDepth
@kwdef mutable struct FriersonOpticalDepth{NF, MatrixType} <: AbstractOpticalDepth
    
    "[OPTION] Spectral band to use"
    band::Int = 1

    nlat::Int
    nlayers::Int

    "[OPTION] Optical depth at the equator"
    τ₀_equator::NF = 6

    "[OPTION] Optical depth at the poles"
    τ₀_pole::NF = 1.5

    "[OPTION] Fraction to mix linear and quadratic profile"
    fₗ::NF = 0.1

    transmittance::MatrixType = zeros(NF, nlayers, nlat)
end

FriersonOpticalDepth(SG::SpectralGrid; kwargs...) = FriersonOpticalDepth{SG.NF, SG.MatrixType}(nlat=SG.nlat, nlayers=SG.nlayers; kwargs...)
function initialize!(od::FriersonOpticalDepth, model::AbstractModel)

    (; τ₀_equator, τ₀_pole, fₗ, transmittance, nlayers) = od
    σ = model.geometry.σ_levels_half
    sinφ = model.geometry.sinlat

    # optical depth τ
    # Frierson 2006, eq. (4), (5) but in a differential form, computing dτ between half levels below and above
    # --- τ(k=1/2)                  # half level above
    # dt = τ(k=1+1/2) - τ(k=1/2)    # differential optical depth on layer k
    # --- τ(k=1+1/2)                # half level below

    NF = eltype(transmittance)
    local τ_above::NF = 0
    for j in eachindex(sinφ)
        τ₀ = τ₀_equator + (τ₀_pole - τ₀_equator)*sinφ[j]^2
        for k in 2:nlayers+1     # loop over half levels below
            τ_below = τ₀*(fₗ*σ[k] + (1 - fₗ)*σ[k]^4)
            dτ = τ_below - τ_above
            transmittance[k-1, j] = exp(-dτ)
            τ_above = τ_below
        end
        τ_above = 0
    end
end

function optical_depth!(
    column::ColumnVariables,
    od::FriersonOpticalDepth,
    model::AbstractModel,
)
    # escape immediately if fewer bands defined in longwave radiation scheme
    od.band > column.nbands_longwave && return nothing

    # Frierson et al. 2006 uses a transparent atmosphere for shortwave radiation
    column.transmittance_shortwave .= 1

    # but the longwave optical depth follows some idealised humidity lat-vert distribution
    # set the transmittance=exp(-optical_deoth) for the longwave band
    column.transmittance_longwave[:, od.band] .= od.transmittance[:, column.jring]
    return nothing
end