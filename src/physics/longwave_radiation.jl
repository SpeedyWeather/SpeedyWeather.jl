abstract type AbstractLongwave <: AbstractRadiation end

export NoLongwave
struct NoLongwave <: AbstractLongwave end
NoLongwave(SG::SpectralGrid) = NoLongwave()
initialize!(::NoLongwave, ::PrimitiveEquation) = nothing
longwave_radiation!(::ColumnVariables, ::NoLongwave, ::PrimitiveEquation) = nothing

# function barrier for all AbstractLongwave
function longwave_radiation!(column::ColumnVariables, model::PrimitiveEquation)
    longwave_radiation!(column, model.longwave_radiation, model)
end

## UNIFORM COOLING

export UniformCooling
Base.@kwdef struct UniformCooling{NF} <: AbstractLongwave
    time_scale::Second = Hour(16)
    temp_min::NF = 207.5
    temp_stratosphere::NF = 200
    time_scale_stratosphere::Second = Day(5)
end

UniformCooling(SG::SpectralGrid; kwargs...) = UniformCooling{SG.NF}(; kwargs...)
initialize!(scheme::UniformCooling, model::PrimitiveEquation) = nothing

function longwave_radiation!(
    column::ColumnVariables,
    scheme::UniformCooling,
    model::PrimitiveEquation,
)
    longwave_radiation!(column, scheme)
end

function longwave_radiation!(
    column::ColumnVariables{NF},
    scheme::UniformCooling,
) where NF
    (; temp, temp_tend) = column
    (; temp_min, temp_stratosphere) = scheme
    
    cooling = -inv(convert(NF, scheme.time_scale.value))
    τ⁻¹ = inv(scheme.time_scale_stratosphere.value)

    @inbounds for k in eachlayer(column)
        temp_tend[k] += temp[k] > temp_min ? cooling : (temp_stratosphere - temp[k])*τ⁻¹
    end
end

## JEEVANJEE TEMPERATURE FLUX

export JeevanjeeRadiation
"""
Temperature flux longwave radiation from Jeevanjee and Romps, 2018,
following Seeley and Wordsworth, 2023, eq (1)

    dF/dT = α*(T_t - T)

with `F` the upward temperature flux between two layers with temperature
difference `dT`, `α = 0.025 W/m²/K²` and `T_t = 200K` a prescribed tropopause
temperature. Flux into the lowermost layer is 0. Flux out of uppermost
layer also 0, but `dT/dt = (T_t - T)/τ` is added to relax the uppermost
layer towards the tropopause temperature `T_t` with time scale `τ = 24h`
(Seeley and Wordsworth, 2023 use 6h, which is unstable a low resolutions here).
Fields are
$(TYPEDFIELDS)"""
Base.@kwdef struct JeevanjeeRadiation{NF} <: AbstractLongwave
    "Radiative forcing constant (W/m²/K²)"
    α::NF = 0.025

    "Tropopause temperature [K]"
    temp_tropopause::NF = 200

    "Tropopause relaxation time scale to temp_tropopause"
    time_scale::Second = Hour(24)
end

JeevanjeeRadiation(SG::SpectralGrid; kwargs...) = JeevanjeeRadiation{SG.NF}(; kwargs...)
initialize!(scheme::JeevanjeeRadiation, model::PrimitiveEquation) = nothing

# function barrier
function longwave_radiation!(
    column::ColumnVariables,
    scheme::JeevanjeeRadiation,
    model::PrimitiveEquation,
)
    longwave_radiation!(column, scheme)
end

function longwave_radiation!(
    column::ColumnVariables{NF},
    scheme::JeevanjeeRadiation,
) where NF

    (; nlev, temp_tend) = column
    T = column.temp                 # to match Seeley, 2023 notation
    F = column.flux_temp_upward
    (; α, time_scale) = scheme
    Tₜ = scheme.temp_tropopause
    
    Fₖ::NF = 0                      # flux into lowermost layer from below

    # integrate from surface up, skipping surface (k=nlev+1) and top-of-atmosphere flux (k=1)
    @inbounds for k in nlev:-1:2    
        # Seeley and Wordsworth, 2023 eq (1)
        Fₖ += (T[k-1] - T[k]) * α * (Tₜ - T[k]) # upward flux from layer k into k-1
        F[k] += Fₖ                              # accumulate fluxes
    end

    # Relax the uppermost level towards prescribed "tropopause temperature"
    temp_tend[1] += (Tₜ - T[1])/time_scale.value
end