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

"""
Uniform cooling following Pauluis and Garner, 2006. JAS. https://doi.org/10.1175/JAS3705.1
imposing a default temperature tendency of -1.5K/day (=1K/16hours for a `time_scale` of 16 hours)
on every level except for the stratosphere (diagnosed as `temp < temp_min`) where a
relaxation term with `time_scale_stratosphere` towards `temp_stratosphere` is applied.

    dT/dt = -1.5K/day for T > 207.5K else (200K-T) / 5 days

Fields are
$(TYPEDFIELDS)"""
@kwdef struct UniformCooling{NF} <: AbstractLongwave
    "[OPTION] time scale of cooling, default = -1.5K/day = -1K/16hrs"
    time_scale::Second = Hour(16)

    "[OPTION] temperature [K] below which stratospheric relaxation is applied"
    temp_min::NF = 207.5

    "[OPTION] target temperature [K] of stratospheric relaxation"
    temp_stratosphere::NF = 200

    "[OPTION] time scale of stratospheric relaxation"
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
        # Pauluis and Garner, 2006, eq (1) and (2)
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
@kwdef struct JeevanjeeRadiation{NF} <: AbstractLongwave
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

    (; nlayers, temp_tend) = column
    T = column.temp                 # to match Seeley, 2023 notation
    F = column.flux_temp_upward
    (; α, time_scale) = scheme
    Tₜ = scheme.temp_tropopause
    
    (; skin_temperature_sea, skin_temperature_land, land_fraction) = column

    # mix radiative fluxes for fractional land-sea mask if both available 
    land_available = isfinite(skin_temperature_land)
    sea_available = isfinite(skin_temperature_sea)

    # land fraction-weighted average of surface temperatures from land/sea
    local surface_temperature::NF = 0
    
    if land_available && sea_available
        surface_temperature = land_fraction*skin_temperature_land + (1-land_fraction)*skin_temperature_sea
    elseif land_available
        surface_temperature = skin_temperature_land
    elseif sea_available
        surface_temperature = skin_temperature_sea
    else    # fallback: use surface air temperature to have zero flux
        surface_temperature = T[end]
    end

    # extension to Jeevanjee: Include temperature flux between surface and lowermost air temperature
    local Fₖ::NF    # flux into lowermost layer from surface land/sea below
    Fₖ = (T[end] - surface_temperature) * α * (Tₜ - surface_temperature)
    F[end] += Fₖ

    # integrate from surface up, skipping surface (k=nlayers+1) and top-of-atmosphere flux (k=1)
    @inbounds for k in nlayers:-1:2
        # Seeley and Wordsworth, 2023 eq (1)
        Fₖ += (T[k-1] - T[k]) * α * (Tₜ - T[k]) # upward flux from layer k into k-1
        F[k] += Fₖ                              # accumulate fluxes
    end

    # Relax the uppermost level towards prescribed "tropopause temperature"
    temp_tend[1] += (Tₜ - T[1])/time_scale.value
end

# dummy one band radiation for now
export OneBandRadiation
@kwdef struct OneBandRadiation <: AbstractLongwave
    nbands::Int = 1
end

OneBandRadiation(SG::SpectralGrid; kwargs...) = OneBandRadiation(; kwargs...)
initialize!(scheme::OneBandRadiation, model::PrimitiveEquation) = nothing
longwave_radiation!(column::ColumnVariables, scheme::OneBandRadiation, model::PrimitiveEquation) = nothing