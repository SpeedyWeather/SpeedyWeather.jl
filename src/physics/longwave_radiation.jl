abstract type AbstractRadiation <: AbstractParameterization end
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
Uniform cooling following Paulius and Garner, 2006. JAS. https://doi.org/10.1175/JAS3705.1
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

    # extension to Jeevanjee: Include temperature flux between surface and lowermost air temperature
    # but zero flux if land/sea not available
    Fₖ_ocean = isfinite(skin_temperature_sea) ?
        (T[end] - skin_temperature_sea) * α * (Tₜ - skin_temperature_sea) : 
        zero(skin_temperature_sea)
    column.surface_longwave_up_ocean = Fₖ_ocean

    Fₖ_land = isfinite(skin_temperature_land) ?
        (T[end] - skin_temperature_land) * α * (Tₜ - skin_temperature_land) : 
        zero(skin_temperature_land)
    column.surface_longwave_up_land = Fₖ_land

    # land-sea mask weighted combined flux from land and ocean
    local Fₖ::NF
    Fₖ = (1-land_fraction)*Fₖ_ocean + land_fraction*Fₖ_land
    F[end] += Fₖ                        # accumulate into total flux
    column.surface_longwave_up = Fₖ     # store for output/diagnostics

    # integrate from surface up, skipping surface (k=nlayers+1) and top-of-atmosphere flux (k=1)
    @inbounds for k in nlayers:-1:2
        # Seeley and Wordsworth, 2023 eq (1)
        Fₖ += (T[k-1] - T[k]) * α * (Tₜ - T[k]) # upward flux from layer k into k-1
        F[k] += Fₖ                              # accumulate fluxes
    end

    # Relax the uppermost level towards prescribed "tropopause temperature"
    temp_tend[1] += (Tₜ - T[1])/time_scale.value

    # for diagnostic, use Fₖ as the outgoing longwave radiation although it's technically into the
    # uppermost layer from below (not out of it)
    column.outgoing_longwave_radiation = Fₖ
end

export NBandRadiation
@kwdef struct NBandRadiation{NF} <: AbstractRadiation
    nbands::Int = 1
    surface_longwave_emissivity::NF = 0.98
end

NBandRadiation(SG::SpectralGrid; kwargs...) = NBandRadiation{SG.NF}(; kwargs...)
initialize!(scheme::NBandRadiation, model::PrimitiveEquation) = nothing

function longwave_radiation!(
    column::ColumnVariables,
    radiation::NBandRadiation,
    model::PrimitiveEquation,
)

    (; nlayers) = column
    nbands = column.nbands_longwave                 # number of spectral bands
    ϵₛ = radiation.surface_longwave_emissivity      # surface emissivity
    NF = eltype(column)

    # precompute Stefan-Boltzmann emissions σT^4
    σ = model.atmosphere.stefan_boltzmann
    B = column.b                                    # reuse scratch vector b
    (; temp) = column
    @. B = σ*temp^4

    @inbounds for band in 1:nbands                  # loop over spectral bands

        # transmittance=exp(-optical_depth) of this band
        t = view(column.transmittance_longwave, :, band)

        # Dowward flux D
        D = zero(NF)                                # top boundary condition of longwave flux downward (no LW from space)
        for k in 1:nlayers
            # Radiative transfer, Frierson et al. 2006, equation 7 but use transmittance t
            # or Fortran SPEEDY documentation, eq. 44
            D = D*t[k] + (1 - t[k])*B[k]
            column.flux_temp_downward[k+1] += D
        end
        column.surface_longwave_down = D            # store for output/diagnostics

        # UPWARD flux U
        (; skin_temperature_sea, skin_temperature_land, land_fraction) = column

        # Flux between surface and lowermost air temperature but zero flux if land/sea not available
        R = (1 - ϵₛ)*D                              # surface longwave reflection
        U_ocean = isfinite(skin_temperature_sea) ? ϵₛ*σ*skin_temperature_sea^4 : zero(skin_temperature_sea)
        U_ocean += R
        column.surface_longwave_up_ocean = U_ocean
    
        U_land = isfinite(skin_temperature_land) ? ϵₛ*σ*skin_temperature_land^4 : zero(skin_temperature_land)
        U_land += R
        column.surface_longwave_up_land = U_land
    
        # land-sea mask weighted combined flux from land and ocean
        U = (1-land_fraction)*U_ocean + land_fraction*U_land
        column.flux_temp_upward[end] += U           # accumulate into total flux
        column.surface_longwave_up = U              # store for output/diagnostics

        for k in nlayers:-1:1                       # integrate from surface up
            # Radiative transfer, e.g. Frierson et al. 2006, equation 6
            U = U*t[k] + (1-t[k])*B[k]
            column.flux_temp_upward[k] += U         # accumulate that flux
        end

        # store outgoing longwave radiation (OLR) for diagnostics, accumulate over bands (reset when column is reset)
        column.outgoing_longwave_radiation += U
    end
end