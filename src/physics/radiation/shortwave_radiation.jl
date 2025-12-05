abstract type AbstractRadiation <: AbstractParameterization end
abstract type AbstractShortwave <: AbstractRadiation end
abstract type AbstractShortwaveRadiativeTransfer <: AbstractShortwave end

"""$(TYPEDSIGNATURES)
Get the number of spectral bands for a radiation scheme. Returns the `nbands` field
if it exists on the radiation scheme type, otherwise returns 0."""
function get_nbands(radiation::Union{AbstractRadiation, Nothing})
    hasfield(typeof(radiation), :nbands) && return radiation.nbands
    return 0
end

"""$(TYPEDSIGNATURES)
Function barrier for shortwave radiation. Dispatches to the appropriate shortwave
radiation calculation based on the scheme type."""
function shortwave_radiation!(column::ColumnVariables, model::PrimitiveEquation)
    shortwave_radiation!(column, model.shortwave_radiation, model)
end

"""$(TYPEDSIGNATURES)
No-op for cases where shortwave radiation is not included in the model."""
shortwave_radiation!(::ColumnVariables, ::Nothing, ::PrimitiveEquation) = nothing

## SHORTWAVE RADIATION FOR A FULLY TRANSPARENT ATMOSPHERE
export TransparentShortwave

"""
    TransparentShortwave <: AbstractShortwave

A shortwave radiation scheme for a fully transparent atmosphere where radiation
passes through without absorption or scattering."""
struct TransparentShortwave <: AbstractShortwave end
TransparentShortwave(SG::SpectralGrid) = TransparentShortwave()
initialize!(::TransparentShortwave, ::PrimitiveEquation) = nothing

"""$(TYPEDSIGNATURES)
Dispatch to shortwave radiation calculation with planet information."""
function shortwave_radiation!(
    column::ColumnVariables,
    scheme::TransparentShortwave,
    model::PrimitiveEquation,
)
    shortwave_radiation!(column, scheme, model.planet)
end

"""$(TYPEDSIGNATURES)
Calculate shortwave radiation for a transparent atmosphere. Radiation equals
solar constant times cosine of zenith angle at the top of the atmosphere and
surface. Surface reflection is determined by ocean and land albedos."""
function shortwave_radiation!(
    column::ColumnVariables,
    scheme::TransparentShortwave,
    planet::AbstractPlanet,
)
    (; cos_zenith, land_fraction, albedo_ocean, albedo_land) = column
    (; solar_constant) = planet

    D = solar_constant * cos_zenith         # top of atmosphere downward radiation
    column.surface_shortwave_down = D       # transparent atmosphere so same at surface (before albedo)

    # shortwave up is after albedo reflection, separated by ocean/land
    column.surface_shortwave_up_ocean = albedo_ocean * D
    column.surface_shortwave_up_land = albedo_land * D

    # land-sea mask-weighted, transparent shortwave so surface = outgoing
    column.surface_shortwave_up = (1 - land_fraction)*column.surface_shortwave_up_ocean +
                                            land_fraction*column.surface_shortwave_up_land
    column.outgoing_shortwave_radiation = column.surface_shortwave_up

    return nothing
end

export OneBandShortwave

"""
    OneBandShortwave <: AbstractShortwave

A one-band shortwave radiation scheme with diagnostic clouds following Fortran SPEEDY
documentation, section B4. Implements cloud top detection and cloud fraction calculation
based on relative humidity and precipitation, with radiative transfer through clear sky and cloudy layers.

Cloud cover is calculated as a combination of relative humidity and precipitation contributions,
and a cloud albedo is applied to the downward beam. Fields and options are

$(TYPEDFIELDS)"""
@parameterized struct OneBandShortwave{C, T, R} <: AbstractShortwave
    @component clouds::C
    @component transmittance::T
    @component radiative_transfer::R
end

# primitive wet model version
OneBandShortwave(SG::SpectralGrid) = OneBandShortwave(
    DiagnosticClouds(SG),
    BackgroundShortwaveTransmittance(SG),
    OneBandShortwaveRadiativeTransfer(SG),
)

# primitive dry model version
export OneBandGreyShortwave
OneBandGreyShortwave(SG::SpectralGrid) = OneBandShortwave(
    NoClouds(SG),
    TransparentShortwaveTransmittance(SG),
    OneBandShortwaveRadiativeTransfer(SG),
)

get_nbands(::OneBandShortwave) = 1

function Base.show(io::IO, M::OneBandShortwave)
    println(io, "OneBandShortwave <: AbstractShortwave")
    properties = propertynames(M)
    n = length(properties)
    for (i, key) in enumerate(properties)
        val = getfield(M, key)
        s = i == n ? "└" : "├"  # choose ending └ for last property
        p = i == n ? print : println
        p(io, "$s $key: $(typeof(val))")
    end
end

# initialize one after another
function initialize!(radiation::OneBandShortwave, model::PrimitiveEquation)
    initialize!(radiation.clouds, model)
    initialize!(radiation.transmittance, model)
    initialize!(radiation.radiative_transfer, model)
end

"""$(TYPEDSIGNATURES)
Calculate shortwave radiation using the one-band scheme with diagnostic clouds.
Computes cloud cover fraction from relative humidity and precipitation, then
integrates downward and upward radiative fluxes accounting for cloud albedo effects."""
function shortwave_radiation!(
    column::ColumnVariables,
    radiation::OneBandShortwave,
    model::PrimitiveEquation,
)
    clouds = clouds!(column, radiation.clouds, model)
    t = transmittance!(column, clouds, radiation.transmittance, model)
    shortwave_radiative_transfer!(column, t, clouds, radiation.radiative_transfer, model)
end

export OneBandShortwaveRadiativeTransfer
@parameterized @kwdef struct OneBandShortwaveRadiativeTransfer{NF} <: AbstractShortwaveRadiativeTransfer
    "[OPTION] Ozone absorption in upper stratosphere (W/m^2)"
    @param ozone_absorp_upper::NF = 0 (bounds=Nonnegative,)
    
    "[OPTION] Ozone absorption in lower stratosphere (W/m^2)"
    @param ozone_absorp_lower::NF = 0 (bounds=Nonnegative,)
end

# generator function
OneBandShortwaveRadiativeTransfer(SG::SpectralGrid; kwargs...) = OneBandShortwaveRadiativeTransfer{SG.NF}(; kwargs...)
initialize!(::OneBandShortwaveRadiativeTransfer, ::PrimitiveEquation) = nothing

# function barrier to unpack model
function shortwave_radiative_transfer!(
    column::ColumnVariables,
    t,          # Transmittance array
    clouds,     # NamedTuple from clouds!
    radiation::OneBandShortwaveRadiativeTransfer,
    model::PrimitiveEquation
)
    shortwave_radiative_transfer!(column, t, clouds, radiation, model.planet)
end

"""$(TYPEDSIGNATURES)
One-band shortwave radiative transfer with cloud reflection and ozone absorption."""
function shortwave_radiative_transfer!(
    column::ColumnVariables,
    t,          # Transmittance array
    clouds,     # NamedTuple from clouds!
    radiation::OneBandShortwaveRadiativeTransfer,
    planet::AbstractPlanet,
)
    (; ozone_absorp_upper, ozone_absorp_lower) = radiation
    (; albedo_ocean, albedo_land) = column
    (; cloud_cover, cloud_top, stratocumulus_cover, cloud_albedo, stratocumulus_albedo) = clouds
    (; cos_zenith, land_fraction, nlayers, flux_temp_downward, flux_temp_upward) = column
   
    # Apply ozone absorption at TOA
    D_TOA = planet.solar_constant * cos_zenith
    D = max(zero(D_TOA), D_TOA - ozone_absorp_upper - ozone_absorp_lower)
    flux_temp_downward[1] += D

    # Clear sky portion until cloud top
    for k in 1:(cloud_top - 1)
        D *= t[k]
        flux_temp_downward[k+1] += D
    end

    # Cloud reflection at cloud top
    U_reflected = zero(D)
    if cloud_top <= nlayers
        R = cloud_albedo * cloud_cover
        U_reflected = D * R
        D *= (1 - R)

        for k in cloud_top:nlayers
            D *= t[k]
            flux_temp_downward[k+1] += D
        end
    end

    # Surface stratocumulus reflection
    U_stratocumulus = D * stratocumulus_albedo * stratocumulus_cover
    D_surface = D - U_stratocumulus
    column.surface_shortwave_down = D_surface

    # Surface albedo reflections (can handle both scalar and band-specific albedos)
    up_ocean = albedo_ocean * D_surface
    up_land  = albedo_land * D_surface
    column.surface_shortwave_up_ocean = up_ocean
    column.surface_shortwave_up_land  = up_land

    # Weighted surface albedo and reflected flux
    albedo = (1 - land_fraction) * albedo_ocean + land_fraction * albedo_land
    U_surface_albedo = albedo * D_surface
    column.surface_shortwave_up = U_surface_albedo

    U = U_surface_albedo + U_stratocumulus
    
    # Upward beam
    flux_temp_upward[nlayers+1] += U
    for k in nlayers:-1:1
        U *= t[k]
        U += k == cloud_top ? U_reflected : zero(U)
        flux_temp_upward[k] += U
    end

    column.outgoing_shortwave_radiation = U
    return nothing
end