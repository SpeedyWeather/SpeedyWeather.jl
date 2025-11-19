abstract type AbstractRadiation <: AbstractParameterization end
abstract type AbstractShortwave <: AbstractRadiation end
abstract type AbstractShortwaveRadiativeTransfer <: AbstractShortwave end

function get_nbands(radiation::Union{AbstractRadiation, Nothing})
    hasfield(typeof(radiation), :nbands) && return radiation.nbands
    return 0
end

## SHORTWAVE RADIATION FOR A FULLY TRANSPARENT ATMOSPHERE
export TransparentShortwave
struct TransparentShortwave <: AbstractShortwave end
TransparentShortwave(SG::SpectralGrid) = TransparentShortwave()

function variables(::TransparentShortwave)
    return (
        DiagnosticVariable(name=:surface_shortwave_down, dims=Grid2D(), desc="Surface shortwave radiation down", units="W/m^2"),
        DiagnosticVariable(name=:surface_shortwave_down, dims=Grid2D(), desc="Surface shortwave radiation down over ocean", units="W/m^2", namespace=:ocean),
        DiagnosticVariable(name=:surface_shortwave_down, dims=Grid2D(), desc="Surface shortwave radiation down over land", units="W/m^2", namespace=:land),
        DiagnosticVariable(name=:surface_shortwave_up,   dims=Grid2D(), desc="Surface shortwave radiation up",   units="W/m^2"),
        DiagnosticVariable(name=:outgoing_shortwave,     dims=Grid2D(), desc="TOA Shortwave radiation up",       units="W/m^2"),
        DiagnosticVariable(name=:cos_zenith,             dims=Grid2D(), desc="Cos zenith angle",                 units="1"),
        DiagnosticVariable(name=:albedo, dims=Grid2D(), desc="Albedo over ocean", units="1", namespace=:ocean),
        DiagnosticVariable(name=:albedo, dims=Grid2D(), desc="Albedo over land", units="1", namespace=:land),
    )
end

initialize!(::TransparentShortwave, ::PrimitiveEquation) = nothing

# function barrier
parameterization!(ij, diagn, progn, shortwave::TransparentShortwave, model) =
    shortwave_radiation!(ij, diagn, shortwave, model.planet, model.land_sea_mask)

function shortwave_radiation!(ij, diagn, shortwave::TransparentShortwave, planet, land_sea_mask)

    (; surface_shortwave_down, surface_shortwave_up) = diagn.physics
    (; outgoing_shortwave) = diagn.physics
    ssrd_ocean = diagn.physics.ocean.surface_shortwave_down
    ssrd_land = diagn.physics.land.surface_shortwave_down

    cos_zenith = diagn.physics.cos_zenith[ij]
    land_fraction = land_sea_mask[ij]
    albedo_ocean = diagn.physics.ocean.albedo[ij]
    albedo_land = diagn.physics.land.albedo[ij]
    S₀ = planet.solar_constant

    D = S₀ * cos_zenith             # top of atmosphere downward radiation
    surface_shortwave_down[ij] = D  # transparent atmosphere so same at surface (before albedo)

    # shortwave up is after albedo reflection, separated by ocean/land
    ssrd_ocean[ij] = albedo_ocean * D
    ssrd_land[ij]  = albedo_land * D

    # land-sea mask-weighted
    albedo = (1 - land_fraction)*albedo_ocean + land_fraction*albedo_land
    surface_shortwave_up[ij] = albedo * D

    # transparent also for reflected shortwave radiation travelling up
    outgoing_shortwave[ij] = surface_shortwave_up[ij]
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
struct OneBandShortwave{C, T, R} <: AbstractShortwave
    clouds::C
    transmittance::T
    radiative_transfer::R
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
@kwdef struct OneBandShortwaveRadiativeTransfer{NF} <: AbstractShortwaveRadiativeTransfer
    "[OPTION] Ozone absorption in upper stratosphere (W/m^2)"
    ozone_absorp_upper::NF = 0
    
    "[OPTION] Ozone absorption in lower stratosphere (W/m^2)"
    ozone_absorp_lower::NF = 0
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