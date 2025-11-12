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

"""
Core shortwave radiative transfer algorithm shared by OneBand and Spectral schemes.
Takes explicit albedo and ozone parameters to handle both single-value and band-specific cases.
"""

function _shortwave_radiative_transfer_core!(
    column::ColumnVariables,
    t,          # Transmittance array
    clouds,     # NamedTuple from clouds!
    ozone_upper, ozone_lower,  # Can be scalars or band-specific
    albedo_ocean, albedo_land, # Can be scalars or band-specific
    model::PrimitiveEquation,
)
    (; cloud_cover, cloud_top, stratocumulus_cover, cloud_albedo, stratocumulus_albedo) = clouds
    (; cos_zenith, land_fraction, nlayers, flux_temp_downward, flux_temp_upward) = column
   
    # Apply ozone absorption at TOA
    D_TOA = model.planet.solar_constant * cos_zenith
    D = D_TOA - ozone_upper - ozone_lower
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

@kwdef struct OneBandShortwaveRadiativeTransfer{NF} <: AbstractShortwaveRadiativeTransfer
    "[OPTION] Ozone absorption in upper stratosphere (W/m^2)"
    ozone_absorp_upper::NF = 0
    
    "[OPTION] Ozone absorption in lower stratosphere (W/m^2)"
    ozone_absorp_lower::NF = 0
end

# generator function
OneBandShortwaveRadiativeTransfer(SG::SpectralGrid; kwargs...) = OneBandShortwaveRadiativeTransfer{SG.NF}(; kwargs...)
initialize!(::OneBandShortwaveRadiativeTransfer, ::PrimitiveEquation) = nothing

"""$(TYPEDSIGNATURES)
One-band shortwave radiative transfer with cloud reflection and ozone absorption.
"""

function shortwave_radiative_transfer!(
    column::ColumnVariables,
    t,          # Transmittance array
    clouds,     # NamedTuple from clouds!
    radiation::OneBandShortwaveRadiativeTransfer,
    model::PrimitiveEquation
)
    # Use shared core with scalar albedos from column
    _shortwave_radiative_transfer_core!(
        column, t, clouds,
        radiation.ozone_absorp_upper,
        radiation.ozone_absorp_lower,
        column.albedo_ocean,  # scalar
        column.albedo_land,   # scalar
        model
    )
end

### DRAFT FOR N BAND SHORTWAVE RADIATION SCHEME
## NBAND SHORTWAVE RADIATION
export NBandShortwave

"""
    NBandShortwave <: AbstractShortwave

An N-band shortwave radiation scheme that extends OneBandShortwave to multiple spectral bands.
Each band is treated independently with its own transmittance and radiative transfer properties,
then weighted by solar spectrum fractions. This allows for realistic spectral dependence of
atmospheric absorption, cloud optical properties, and surface albedo.

Follows the concatenation approach where each band uses the same radiative transfer algorithm
but with band-specific optical properties. The total flux is computed as the weighted sum
over all spectral bands.

Fields are:
$(TYPEDFIELDS)"""
struct NBandShortwave{N, C, T, R, NF<:AbstractFloat} <: AbstractShortwave
    clouds::C
    transmittance::NTuple{N, T}
    radiative_transfer::NTuple{N, R}
    band_weights::NTuple{N, NF}
    wavelengths::NTuple{N, NF}
end

get_nbands(::NBandShortwave{N}) where N = N

"""$(TYPEDSIGNATURES)
Show NBandShortwave information in a structured format."""
function Base.show(io::IO, M::NBandShortwave{N}) where N
    println(io, "NBandShortwave{$N} <: AbstractShortwave")
    properties = propertynames(M)
    n = length(properties)
    for (i, key) in enumerate(properties)
        val = getfield(M, key)
        s = i == n ? "└" : "├"  # choose ending └ for last property
        p = i == n ? print : println
        if key == :transmittance || key == :radiative_transfer
            p(io, "$s $key: NTuple{$N, $(eltype(val))}")
        elseif key == :band_weights || key == :wavelengths
            p(io, "$s $key: $(typeof(val))")
        else
            p(io, "$s $key: $(typeof(val))")
        end
    end
end

# Basic constructor
NBandShortwave(SG::SpectralGrid; kwargs...) = NBandShortwave(SG, 1; kwargs...)

# Single-band constructor - uses N-band infrastructure for consistency
function NBandShortwave{1}(SG::SpectralGrid; 
    clouds = SpectralDiagnosticClouds{1}(SG),  # Use spectral clouds for consistency
    band_weights = (SG.NF(1.0),),  # All radiation in single band
    wavelengths = (SG.NF(0.5),))   # Representative solar wavelength
    
    # Use N-band infrastructure consistently
    transmittance = BackgroundShortwaveTransmittance(SG)
    radiative_transfer = SpectralShortwaveRadiativeTransfer{1}(SG)
    
    return NBandShortwave{1, typeof(clouds), typeof(transmittance), typeof(radiative_transfer), SG.NF}(
        clouds, (transmittance,), (radiative_transfer,), band_weights, wavelengths
    )
end

# Two-band constructor for visible + near-IR bands with full spectral dependence
function NBandShortwave{2}(SG::SpectralGrid; 
    clouds = SpectralDiagnosticClouds{2}(SG),  # Now spectral!
    band_weights = (SG.NF(0.5), SG.NF(0.5)),
    wavelengths = (SG.NF(0.5), SG.NF(1.5)))
    
    # Band 1: Visible (0.2-0.7 μm) - higher Rayleigh scattering, lower water vapor absorption
    transmittance_vis = BackgroundShortwaveTransmittance(SG;
        absorptivity_dry_air = 0.040,      # Higher Rayleigh scattering in visible
        absorptivity_water_vapor = 0.00005, # Lower water vapor absorption in visible
        absorptivity_aerosol = 0.040,      # Higher aerosol scattering in visible
        absorptivity_cloud_base = 0.000020, # Higher cloud scattering in visible
        absorptivity_cloud_limit = 0.20
    )
    
    # Band 2: Near-IR (0.7-4.0 μm) - lower Rayleigh scattering, higher water vapor absorption
    transmittance_nir = BackgroundShortwaveTransmittance(SG;
        absorptivity_dry_air = 0.020,      # Lower Rayleigh scattering in near-IR
        absorptivity_water_vapor = 0.0003, # Much higher water vapor absorption in near-IR
        absorptivity_aerosol = 0.025,      # Lower aerosol scattering in near-IR
        absorptivity_cloud_base = 0.000010, # Lower cloud scattering in near-IR
        absorptivity_cloud_limit = 0.08
    )
    
    # Use spectral radiative transfer with band-dependent surface albedo
    radiative_transfer_vis = SpectralShortwaveRadiativeTransfer{2}(SG)
    radiative_transfer_nir = SpectralShortwaveRadiativeTransfer{2}(SG)
    
    return NBandShortwave{2, typeof(clouds), typeof(transmittance_vis), typeof(radiative_transfer_vis), SG.NF}(
        clouds,
        (transmittance_vis, transmittance_nir),
        (radiative_transfer_vis, radiative_transfer_nir),
        band_weights,
        wavelengths
    )
end

# Generic constructor for any number of bands
function NBandShortwave(SG::SpectralGrid, N::Int; 
    clouds = DiagnosticClouds(SG),
    band_weights = ntuple(i -> SG.NF(1.0/N), N),  # Equal weights by default
    wavelengths = ntuple(i -> SG.NF(0.5), N))     # Default wavelengths
    
    # Create N identical components (user can customize later)
    transmittance = ntuple(_ -> BackgroundShortwaveTransmittance(SG), N)
    radiative_transfer = ntuple(_ -> OneBandShortwaveRadiativeTransfer(SG), N)
    
    return NBandShortwave{N, typeof(clouds), eltype(transmittance), eltype(radiative_transfer), SG.NF}(
        clouds, transmittance, radiative_transfer, band_weights, wavelengths
    )
end

function initialize!(radiation::NBandShortwave{N}, model::PrimitiveEquation) where N
    initialize!(radiation.clouds, model)
    for i in 1:N
        initialize!(radiation.transmittance[i], model) 
        initialize!(radiation.radiative_transfer[i], model)
    end
end

"""$(TYPEDSIGNATURES)
Calculate N-band shortwave radiation using the concatenation approach. Each spectral band
is computed independently using the same radiative transfer algorithm but with band-specific
optical properties. The total surface and TOA fluxes are computed as weighted sums over all bands.

The algorithm:
1. Diagnose clouds once (shared across all bands) 
2. For each band: compute transmittance and radiative transfer with band-specific properties
3. Weight band results by solar spectrum fractions
4. Sum to get total broadband fluxes"""
function shortwave_radiation!(
    column::ColumnVariables,
    radiation::NBandShortwave{N},
    model::PrimitiveEquation,
) where N
    
    # Initialize total fluxes - use proper types
    T = eltype(column.surface_shortwave_down)
    total_surface_down = zero(T)
    total_surface_up = zero(T)
    total_outgoing = zero(T)
    
    # Initialize flux profile accumulators
    flux_down_total = zeros(T, length(column.flux_temp_downward))
    flux_up_total = zeros(T, length(column.flux_temp_upward))
    
    # Loop over spectral bands
    for band in 1:N
        # Reset column values for clean band calculation
        column.surface_shortwave_down = zero(T)
        column.surface_shortwave_up = zero(T)
        column.outgoing_shortwave_radiation = zero(T)
        column.flux_temp_downward .= zero(T)
        column.flux_temp_upward .= zero(T)
        
        # Band-specific cloud properties
        clouds = clouds!(column, radiation.clouds, model, band)
        
        # Compute band-specific transmittance
        transmittance = transmittance!(column, clouds, radiation.transmittance[band], model)
        
        # Compute band-specific radiative transfer
        shortwave_radiative_transfer_bands!(column, transmittance, clouds, radiation.radiative_transfer[band], model, band)
        
        # Weight by solar spectrum fraction and accumulate
        weight = radiation.band_weights[band]
        total_surface_down += weight * column.surface_shortwave_down
        total_surface_up += weight * column.surface_shortwave_up
        total_outgoing += weight * column.outgoing_shortwave_radiation
        
        # Accumulate weighted flux profiles  
        @. flux_down_total += weight * column.flux_temp_downward
        @. flux_up_total += weight * column.flux_temp_upward
    end
    
    # Store final totals
    column.surface_shortwave_down = total_surface_down
    column.surface_shortwave_up = total_surface_up
    column.outgoing_shortwave_radiation = total_outgoing
    column.flux_temp_downward .= flux_down_total
    column.flux_temp_upward .= flux_up_total
    
    return nothing
end

"""
Spectral radiative transfer with band-dependent surface albedo
"""

export SpectralShortwaveRadiativeTransfer

"""
    SpectralShortwaveRadiativeTransfer{N} <: AbstractShortwaveRadiativeTransfer

N-band shortwave radiative transfer scheme with wavelength-dependent surface albedos
and ozone absorption. Allows realistic treatment of spectral dependence in surface
reflection properties.

$(TYPEDFIELDS)
"""
@kwdef struct SpectralShortwaveRadiativeTransfer{N, NF} <: AbstractShortwaveRadiativeTransfer
    "[OPTION] Ozone absorption in upper stratosphere for each band (W/m^2)"
    ozone_absorp_upper::NTuple{N, NF} = ntuple(i -> NF(0), N)
    
    "[OPTION] Ozone absorption in lower stratosphere for each band (W/m^2)"
    ozone_absorp_lower::NTuple{N, NF} = ntuple(i -> NF(0), N)
    
    "[OPTION] Ocean albedo for each spectral band [1]"
    albedo_ocean::NTuple{N, NF} = ntuple(i -> NF(0.06), N)
    
    "[OPTION] Land albedo for each spectral band [1]"
    albedo_land::NTuple{N, NF} = ntuple(i -> NF(0.15), N)
end

# Constructor with realistic spectral albedos
function SpectralShortwaveRadiativeTransfer{2}(SG::SpectralGrid;
    # Ocean: lower albedo in visible, higher in near-IR (typical for water)
    albedo_ocean = (SG.NF(0.05), SG.NF(0.08)),    # [visible, near-IR]
    # Land: much higher near-IR albedo due to vegetation reflectance
    albedo_land = (SG.NF(0.12), SG.NF(0.25)),     # [visible, near-IR]
    kwargs...)
    
    return SpectralShortwaveRadiativeTransfer{2, SG.NF}(;
        albedo_ocean=albedo_ocean,
        albedo_land=albedo_land,
        kwargs...
    )
end

SpectralShortwaveRadiativeTransfer(SG::SpectralGrid; kwargs...) = SpectralShortwaveRadiativeTransfer{1}(SG; kwargs...)
SpectralShortwaveRadiativeTransfer{N}(SG::SpectralGrid; kwargs...) where N = SpectralShortwaveRadiativeTransfer{N, SG.NF}(; kwargs...)

initialize!(::SpectralShortwaveRadiativeTransfer, ::PrimitiveEquation) = nothing

"""$(TYPEDSIGNATURES)
Band-specific shortwave radiative transfer with wavelength-dependent surface albedo."""
function shortwave_radiative_transfer_bands!(
    column::ColumnVariables,
    t,          # Transmittance array for this band
    clouds,     # NamedTuple from clouds! (already band-specific)
    radiation::SpectralShortwaveRadiativeTransfer{N},
    model::PrimitiveEquation,
    band::Int = 1  # Which spectral band to compute
) where N
    # Use shared core with band-specific parameters
    _shortwave_radiative_transfer_core!(
        column, t, clouds,
        radiation.ozone_absorp_upper[band],  # Band-specific
        radiation.ozone_absorp_lower[band],  # Band-specific
        radiation.albedo_ocean[band],        # Band-specific
        radiation.albedo_land[band],         # Band-specific
        model
    )
end