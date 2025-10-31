abstract type AbstractShortwave <: AbstractRadiation end

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

"""$(TYPEDSIGNATURES)
Calculate shortwave radiation for an N-band (multi-spectral) radiation scheme.
Computes downward and upward radiative fluxes through each layer using optical depth
and radiative transfer equations, accumulating fluxes and outgoing shortwave radiation."""
function shortwave_radiation!(
    column::ColumnVariables,
    scheme::NBandRadiation,
    model::PrimitiveEquation,
)

    (; nlayers, cos_zenith, albedo) = column
    nbands = column.nbands_shortwave                # number of spectral bands
    (; solar_constant) = model.planet

    @inbounds for band in 1:nbands                  # loop over spectral bands
        dτ = view(column.optical_depth_shortwave, :, band)   # differential optical depth per layer of that band

        # DOWNWARD flux D
        D = solar_constant * cos_zenith             # top boundary condition of longwave flux
        column.flux_temp_downward[1] += D           # accumulate the top downward flux

        for k in 1:nlayers
            D -= dτ[k]*D                            # flux through layer k with optical depth dτ, radiative transfer
            column.flux_temp_downward[k+1] += D
        end

        # UPWARD flux U
        U = D*albedo                                # boundary condition at surface, reflection from albedo
        column.flux_temp_upward[nlayers+1] += U     # accumulate fluxes

        for k in nlayers:-1:1                       # integrate from surface up
            # Radiative transfer, e.g. Frierson et al. 2006, equation 6
            U -= dτ[k]*U                            # negative because we integrate from surface up in -τ direction
            column.flux_temp_upward[k] += U         # accumulate that flux
        end

        # store outgoing shortwave radiation (OSR) for diagnostics, accumulate over bands (reset when column is reset)
        column.outgoing_shortwave_radiation += U
    end
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
@kwdef struct OneBandShortwave{NF} <: AbstractShortwave
    "[OPTION] Relative humidity threshold for cloud cover = 0 [1]"
    relative_humidity_threshold_min::NF = 0.3

    "[OPTION] Relative humidity threshold for cloud cover = 1 [1]"
    relative_humidity_threshold_max::NF = 1

    "[OPTION] Specific humidity threshold for cloud cover [Kg/kg]"
    specific_humidity_threshold_min::NF = 0.0002

    "[OPTION] Weight for √precip term [1] at 1 mm/day"
    precipitation_weight::NF = 0.2

    "[OPTION] Cap on precip contributing to cloud cover [mm/day]"
    precipitation_max::NF = 10

    "[OPTION] Cloud albedo for visible band at CLC=1 [1]"
    cloud_albedo::NF = 0.43

    "[OPTION] Use cloud top reflection?"
    cloud_top_reflection::Bool = true

    "[OPTION] Ozone absorption in upper stratosphere (W/m^2)"
    ozone_absorp_upper::NF = 0
    
    "[OPTION] Ozone absorption in lower stratosphere (W/m^2)"
    ozone_absorp_lower::NF = 0

    "[OPTION] Stratocumulus cloud albedo (surface reflection/absorption) [1]"
    stratocumulus_cloud_albedo::NF = 0.5

    "[OPTION] Static stability lower threshold for stratocumulus (GSEN, called GSES0 in the paper) [J/kg]"
    stratocumulus_stability_min::NF = 0.25

    "[OPTION] Static stability upper threshold for stratocumulus (GSEN, called GSES1 in the paper) [J/kg]"
    stratocumulus_stability_max::NF = 0.40

    "[OPTION] Maximum stratocumulus cloud cover (called CLSMAX in the paper) [1]"
    stratocumulus_cover_max::NF = 0.60

    "[OPTION] Enable stratocumulus cloud parameterization?"
    use_stratocumulus::Bool = true

    "[OPTION] Absorptivity of dry air (effective OneBand) [1]"
    # (0.95 * 0.033) + (0.05 * 0.0) = 0.03135
    absorptivity_dry_air::NF = 0.03135

    "[OPTION] Absorptivity of aerosols (effective OneBand) [1]"
    # (0.95 * 0.033) + (0.05 * 0.0) = 0.03135
    absorptivity_aerosol::NF = 0.03135

    "[OPTION] Absorptivity of water vapor (effective OneBand) [1]"
    # (0.95 * 22.0) + (0.05 * 15000.0) = 770.9
    # Scaled by 1000 for kg/kg units
    absorptivity_water_vapor::NF = 770.9

    "[OPTION] Base cloud absorptivity (effective OneBand) [1]"
    # (0.95 * 15.0) + (0.05 * 0.0) = 14.25
    # Scaled by 1000 for kg/kg units
    absorptivity_cloud_base::NF = 14.25

    "[OPTION] Maximum cloud absorptivity (effective OneBand) [1]"
    # (0.95 * 0.15) + (0.05 * 0.0) = 0.1425
    # This parameter (abscl2) is unitless, so it stays the same
    absorptivity_cloud_limit::NF = 0.1425

    "[OPTION] Use SPEEDY-style shortwave transmittance scheme?"
    use_speedy_transmittance::Bool = true
end

# generator function
OneBandShortwave(SG::SpectralGrid; kwargs...) = OneBandShortwave{SG.NF}(; kwargs...)
get_nbands(::OneBandShortwave) = 1
initialize!(::OneBandShortwave, ::PrimitiveEquation) = nothing

"""$(TYPEDSIGNATURES)
Calculate shortwave radiation using the one-band scheme with diagnostic clouds.
Computes cloud cover fraction from relative humidity and precipitation, then
integrates downward and upward radiative fluxes accounting for cloud albedo effects."""
function shortwave_radiation!(
    column::ColumnVariables,
    radiation::OneBandShortwave,
    model::PrimitiveEquation,
)
    # ---- pull only fields that exist on ColumnVariables ----
    (; nlayers, cos_zenith, land_fraction,
       humid, sat_humid, albedo_ocean, albedo_land,
       transmittance_shortwave, flux_temp_downward, flux_temp_upward, dry_static_energy) = column

     
    # Weight and max precip for cloud cover calculation
    wpcl = radiation.precipitation_weight
    pmcl = radiation.precipitation_max
    # Contribution to cloud cover from precipitation, convert m/s to mm/day
    P = wpcl*sqrt(min(pmcl, (86400 * (column.rain_rate_large_scale + column.rain_rate_convection) / 1000)))

    # Get Precipitation Top (kprtop) from the condensation scheme
    # This was previously calculated by other schemes and stored in column.cloud_top
    kprtop = column.cloud_top

    # Diagnose RH_term and Humidity Cloud Top (kcltop)
    rh_min = radiation.relative_humidity_threshold_min
    rh_max = radiation.relative_humidity_threshold_max
    q_min_kg_kg = radiation.specific_humidity_threshold_min

    humidity_term_max = zero(P)
    kcltop = nlayers + 1 # Default to no cloud (below surface)

    # Loop from top layer (k=1) to the layer ABOVE the surface (nlayers-1)
    for k in 1:(nlayers-1)
        q = humid[k]
        qsat = sat_humid[k]

        # Check if specific humidity is above the absolute threshold
        if q > q_min_kg_kg
            # Calculate this layer's RH contribution
            rh = q / qsat
            rh_norm = max(0, (rh - rh_min) / (rh_max - rh_min))
            humidity_term_k = min(1, rh_norm)^2

            # If this is the new maximum, update the max term and set humidity top
            if humidity_term_k > humidity_term_max
                humidity_term_max = humidity_term_k
                kcltop = k  # This layer is now the humidity cloud top
            end
        end
    end
    
    # Calculate the Single Column Cloud Cover (CLC)
    # This uses the humidity term from the kcltop layer
    cloud_cover = min(1, P + humidity_term_max)

    # The final cloud top is the minimum (highest in alt) of the two
    cloud_top = min(kcltop, kprtop)

    # Update the column's cloud_top with this final, correct value
    column.cloud_top = cloud_top

    # Local variable for reflection (can be disabled)
    reflection_cloud_top = radiation.cloud_top_reflection ? cloud_top : nlayers + 1

    # --- Stratocumulus parameterization (CLS) ---
    CLS = zero(P)
    if radiation.use_stratocumulus
        # Compute static stability (GSEN)
        GSEN = dry_static_energy[nlayers] - dry_static_energy[nlayers-1]

        # Use tunable parameters from struct
        stability_min = radiation.stratocumulus_stability_min
        stability_max = radiation.stratocumulus_stability_max
        cover_max = radiation.stratocumulus_cover_max
        clfact = 1.2  # Added missing factor from SPEEDY

        # F_ST: stability factor
        F_ST = max(0, min(1, (GSEN - stability_min) / (stability_max - stability_min)))

        # Stratocumulus cloud cover (CLS) over ocean
        # Uses the single column CLC
        CLS_ocean = F_ST * max(cover_max - clfact * cloud_cover, 0)
        
        # Over land, further modulate by surface RH
        RH_N = humid[nlayers] / sat_humid[nlayers]
        CLS_land = CLS_ocean * RH_N
        
        # Land-sea mask weighted stratocumulus cover
        CLS = (1 - land_fraction) * CLS_ocean + land_fraction * CLS_land
    end

    # ---  Transmittance Calculation ---
    t = view(transmittance_shortwave, :, 1)
    if radiation.use_speedy_transmittance
        
        (; absorptivity_dry_air, absorptivity_aerosol, absorptivity_water_vapor,
        absorptivity_cloud_base, absorptivity_cloud_limit) = radiation

        sigma_levels   = model.geometry.σ_levels_half
        surface_pressure = column.pres[end] 
        reference_surface_pressure = model.atmosphere.pres_ref # Corrected path

        # Cloud absorption term based on cloud base humidity (SPEEDY logic)
        q_base = humid[nlayers-1] # Assumes nlayers-1 is the cloud base layer
        cloud_absorptivity_term = min(absorptivity_cloud_base * q_base,
                                    absorptivity_cloud_limit)

        for k in 1:nlayers
            aerosol_factor = sigma_levels[k]^0.5
            layer_absorptivity = absorptivity_dry_air +
                                absorptivity_aerosol * aerosol_factor^2 +
                                absorptivity_water_vapor * humid[k]

            # Add cloud absorption below the FINAL cloud top
            if k >= cloud_top
                # Use the single CLC value and the pre-calculated term
                layer_absorptivity += cloud_absorptivity_term * cloud_cover
            end

            # compute differential optical depth
            delta_sigma = sigma_levels[k+1] - sigma_levels[k]
            optical_depth = layer_absorptivity * delta_sigma * surface_pressure / reference_surface_pressure
        

            # transmittance through layer k
            t[k] = exp(-optical_depth)
        end
    end

    # Downward beam
    (; cloud_albedo, ozone_absorp_upper, ozone_absorp_lower, stratocumulus_cloud_albedo) = radiation
   
    # Apply ozone absorption at TOA as subtracted fluxes 
    D_TOA = model.planet.solar_constant * cos_zenith
    D = D_TOA - ozone_absorp_upper - ozone_absorp_lower
    flux_temp_downward[1] += D

    # Clear sky portion until cloud top
    for k in 1:(reflection_cloud_top - 1)
        D *= t[k]
        flux_temp_downward[k+1] += D
    end

    # Cloud reflection (cloud top)
    U_reflected = zero(D)
    if reflection_cloud_top < nlayers+1

        # Use single CLC value for reflection
        R = cloud_albedo * cloud_cover
        U_reflected = D * R
        D *= (1 - R)
        for k in reflection_cloud_top:nlayers
            D *= t[k]
            flux_temp_downward[k+1] += D
        end
    end

    # --- Surface stratocumulus reflection ---
    # At the surface, apply stratocumulus reflection using CLS and stratocumulus_cloud_albedo
    # Compute the reflected (upward) flux from stratocumulus at the surface
    U_stratocumulus = D * stratocumulus_cloud_albedo * CLS
    D_surface = D - U_stratocumulus
    column.surface_shortwave_down = D_surface

    up_ocean = albedo_ocean * D_surface
    up_land  = albedo_land  * D_surface
    column.surface_shortwave_up_ocean = up_ocean
    column.surface_shortwave_up_land  = up_land

    # Computes grid-cell-average surface albedo and reflected shortwave flux
    albedo = (1 - land_fraction) * albedo_ocean + land_fraction * albedo_land
    U_surface_albedo = albedo * D_surface
    column.surface_shortwave_up = U_surface_albedo

    U = U_surface_albedo + U_stratocumulus
    
    # Upward beam
    flux_temp_upward[nlayers+1] += U

    for k in nlayers:-1:1
        U *= t[k]
        # Add reflected flux at the correct layer
        U += k == reflection_cloud_top ? U_reflected : zero(U) # <--- THIS IS CORRECT
        flux_temp_upward[k] += U
    end

    # TOA outgoing shortwave
    column.outgoing_shortwave_radiation = U
    return nothing
end