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
Computes downward and upward radiative fluxes through each layer using transmittance
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
        t = view(column.transmittance_shortwave, :, band)   # transmittance per layer in that band

        # DOWNWARD flux D
        D = solar_constant * cos_zenith             # top boundary condition of longwave flux
        column.flux_temp_downward[1] += D           # accumulate the top downward flux

        for k in 1:nlayers
            D *= t[k]                               # flux through layer k with transmittance t, radiative transfer
            column.flux_temp_downward[k+1] += D
        end

        # UPWARD flux U
        U = D*albedo                                # boundary condition at surface, reflection from albedo
        column.flux_temp_upward[nlayers+1] += U     # accumulate fluxes

        for k in nlayers:-1:1                       # integrate from surface up
            # Radiative transfer, e.g. Frierson et al. 2006, equation 6
            U *= t[k]                               # transmittance through layer k
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

    "[OPTION] Weight for √precip term [1]"
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

    "[OPTION] Stratocumulus cloud factor (SPEEDY clfact) [1]"
    stratocumulus_clfact::NF = 1.2

    "[OPTION] Zenith correction amplitude (SPEEDY azen) [1]"
    zenith_amplitude::NF = 1

    "[OPTION] Zenith correction exponent (SPEEDY nzen)"
    zenith_exponent::NF = 2

    "[OPTION] Absorptivity of dry air [per 10^5 Pa]"
    # Weighted visible + near-IR: 0.95*0.033 + 0.05*0.0 = 0.03135 (SPEEDY absdry, fband weights)
    absorptivity_dry_air::NF = 0.03135

    "[OPTION] Absorptivity of aerosols [per 10^5 Pa]"
    # Weighted visible + near-IR: 0.95*0.033 + 0.05*0.0 = 0.03135 (SPEEDY absaer, fband weights)
    absorptivity_aerosol::NF = 0.03135

    "[OPTION] Absorptivity of water vapor [per kg/kg per 10^5 Pa]"
    # Weighted visible + near-IR: 0.95*0.022 + 0.05*15.0*0.2 = 0.171 per g/kg → 1.7e-4 per kg/kg (SPEEDY abswv1, abswv2)
    absorptivity_water_vapor::NF = 0.00017

    "[OPTION] Base cloud absorptivity [per kg/kg per 10^5 Pa]"
    # Weighted visible band: 0.95*0.015 = 0.014 per g/kg → 1.4e-5 per kg/kg (SPEEDY abscl1)
    absorptivity_cloud_base::NF = 0.000014

    "[OPTION] Maximum cloud absorptivity [per 10^5 Pa]"
    # Weighted one-band scaling: 0.95*0.15 = 0.1425 → rounded to 0.14 (SPEEDY abscl2)
    absorptivity_cloud_limit::NF = 0.14

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
    # Contribution to cloud cover from precipitation, convert m/s to mm/day (Speedy multiplies by 86.4 because its rates are already in kg m⁻² s⁻¹ (= mm s⁻¹))
    P = wpcl*sqrt(min(pmcl, (86400 * (column.rain_rate_large_scale + column.rain_rate_convection) / 1000)))

    # Get Precipitation Top (kprtop) from the condensation scheme
    # This was previously calculated by other schemes and stored in column.cloud_top
    kprtop = column.cloud_top

    # Diagnose RH_term and Humidity Cloud Top (kcltop)
    (;relative_humidity_threshold_min, relative_humidity_threshold_max, specific_humidity_threshold_min) = radiation

    humidity_term_max = zero(P)
    
    kcltop = nlayers + 1 # Default to no cloud (below surface)

    # Loop from top layer (k=1) to the layer ABOVE the surface (nlayers-1)
    # to find the maximum humidity term and corresponding cloud top 
    for k in 1:(nlayers-1)
        humidity_k = humid[k]
        qsat = sat_humid[k]

        # Check if specific humidity is above the absolute threshold
        if humidity_k > specific_humidity_threshold_min
            # Calculate this layer's RH contribution
            relative_humidity_k = humidity_k / qsat
            rh_norm = max(0, (relative_humidity_k - relative_humidity_threshold_min) / (relative_humidity_threshold_max - relative_humidity_threshold_min))
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

    # --- Stratocumulus parameterization (stratocumulus_cloudcover) ---
    stratocumulus_cloudcover = zero(P)
    if radiation.use_stratocumulus
        # Compute static stability (GSEN)
        GSEN = dry_static_energy[nlayers] - dry_static_energy[nlayers-1]

        # Use tunable parameters from struct
        stability_min = radiation.stratocumulus_stability_min
        stability_max = radiation.stratocumulus_stability_max
        cover_max = radiation.stratocumulus_cover_max
        clfact = radiation.stratocumulus_clfact

        # F_ST: stability factor
        F_ST = max(0, min(1, (GSEN - stability_min) / (stability_max - stability_min)))

        # Stratocumulus cloud cover (stratocumulus_cloudcover) over ocean
        # Uses the single column CLC
        stratocumulus_cloudcover_ocean = F_ST * max(cover_max - clfact * cloud_cover, 0)
        
        # Over land, further modulate by surface RH
        RH_N = humid[nlayers] / sat_humid[nlayers]
        stratocumulus_cloudcover_land = stratocumulus_cloudcover_ocean * RH_N
        
        # Land-sea mask weighted stratocumulus cover
        stratocumulus_cloudcover = (1 - land_fraction) * stratocumulus_cloudcover_ocean + land_fraction * stratocumulus_cloudcover_land
    end

    # ---  Transmittance Calculation ---
    t = view(transmittance_shortwave, :, 1)
    if radiation.use_speedy_transmittance
        
        (; absorptivity_dry_air, absorptivity_aerosol, absorptivity_water_vapor,
        absorptivity_cloud_base, absorptivity_cloud_limit) = radiation

        sigma_levels = model.geometry.σ_levels_half
        surface_pressure = column.pres[end]  # This is in Pa

        # Zenith angle correction factor 
        azen = radiation.zenith_amplitude
        nzen = radiation.zenith_exponent

        # Zenith angle correction to (downward) absorptivity
        zenit_factor = 1 + azen * (1 - cos_zenith)^nzen

        # Cloud absorption term based on cloud base humidity (SPEEDY logic)
        q_base = nlayers > 1 ? humid[nlayers-1] : humid[nlayers]
        cloud_absorptivity_term = min(absorptivity_cloud_base * q_base,
                                    absorptivity_cloud_limit)

        for k in 1:nlayers
            q = humid[k]

            # Aerosol factor: use mid-level sigma, squared
            sigma_mid =  (sigma_levels[k] + sigma_levels[k+1]) / 2
            aerosol_factor = sigma_mid^2
            
            # Layer absorptivity (all humidity-based parameters are per kg/kg per 10^5 Pa)
            # Aerosol loading increases toward surface (proportional to σ²).
            layer_absorptivity = (absorptivity_dry_air +
                                absorptivity_aerosol * aerosol_factor +
                                absorptivity_water_vapor * q)

            # Add cloud absorption below the FINAL cloud top
            if k >= cloud_top
                layer_absorptivity += cloud_absorptivity_term * cloud_cover
            end

            # Compute differential optical depth with zenith correction
            # CRITICAL: Normalize pressure to 10^5 Pa since absorptivities are per 10^5 Pa
            delta_sigma = sigma_levels[k+1] - sigma_levels[k]
            normalized_pressure = surface_pressure / 100000   # Convert Pa to units of 10^5 Pa
            optical_depth = layer_absorptivity * delta_sigma * normalized_pressure * zenit_factor

            # Transmittance through layer k
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

    # At cloud top: reflect R fraction upward, transmit (1-R) downward
    if reflection_cloud_top < nlayers+1
        # Use single CLC value for reflection
        R = cloud_albedo * cloud_cover
        U_reflected = D * R
        D *= (1 - R)

        # Continue propagating downward below cloud top.
        for k in reflection_cloud_top:nlayers
            D *= t[k]
            flux_temp_downward[k+1] += D
        end
    end

    # --- Surface stratocumulus reflection ---
    # At the surface, apply stratocumulus reflection using stratocumulus_cloudcover and stratocumulus_cloud_albedo
    # Compute the reflected (upward) flux from stratocumulus at the surface
    U_stratocumulus = D * stratocumulus_cloud_albedo * stratocumulus_cloudcover
    D_surface = D - U_stratocumulus
    column.surface_shortwave_down = D_surface

    # Separate ocean and land albedo reflections.
    up_ocean = albedo_ocean * D_surface
    up_land  = albedo_land  * D_surface
    column.surface_shortwave_up_ocean = up_ocean
    column.surface_shortwave_up_land  = up_land

    # Computes grid-cell-average surface albedo and reflected shortwave flux
    albedo = (1 - land_fraction) * albedo_ocean + land_fraction * albedo_land
    U_surface_albedo = albedo * D_surface
    column.surface_shortwave_up = U_surface_albedo

    # Total upward flux = surface reflection + stratocumulus reflection.
    U = U_surface_albedo + U_stratocumulus
    
    # Upward beam
    flux_temp_upward[nlayers+1] += U
    
    # Propagate upward, applying transmittance. At cloud top, add the reflected flux.
    for k in nlayers:-1:1
        U *= t[k]
        # Add reflected flux at the correct layer
        U += k == reflection_cloud_top ? U_reflected : zero(U)
        flux_temp_upward[k] += U
    end

    # TOA outgoing shortwave
    column.outgoing_shortwave_radiation = U
    return nothing
end