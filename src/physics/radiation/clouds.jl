abstract type AbstractShortwaveClouds <: AbstractShortwave end

export NoClouds
struct NoClouds <: AbstractShortwaveClouds end
NoClouds(SG::SpectralGrid) = NoClouds()
initialize!(clouds::NoClouds, ::AbstractModel) = nothing
function clouds!(
    column::ColumnVariables{NF},
    ::NoClouds,
    model::AbstractModel,
) where NF
    cloud_cover = zero(NF)          # no cloud cover
    cloud_top = column.nlayers+1    # below surface

    return (    # NamedTuple
        cloud_cover=cloud_cover,
        cloud_albedo=zero(NF),
        cloud_top=cloud_top,
        stratocumulus_cover=zero(NF),
        stratocumulus_albedo=zero(NF),
    )
end

export DiagnosticClouds
@kwdef struct DiagnosticClouds{NF} <: AbstractShortwaveClouds
    "[OPTION] Relative humidity threshold for cloud cover = 0 [1]."
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

    "[OPTION] Stratocumulus cloud albedo (surface reflection/absorption) [1]"
    stratocumulus_albedo::NF = 0.5

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
end

DiagnosticClouds(SG::SpectralGrid; kwargs...) = DiagnosticClouds{SG.NF}(; kwargs...)
initialize!(clouds::DiagnosticClouds, model::AbstractModel) = nothing
function clouds!(
    column::ColumnVariables,
    clouds::DiagnosticClouds,
    model::AbstractModel,
)

    # ---- pull only fields that exist on ColumnVariables ----
    (; nlayers, land_fraction, humid, sat_humid, dry_static_energy) = column

    # Weight and max precip for cloud cover calculation
    wpcl = clouds.precipitation_weight
    pmcl = clouds.precipitation_max
    # Contribution to cloud cover from precipitation, convert m/s to mm/day (Speedy multiplies by 86.4 because its rates are already in kg m⁻² s⁻¹ (= mm s⁻¹))
    P = wpcl*sqrt(min(pmcl, (86400 * (column.rain_rate_large_scale + column.rain_rate_convection) / 1000)))

    # Get Precipitation Top (kprtop) from the condensation scheme
    # This was previously calculated by other schemes and stored in column.cloud_top
    kprtop = column.cloud_top

    # Diagnose RH_term and Humidity Cloud Top (kcltop)
    (;relative_humidity_threshold_min, relative_humidity_threshold_max, specific_humidity_threshold_min) = clouds

    humidity_term_kcltop = zero(P)
    
    kcltop = nlayers + 1 # Default to no cloud (below surface)

    # Following SPEEDY documentation B.4: cloud_top starts at nlayers+1 (no cloud),
    # then is moved up to the highest layer k where condensation occurs.
    # Loop from top layer (k=1) to find the TOPMOST layer where RH ≥ 95% (condensation threshold).
    # This matches ImplicitCondensation which also triggers at 95% RH.
    for k in 1:(nlayers-1)
        humidity_k = humid[k]
        qsat = sat_humid[k]

        # Check if specific humidity is above the absolute threshold
        if humidity_k > specific_humidity_threshold_min && qsat > 0
            # Calculate this layer's RH (avoid division by zero in dry models)
            relative_humidity_k = humidity_k / qsat
            
            # Check if RH exceeds the condensation threshold 
            if relative_humidity_k >= relative_humidity_threshold_min
                # Cloud cover increases smoothly from RH=95% to RH=100%
                rh_norm = max(0, (relative_humidity_k - relative_humidity_threshold_min) / (relative_humidity_threshold_max - relative_humidity_threshold_min))
                humidity_term_k = min(1, rh_norm)^2

                # Set kcltop to the TOPMOST (smallest k) layer with condensation
                kcltop = k
                humidity_term_kcltop = humidity_term_k
                break
            end
        end
    end
    
    # Calculate the Single Column Cloud Cover (CLC)
    # This uses the humidity term from the cloud top layer (SPEEDY B.4)
    cloud_cover = min(1, P + humidity_term_kcltop)

    # The final cloud top is the minimum (highest in alt) of the two
    cloud_top = min(kcltop, kprtop)

    # Update the column's cloud_top with this final, correct value
    column.cloud_top = cloud_top

    # --- Stratocumulus parameterization (stratocumulus_cover) ---
    stratocumulus_cover = zero(P)
    if clouds.use_stratocumulus
        # Compute static stability (GSEN)
        GSEN = dry_static_energy[nlayers] - dry_static_energy[nlayers-1]

        # Use tunable parameters from struct
        stability_min = clouds.stratocumulus_stability_min
        stability_max = clouds.stratocumulus_stability_max
        cover_max = clouds.stratocumulus_cover_max
        clfact = clouds.stratocumulus_clfact

        # F_ST: stability factor
        F_ST = max(0, min(1, (GSEN - stability_min) / (stability_max - stability_min)))

        # Stratocumulus cloud cover (stratocumulus_cover) over ocean
        # Uses the single column CLC
        stratocumulus_cover_ocean = F_ST * max(cover_max - clfact * cloud_cover, 0)
        
        # Over land, further modulate by surface RH
        RH_N = sat_humid[nlayers] > 0 ? humid[nlayers] / sat_humid[nlayers] : 0
        stratocumulus_cover_land = stratocumulus_cover_ocean * RH_N
        
        # Land-sea mask weighted stratocumulus cover
        stratocumulus_cover = (1 - land_fraction) * stratocumulus_cover_ocean + land_fraction * stratocumulus_cover_land
    end

    return (    # NamedTuple
        cloud_cover=cloud_cover,
        cloud_albedo=clouds.cloud_albedo,
        cloud_top=cloud_top,
        stratocumulus_cover=stratocumulus_cover,
        stratocumulus_albedo=clouds.stratocumulus_albedo,
    )
end