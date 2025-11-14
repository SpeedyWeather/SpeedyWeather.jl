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

"""
Core cloud diagnosis algorithm shared by DiagnosticClouds and SpectralDiagnosticClouds.
Returns (cloud_cover, cloud_top, stratocumulus_cover) tuple.
"""
function _diagnose_cloud_properties(
    column::ColumnVariables,
    rh_min, rh_max, q_min, precip_weight, precip_max,
    use_strat, stab_min, stab_max, cover_max, clfact
)
    (; nlayers, land_fraction, humid, sat_humid, dry_static_energy) = column

    # Precipitation contribution
    P = precip_weight * sqrt(min(precip_max, (86400 * (column.rain_rate_large_scale + column.rain_rate_convection) / 1000)))
    
    kprtop = column.cloud_top
    humidity_term_kcltop = zero(P)
    kcltop = nlayers + 1

    # Find cloud top from RH threshold
    for k in 1:(nlayers-1)
        humidity_k = humid[k]
        qsat = sat_humid[k]

        if humidity_k > q_min && qsat > 0
            relative_humidity_k = humidity_k / qsat
            
            if relative_humidity_k >= rh_min
                rh_norm = max(0, (relative_humidity_k - rh_min) / (rh_max - rh_min))
                humidity_term_kcltop = min(1, rh_norm)^2
                kcltop = k
                break
            end
        end
    end
    
    # Combined cloud cover
    cloud_cover = min(1, P + humidity_term_kcltop)
    cloud_top = min(kcltop, kprtop)
    column.cloud_top = cloud_top

    # Stratocumulus parameterization
    stratocumulus_cover = zero(P)
    if use_strat
        GSEN = dry_static_energy[nlayers] - dry_static_energy[nlayers-1]
        F_ST = max(0, min(1, (GSEN - stab_min) / (stab_max - stab_min)))
        
        stratocumulus_cover_ocean = F_ST * max(cover_max - clfact * cloud_cover, 0)
        RH_N = sat_humid[nlayers] > 0 ? humid[nlayers] / sat_humid[nlayers] : 0
        stratocumulus_cover_land = stratocumulus_cover_ocean * RH_N
        
        stratocumulus_cover = (1 - land_fraction) * stratocumulus_cover_ocean + land_fraction * stratocumulus_cover_land
    end

    return (cloud_cover, cloud_top, stratocumulus_cover)
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

    # Diagnose cloud cover and cloud top
    (cloud_cover, cloud_top, stratocumulus_cover) = _diagnose_cloud_properties(
        column,
        clouds.relative_humidity_threshold_min,
        clouds.relative_humidity_threshold_max,
        clouds.specific_humidity_threshold_min,
        clouds.precipitation_weight,    
        clouds.precipitation_max,
        clouds.use_stratocumulus,
        clouds.stratocumulus_stability_min,
        clouds.stratocumulus_stability_max,
        clouds.stratocumulus_cover_max,
        clouds.stratocumulus_clfact,
    )

    return (    # NamedTuple
        cloud_cover=cloud_cover,
        cloud_albedo=clouds.cloud_albedo,
        cloud_top=cloud_top,
        stratocumulus_cover=stratocumulus_cover,
        stratocumulus_albedo=clouds.stratocumulus_albedo,
    )
end

# Add clouds! method for regular DiagnosticClouds with band parameter (ignores band)
clouds!(column, clouds::DiagnosticClouds, model, band) = clouds!(column, clouds, model)


"""
Spectral cloud schemes with wavelength-dependent optical properties
"""

export SpectralDiagnosticClouds

"""
    SpectralDiagnosticClouds{N} <: AbstractShortwaveClouds

Diagnostic cloud scheme with band-dependent optical properties for N-band radiation.
Cloud cover and cloud top are diagnosed once but cloud albedos vary by spectral band
to represent realistic wavelength dependence of cloud droplet scattering.

$(TYPEDFIELDS)
"""
@kwdef struct SpectralDiagnosticClouds{N, NF} <: AbstractShortwaveClouds
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

    "[OPTION] Cloud albedo for each spectral band at CLC=1 [1]"
    cloud_albedo::NTuple{N, NF} = ntuple(i -> NF(0.43), N)

    "[OPTION] Stratocumulus cloud albedo for each spectral band [1]"
    stratocumulus_albedo::NTuple{N, NF} = ntuple(i -> NF(0.5), N)

    "[OPTION] Static stability lower threshold for stratocumulus [J/kg]"
    stratocumulus_stability_min::NF = 0.25

    "[OPTION] Static stability upper threshold for stratocumulus [J/kg]"
    stratocumulus_stability_max::NF = 0.40

    "[OPTION] Maximum stratocumulus cloud cover [1]"
    stratocumulus_cover_max::NF = 0.60

    "[OPTION] Enable stratocumulus cloud parameterization?"
    use_stratocumulus::Bool = true

    "[OPTION] Stratocumulus cloud factor [1]"
    stratocumulus_clfact::NF = 1.2
end

# Constructor with realistic spectral cloud albedos
function SpectralDiagnosticClouds{2}(SG::SpectralGrid; 
    # Realistic cloud albedos: higher in visible, lower in near-IR
    cloud_albedo = (SG.NF(0.50), SG.NF(0.35)),           # [visible, near-IR]
    stratocumulus_albedo = (SG.NF(0.55), SG.NF(0.40)),   # [visible, near-IR]
    kwargs...)
    
    return SpectralDiagnosticClouds{2, SG.NF}(;
        cloud_albedo=cloud_albedo,
        stratocumulus_albedo=stratocumulus_albedo,
        kwargs...
    )
end

SpectralDiagnosticClouds(SG::SpectralGrid; kwargs...) = SpectralDiagnosticClouds{1}(SG; kwargs...)
SpectralDiagnosticClouds{N}(SG::SpectralGrid; kwargs...) where N = SpectralDiagnosticClouds{N, SG.NF}(; kwargs...)

initialize!(clouds::SpectralDiagnosticClouds, model::AbstractModel) = nothing

"""$(TYPEDSIGNATURES)
Band-specific cloud calculation. Returns cloud properties for a specific spectral band."""
function clouds!(
    column::ColumnVariables,
    clouds::SpectralDiagnosticClouds{N},
    model::AbstractModel,
    band::Int = 1  # Which spectral band to compute
) where N

    # Diagnose cloud cover and cloud top once (shared across bands)
    (cloud_cover, cloud_top, stratocumulus_cover) = _diagnose_cloud_properties(
        column,
        clouds.relative_humidity_threshold_min,
        clouds.relative_humidity_threshold_max,
        clouds.specific_humidity_threshold_min,
        clouds.precipitation_weight,
        clouds.precipitation_max,
        clouds.use_stratocumulus,
        clouds.stratocumulus_stability_min,
        clouds.stratocumulus_stability_max,
        clouds.stratocumulus_cover_max,
        clouds.stratocumulus_clfact,
    )

    # Return band-specific albedos
    return (
        cloud_cover=cloud_cover,
        cloud_albedo=clouds.cloud_albedo[band],        # Band-specific!
        cloud_top=cloud_top,
        stratocumulus_cover=stratocumulus_cover,
        stratocumulus_albedo=clouds.stratocumulus_albedo[band],  # Band-specific!
    )
end