abstract type AbstractShortwaveClouds <: AbstractShortwave end

export NoClouds
struct NoClouds <: AbstractShortwaveClouds end
NoClouds(SG::SpectralGrid) = NoClouds()
initialize!(clouds::NoClouds, ::AbstractModel) = nothing
function clouds!(
        column::ColumnVariables{NF},
        ::NoClouds,
        model::AbstractModel,
    ) where {NF}
    cloud_cover = zero(NF)          # no cloud cover
    cloud_top = column.nlayers + 1    # below surface

    return (    # NamedTuple
        cloud_cover = cloud_cover,
        cloud_albedo = zero(NF),
        cloud_top = cloud_top,
        stratocumulus_cover = zero(NF),
        stratocumulus_albedo = zero(NF),
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
    stratocumulus_stability_max::NF = 0.4

    "[OPTION] Maximum stratocumulus cloud cover (called CLSMAX in the paper) [1]"
    stratocumulus_cover_max::NF = 0.6

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
    # Diagnose cloud cover and cloud top.
    (cloud_cover, cloud_top, stratocumulus_cover) = diagnose_cloud_properties(column, clouds)

    return (    # NamedTuple
        cloud_cover = cloud_cover,
        cloud_albedo = clouds.cloud_albedo,
        cloud_top = cloud_top,
        stratocumulus_cover = stratocumulus_cover,
        stratocumulus_albedo = clouds.stratocumulus_albedo,
    )
end

"""$(TYPEDSIGNATURES)
Core cloud diagnosis algorithm shared by DiagnosticClouds and SpectralDiagnosticClouds.
Returns (`cloud_cover`, `cloud_top`, `stratocumulus_cover`) tuple."""
function diagnose_cloud_properties(column::ColumnVariables, clouds::DiagnosticClouds)

    # rename for brevity
    rh_min = clouds.relative_humidity_threshold_min
    rh_max = clouds.relative_humidity_threshold_max
    q_min = clouds.specific_humidity_threshold_min
    precip_weight = clouds.precipitation_weight
    precip_max = clouds.precipitation_max
    use_strat = clouds.use_stratocumulus
    stab_min = clouds.stratocumulus_stability_min
    stab_max = clouds.stratocumulus_stability_max
    cover_max = clouds.stratocumulus_cover_max
    clfact = clouds.stratocumulus_clfact

    (; nlayers, land_fraction, humid, sat_humid, dry_static_energy) = column

    # Precipitation contribution
    P = precip_weight * sqrt(min(precip_max, (86400 * (column.rain_rate_large_scale + column.rain_rate_convection) / 1000)))

    kprtop = column.cloud_top
    humidity_term_kcltop = zero(P)
    kcltop = nlayers + 1

    # Find cloud top from RH threshold
    for k in 1:(nlayers - 1)
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
        GSEN = dry_static_energy[nlayers] - dry_static_energy[nlayers - 1]
        F_ST = max(0, min(1, (GSEN - stab_min) / (stab_max - stab_min)))

        stratocumulus_cover_ocean = F_ST * max(cover_max - clfact * cloud_cover, 0)
        RH_N = sat_humid[nlayers] > 0 ? humid[nlayers] / sat_humid[nlayers] : 0
        stratocumulus_cover_land = stratocumulus_cover_ocean * RH_N

        stratocumulus_cover = (1 - land_fraction) * stratocumulus_cover_ocean + land_fraction * stratocumulus_cover_land
    end

    return (cloud_cover, cloud_top, stratocumulus_cover)
end
