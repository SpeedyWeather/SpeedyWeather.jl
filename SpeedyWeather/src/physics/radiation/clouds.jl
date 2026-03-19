abstract type AbstractShortwaveClouds <: AbstractShortwave end

export NoClouds
struct NoClouds <: AbstractShortwaveClouds end
Adapt.@adapt_structure NoClouds
NoClouds(SG::SpectralGrid) = NoClouds()
initialize!(clouds::NoClouds, ::AbstractModel) = nothing
@propagate_inbounds function clouds!(
        ij,
        diagn,
        progn,
        ::NoClouds,
        model,
    )
    NF = eltype(diagn.grid.temp_grid_prev)
    cloud_cover = zero(NF)           # no cloud cover
    cloud_top = diagn.nlayers + 1    # below surface

    return (    # NamedTuple
        cloud_cover = cloud_cover,
        cloud_albedo = zero(NF),
        cloud_top = cloud_top,
        stratocumulus_cover = zero(NF),
        stratocumulus_albedo = zero(NF),
    )
end

export DiagnosticClouds
@parameterized @kwdef struct DiagnosticClouds{NF} <: AbstractShortwaveClouds
    "[OPTION] Relative humidity threshold for cloud cover = 0 [1]."
    @param relative_humidity_threshold_min::NF = 0.3 (bounds = 0 .. 1,)

    "[OPTION] Relative humidity threshold for cloud cover = 1 [1]"
    @param relative_humidity_threshold_max::NF = 1 (bounds = 0 .. 1,)

    "[OPTION] Specific humidity threshold for cloud cover [Kg/kg]"
    @param specific_humidity_threshold_min::NF = 0.0002 (bounds = 0 .. 1,)

    "[OPTION] Weight for √precip term [1]"
    @param precipitation_weight::NF = 0.2 (bounds = Nonnegative,)

    "[OPTION] Cap on precip contributing to cloud cover [mm/day]"
    @param precipitation_max::NF = 10 (bounds = Positive,)

    "[OPTION] Cloud albedo for visible band at CLC=1 [1]"
    @param cloud_albedo::NF = 0.6 (bounds = 0 .. 1,)

    "[OPTION] Stratocumulus cloud albedo (surface reflection/absorption) [1]"
    @param stratocumulus_albedo::NF = 0.5 (bounds = 0 .. 1,)

    "[OPTION] Static stability lower threshold for stratocumulus (GSEN, called GSES0 in the paper) [J/kg]"
    @param stratocumulus_stability_min::NF = 0.25 (bounds = Nonnegative,)

    "[OPTION] Static stability upper threshold for stratocumulus (GSEN, called GSES1 in the paper) [J/kg]"
    @param stratocumulus_stability_max::NF = 0.4 (bounds = Nonnegative,)

    "[OPTION] Maximum stratocumulus cloud cover (called CLSMAX in the paper) [1]"
    @param stratocumulus_cover_max::NF = 0.6 (bounds = 0 .. 1,)

    "[OPTION] Enable stratocumulus cloud parameterization?"
    use_stratocumulus::Bool = true

    "[OPTION] Stratocumulus cloud factor (SPEEDY clfact) [1]"
    @param stratocumulus_cloud_factor::NF = 1.2 (bounds = Nonnegative,)
end

Adapt.@adapt_structure DiagnosticClouds
DiagnosticClouds(SG::SpectralGrid; kwargs...) = DiagnosticClouds{SG.NF}(; kwargs...)
initialize!(clouds::DiagnosticClouds, model::AbstractModel) = nothing
@propagate_inbounds function clouds!(
        ij,
        diagn,
        progn,
        clouds::DiagnosticClouds,
        model,
    )
    # Diagnose cloud cover and cloud top.
    (; cloud_cover, cloud_top, stratocumulus_cover) = diagnose_cloud_properties(ij, diagn, clouds, model)

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
Returns (cloud_cover, cloud_top, stratocumulus_cover) tuple."""
@propagate_inbounds function diagnose_cloud_properties(ij, diagn, clouds::DiagnosticClouds, model)

    temp = diagn.grid.temp_grid_prev
    humid = diagn.grid.humid_grid_prev
    geopotential = diagn.grid.geopotential
    NF = eltype(temp)
    nlayers = size(temp, 2)

    pₛ = diagn.grid.pres_grid_prev[ij]
    sigma_levels = model.geometry.σ_levels_full
    land_fraction = model.land_sea_mask.mask[ij]
    cₚ = model.atmosphere.heat_capacity

    # rename for brevity
    rh_min = clouds.relative_humidity_threshold_min
    rh_max = clouds.relative_humidity_threshold_max
    q_min = clouds.specific_humidity_threshold_min
    precip_weight = clouds.precipitation_weight
    precip_max = clouds.precipitation_max
    stab_min = clouds.stratocumulus_stability_min
    stab_max = clouds.stratocumulus_stability_max
    cover_max = clouds.stratocumulus_cover_max
    cloud_factor = clouds.stratocumulus_cloud_factor

    # Precipitation contribution (rain rate is in m/s)
    rain_rate = diagn.physics.rain_rate[ij]
    precip_term = min(precip_max, (86400 * rain_rate) / 1000)   # convert to mm/day
    P = precip_weight * sqrt(precip_term)

    cloud_top_precipitation = diagn.physics.cloud_top[ij]       # from convection or large-scale condensation
    humidity_term_cloud_top::NF = 0
    cloud_top_humidity = nlayers + 1

    # Find cloud top from RH threshold
    for k in 1:(nlayers - 1)
        humidity_k = humid[ij, k]
        qsat = saturation_humidity(temp[ij, k], sigma_levels[k] * pₛ, model.atmosphere)
        if humidity_k > q_min && qsat > 0
            relative_humidity_k = humidity_k / qsat

            if relative_humidity_k >= rh_min
                rh_norm = max(0, (relative_humidity_k - rh_min) / (rh_max - rh_min))
                humidity_term_cloud_top = min(1, rh_norm)^2
                cloud_top_humidity = min(k, cloud_top_humidity)
            end
        end
    end

    # Combined cloud cover
    cloud_cover = min(1, P + humidity_term_cloud_top)
    cloud_top = min(cloud_top_humidity, cloud_top_precipitation)
    diagn.physics.cloud_top[ij] = cloud_top

    # Stratocumulus parameterization
    stratocumulus_cover::NF = 0         # fallback for no stratocumulus
    if clouds.use_stratocumulus
        # Vertical difference of dry static energy near surface
        surface = nlayers               # surface index
        above = max(1, nlayers - 1)     # layer above surface, max for 1 layer case, yields GSEN=0
        G = (cₚ * temp[ij, surface] + geopotential[ij, surface]) -
            (cₚ * temp[ij, above] + geopotential[ij, above])

        # normalized static stability factor
        static_stability = clamp((G - stab_min) / (stab_max - stab_min), 0, 1)

        stratocumulus_cover_ocean = static_stability * max(cover_max - cloud_factor * cloud_cover, 0)
        qsat_surface = saturation_humidity(temp[ij, surface], sigma_levels[surface] * pₛ, model.atmosphere)
        rh_surface = humid[ij, surface] / qsat_surface      # relative humidity at surface
        stratocumulus_cover_land = stratocumulus_cover_ocean * rh_surface

        stratocumulus_cover =
            (1 - land_fraction) * stratocumulus_cover_ocean +
            land_fraction * stratocumulus_cover_land
    end

    return (; cloud_cover, cloud_top, stratocumulus_cover)
end
