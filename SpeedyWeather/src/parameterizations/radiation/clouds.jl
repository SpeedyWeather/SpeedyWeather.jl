abstract type AbstractShortwaveClouds <: AbstractShortwave end

export NoClouds
struct NoClouds <: AbstractShortwaveClouds end
Adapt.@adapt_structure NoClouds
NoClouds(SG::SpectralGrid) = NoClouds()
initialize!(clouds::NoClouds, ::AbstractModel) = nothing
@propagate_inbounds function clouds!(
        ij,
        vars,
        clouds::NoClouds,
        model,
    )
    temp = get_prognostic_step(vars.grid.temperature, model.time_stepping, clouds)
    NF = eltype(temp)
    cloud_cover = zero(NF)           # no cloud cover
    nlayers = size(temp, 2)
    cloud_top = nlayers + 1          # below surface

    return (    # NamedTuple
        cloud_cover = cloud_cover,
        cloud_albedo = zero(NF),
        cloud_top = cloud_top,
        cloud_base = nlayers,
        stratocumulus_cover = zero(NF),
        stratocumulus_albedo = zero(NF),
    )
end

export DiagnosticClouds
@parameterized @kwdef struct DiagnosticClouds{NF, RI} <: AbstractShortwaveClouds
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
    @param cloud_albedo::NF = 0.43 (bounds = 0 .. 1,)

    "[OPTION] Stratocumulus cloud albedo (surface reflection/absorption) [1]"
    @param stratocumulus_albedo::NF = 0.5 (bounds = 0 .. 1,)

    "[OPTION] Static stability lower threshold for stratocumulus, as dry static energy gradient Δs/ΔΦ [1]"
    @param stratocumulus_stability_min::NF = 0.25 (bounds = Nonnegative,)

    "[OPTION] Static stability upper threshold for stratocumulus, as dry static energy gradient Δs/ΔΦ [1]"
    @param stratocumulus_stability_max::NF = 0.4 (bounds = Nonnegative,)

    "[OPTION] Maximum stratocumulus cloud cover (called CLSMAX in the paper) [1]"
    @param stratocumulus_cover_max::NF = 0.6 (bounds = 0 .. 1,)

    "[OPTION] Minimum stratocumulus cloud cover over land (for surface RH = 1) [1]"
    @param stratocumulus_cover_min_land::NF = 0.15 (bounds = 0 .. 1,)

    "[OPTION] Clouds only form in layers below the tropopause (σ ≥ σ_tropopause) [1]"
    σ_tropopause::NF = 0.14

    "[OPTION] Top of the boundary layer: clouds form in layers above (σ ≤ σ_boundary_layer), stratocumulus at its top [1]"
    σ_boundary_layer::NF = 0.9

    "[OPTION] Enable stratocumulus cloud parameterization?"
    use_stratocumulus::Bool = true

    "[OPTION] Stratocumulus cloud factor [1]"
    @param stratocumulus_cloud_factor::NF = 1.2 (bounds = Nonnegative,)

    "[OPTION] σ thickness over which relative humidity is averaged (Δσ-weighted, centred on each layer)
    before taking its maximum over the free troposphere, so that cloud cover does not depend on the
    number of layers, e.g. 1/8 ≈ a layer at L8; 0 = the layer's own relative humidity (legacy)"
    @param humidity_smoothing_thickness::NF = 0 (bounds = Nonnegative,)

    "[DERIVED] Highest layer in which clouds can form, first layer below the tropopause, set in `initialize!`"
    layer_top::RI = Ref(0)

    "[DERIVED] Cloud base, lowest layer above the boundary layer, set in `initialize!`"
    cloud_base::RI = Ref(0)
end

Adapt.@adapt_structure DiagnosticClouds
DiagnosticClouds(SG::SpectralGrid; kwargs...) = DiagnosticClouds{SG.NF, Base.RefValue{Int}}(; kwargs...)

# the free troposphere layers only depend on the vertical coordinates, precompute them
function initialize!(clouds::DiagnosticClouds, model::AbstractModel)
    clouds.layer_top[], clouds.cloud_base[] = free_troposphere_layers(clouds, model.geometry.vertical_coordinates)
    return nothing
end

variables(::DiagnosticClouds) = (
    ParameterizationVariable(:cloud_top, Grid2D(), desc = "Cloud top layer index", units = "1"),
    ParameterizationVariable(:cloud_top_height, Grid2D(), desc = "Cloud top height", units = "m"),
    ParameterizationVariable(:cloud_cover, Grid2D(), desc = "Cloud cover", units = "1"),
)

@propagate_inbounds function clouds!(
        ij,
        vars,
        clouds::DiagnosticClouds,
        model,
    )
    # Diagnose cloud cover and cloud top.
    (; cloud_cover, cloud_top, cloud_base, stratocumulus_cover) = diagnose_cloud_properties(ij, vars, clouds, model)

    return (    # NamedTuple
        cloud_cover = cloud_cover,
        cloud_albedo = clouds.cloud_albedo,
        cloud_top = cloud_top,
        cloud_base = cloud_base,
        stratocumulus_cover = stratocumulus_cover,
        stratocumulus_albedo = clouds.stratocumulus_albedo,
    )
end

"""$(TYPEDSIGNATURES)
Core cloud diagnosis algorithm shared by DiagnosticClouds and SpectralDiagnosticClouds.
Returns (cloud_cover, cloud_top, cloud_base, stratocumulus_cover) tuple."""
@propagate_inbounds function diagnose_cloud_properties(ij, vars, clouds::DiagnosticClouds, model)

    temp = get_prognostic_step(vars.grid.temperature, model.time_stepping, clouds)
    humid = get_prognostic_step(vars.grid.humidity, model.time_stepping, clouds)
    geopotential = vars.dynamics.geopotential
    NF = eltype(temp)
    nlayers = size(temp, 2)

    pₛ = vars.parameterizations.surface_pressure[ij]          # surface pressure [Pa]
    (; vertical_coordinates) = model.geometry
    land_fraction = model.land_sea_mask.land_fraction[ij]
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
    cover_min_land = clouds.stratocumulus_cover_min_land

    # Precipitation contribution from rain and snow (both in m/s)
    (; parameterizations) = vars
    snow_rate = haskey(parameterizations, :snow_rate) ? parameterizations.snow_rate[ij] : zero(NF)
    precip_rate = parameterizations.rain_rate[ij] + snow_rate
    precip_term = min(precip_max, max(0, precip_rate) * 86400 * 1000)   # convert m/s to mm/day
    P = precip_weight * sqrt(precip_term)

    # from convection or large-scale condensation
    cloud_top_precipitation = vars.parameterizations.cloud_top[ij]

    # clouds form in the free troposphere, layers between the tropopause and the boundary layer
    layer_top, cloud_base = clouds.layer_top[], clouds.cloud_base[]

    # Find the layer of maximum relative humidity in the free troposphere (with q > q_min,
    # except for the cloud base layer directly above the boundary layer as in SPEEDY, so that
    # clouds can also form in cold and dry polar air), cloud cover is a quadratic function of
    # that maximum, the cloud top is its layer
    rh_excess_max::NF = 0
    cloud_top_humidity = nlayers + 1
    σ = model.geometry.σ_levels_full
    Δσ = model.geometry.σ_levels_thick
    half_window = clouds.humidity_smoothing_thickness / 2
    smooth = half_window > 0

    # with smoothing, only layers whose averaging window lies fully above the boundary layer are
    # candidates (at high vertical resolution the cloud_base layer sits right at the humid boundary-layer
    # top, at L8 it is centred well above); the q_min exemption moves to the lowest such layer
    σ_bl = clouds.σ_boundary_layer
    k_low = cloud_base
    if smooth
        k_low = layer_top
        for k in layer_top:cloud_base
            k_low = ifelse(σ[k] + half_window <= σ_bl, k, k_low)
        end
    end

    for k in layer_top:cloud_base
        humidity_k = humid[ij, k]
        pₖ = pressure(k, pₛ, vertical_coordinates)
        qsat = saturation_humidity(temp[ij, k], pₖ, model.atmosphere)
        eligible = !smooth | (k <= k_low)
        valid = ((humidity_k > q_min) | (k == k_low)) & (qsat > 0) & eligible
        rh_k = humidity_k / ifelse(qsat > 0, qsat, one(NF))         # avoid division by zero

        # optionally: Δσ-weighted mean relative humidity of the free-tropospheric layers within
        # ±half_window of layer k (always including k itself), independent of the layer count
        if smooth
            rh_sum = zero(NF)
            w_sum = zero(NF)
            for j in layer_top:cloud_base
                pⱼ = pressure(j, pₛ, vertical_coordinates)
                qsatⱼ = saturation_humidity(temp[ij, j], pⱼ, model.atmosphere)
                w = ifelse((abs(σ[j] - σ[k]) <= half_window) | (j == k), Δσ[j], zero(NF))
                rh_sum += w * humid[ij, j] / ifelse(qsatⱼ > 0, qsatⱼ, one(NF))
                w_sum += w
            end
            rh_k = rh_sum / w_sum                                   # w_sum > 0 as j == k is included
        end
        rh_excess = rh_k - rh_min
        new_max = valid & (rh_excess > rh_excess_max)
        rh_excess_max = ifelse(new_max, rh_excess, rh_excess_max)
        cloud_top_humidity = ifelse(new_max, k, cloud_top_humidity)
    end
    humidity_term = min(1, rh_excess_max / (rh_max - rh_min))^2

    # Combined cloud cover
    cloud_cover = min(1, P + humidity_term)
    # highest of humidity and precipitation cloud top, but kept within the free troposphere
    cloud_top = min(cloud_top_humidity, cloud_top_precipitation)
    cloud_top = ifelse(cloud_top <= nlayers, clamp(cloud_top, layer_top, cloud_base), cloud_top)
    vars.parameterizations.cloud_top[ij] = cloud_top

    # store cloud cover and cloud top height [m] (0 for no cloud) for output
    vars.parameterizations.cloud_cover[ij] = cloud_cover
    # cloud_top is stored as float (integer-valued) so convert to index, nlayers + 1 = no cloud
    k_top = min(unsafe_trunc(Int, cloud_top), nlayers)
    Φ_top = geopotential[ij, k_top]
    vars.parameterizations.cloud_top_height[ij] = ifelse(cloud_top <= nlayers, Φ_top / model.planet.gravity, zero(NF))

    # Stratocumulus parameterization, computed always and switched off with ifelse (branch-free)
    # Static stability as vertical gradient of dry static energy s = cₚT + Φ between the
    # surface layer and the cloud base layer above the boundary layer, normalized by the
    # geopotential difference, i.e. GSE = Δs/ΔΦ = 1 + cₚΔT/ΔΦ, 0 for dry adiabatic,
    # 1 for isothermal. With 8 layers (as in SPEEDY) these are the lowest two layers
    surface = nlayers
    above = cloud_base
    ΔΦ = geopotential[ij, above] - geopotential[ij, surface]
    Δs = cₚ * (temp[ij, above] - temp[ij, surface]) + ΔΦ
    # Guard against division by zero (ΔΦ = 0 for a single layer, as above == surface) by
    # always dividing with guarded numerator and denominator. Do not write this as
    # ifelse(c, Δs / ifelse(c, ΔΦ, 1), 0): LLVM removes the inner ifelse as c is known to be true
    # in the taken branch, so 0/0 = NaN is computed and propagates into the reverse-mode gradients
    stable_layers = ΔΦ > 0
    GSE = ifelse(stable_layers, Δs, zero(NF)) / ifelse(stable_layers, ΔΦ, one(NF))

    # normalized static stability factor
    static_stability = clamp((GSE - stab_min) / (stab_max - stab_min), 0, 1)

    stratocumulus_cover_ocean = static_stability * max(cover_max - cloud_factor * cloud_cover, 0)
    pN = pressure(surface, pₛ, vertical_coordinates)
    qsat_surface = saturation_humidity(temp[ij, surface], pN, model.atmosphere)
    rh_surface = humid[ij, surface] / ifelse(qsat_surface > 0, qsat_surface, one(NF))  # relative humidity at surface
    stratocumulus_cover_land = max(stratocumulus_cover_ocean, cover_min_land) * rh_surface

    stratocumulus_cover_land_sea =
        (1 - land_fraction) * stratocumulus_cover_ocean +
        land_fraction * stratocumulus_cover_land
    stratocumulus_cover = ifelse(clouds.use_stratocumulus & (cloud_base < nlayers), stratocumulus_cover_land_sea, zero(NF))

    return (; cloud_cover, cloud_top, cloud_base, stratocumulus_cover)
end

"""$(TYPEDSIGNATURES)
Range of layers `layer_top:cloud_base` of the free troposphere in which clouds can form:
from the first layer below the tropopause (σ ≥ `σ_tropopause`) to the cloud base, the lowest
layer above the boundary layer (σ ≤ `σ_boundary_layer`). Clouds never include the surface layer
(unless there is only one layer), the cloud base is at least the first layer below the tropopause
and at most the layer above the surface layer. For 8 equally spaced layers these are layers 2 to 7.
Uses the surface pressure-independent `sigma(k, vertical_coordinates)`, so this is computed once
in `initialize!` for any vertical coordinates."""
function free_troposphere_layers(clouds::DiagnosticClouds, vertical_coordinates::AbstractVerticalCoordinates)
    nlayers = get_nlayers(vertical_coordinates)
    # σ increases monotonically downwards, so layer indices are counts of layers above
    # a threshold, vectorized on the architecture of the vertical coordinates (no scalar indexing)
    σ = sigma_full(vertical_coordinates)
    σ_above_surface = view(σ, 1:(nlayers - 1))                  # exclude the surface layer
    layer_top = count(<(clouds.σ_tropopause), σ_above_surface) + 1
    cloud_base = count(<=(clouds.σ_boundary_layer), σ_above_surface)
    cloud_base = max(cloud_base, min(layer_top, nlayers - 1), 1)
    layer_top = min(layer_top, cloud_base)
    return layer_top, cloud_base
end
