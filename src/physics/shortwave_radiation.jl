abstract type AbstractShortwave <: AbstractRadiation end

function get_nbands(radiation::Union{AbstractRadiation, Nothing})
    hasfield(typeof(radiation), :nbands) && return radiation.nbands
    return 0
end

# function barrier for all AbstractShortwave
function shortwave_radiation!(column::ColumnVariables, model::PrimitiveEquation)
    shortwave_radiation!(column, model.shortwave_radiation, model)
end

## NO SHORTWAVE RADIATION
shortwave_radiation!(::ColumnVariables, ::Nothing, ::PrimitiveEquation) = nothing

## SHORTWAVE RADIATION FOR A FULLY TRANSPARENT ATMOSPHERE
export TransparentShortwave
struct TransparentShortwave <: AbstractShortwave end
TransparentShortwave(SG::SpectralGrid) = TransparentShortwave()
initialize!(::TransparentShortwave, ::PrimitiveEquation) = nothing

function shortwave_radiation!(
    column::ColumnVariables,
    scheme::TransparentShortwave,
    model::PrimitiveEquation,
)
    shortwave_radiation!(column, scheme, model.planet)
end

function shortwave_radiation!(
    column::ColumnVariables,
    scheme::TransparentShortwave,
    planet::AbstractPlanet,
)
    (; cos_zenith, land_fraction, albedo_ocean, albedo_land) = column
    (; solar_constant) = planet

    R = solar_constant * cos_zenith         # top of atmosphere downward radiation
    column.surface_shortwave_down = R       # transparent atmosphere so same at surface (before albedo)

    # shortwave up is after albedo reflection, separated by ocean/land
    column.surface_shortwave_up_ocean = albedo_ocean * R
    column.surface_shortwave_up_land = albedo_land * R

    # land-sea mask-weighted, transparent shortwave so surface = outgoing
    column.surface_shortwave_up = (1 - land_fraction)*column.surface_shortwave_up_ocean +
                                            land_fraction*column.surface_shortwave_up_land
    column.outgoing_shortwave_radiation = column.surface_shortwave_up

    return nothing
end

# NBandRadiation is defined in longwave_radiation.jl
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


############### SPEEDY B.4 — One-Band Shortwave with Diagnostic Clouds ###############

export OneBandShortwave

struct OneBandShortwave <: AbstractShortwave
    # thresholds / weights
    rh_cl::Float64        # RH_cl,  thresholds for RH and cloud water content (Eqs. 26–27).
    rh_clp::Float64       # RH'_cl,     ---------------------''-----------------------
    q_cl::Float64         # Q_cl,       ---------------------''-----------------------
    wp_cl::Float64        # w_pcl   weighting of precipitation term (Eq. 28).
    pmax_cl::Float64      # p_maxcl (mm/day) 
    # optics
    alb_cl::Float64       # A_cl (cloud albedo, visible) (Eq. 33).
end

"""
    OneBandShortwave(; 
        rh_cl=0.30,          # RH threshold for cloud cover = 0 (–)
        rh_clp=1.00,         # RH threshold for cloud cover = 1 (–)
        q_cl=0.20,           # Absolute humidity threshold for cloud cover (g kg⁻¹)
        wp_cl=0.20,          # Weight for √precip term (–) at 1 mm day⁻¹
        pmax_cl=10.0,        # Cap on precip contributing to cloud cover (mm day⁻¹)
        alb_cl=0.43,         # Cloud albedo for visible band at CLC=1 (–)
    )

Defaults mirror SPEEDY’s shortwave constants. 
Notes:
- `abs_wv_vis` expects q in g kg⁻¹. If your model uses kg kg⁻¹, internally scale by 1000 (or rescale the coefficient).
- `abs_dry` is calibrated per layer mass of Δp = 10⁵ Pa (SPEEDY convention). Match your Beer–Lambert exponent accordingly.
- If you don’t implement Eq. 34 stratocumulus, you can ignore `alb_cls`, `cls_*`, `clf`, `gses*`.
"""
OneBandShortwave() = OneBandShortwave(
    0.30, 1.00, 0.20,   # rh_cl (–), rh_clp (–), q_cl (g kg⁻¹)
    0.20, 10.0,         # wp_cl (–), pmax_cl (mm day⁻¹)
    0.43                # alb_cl (–)
)

get_nbands(::OneBandShortwave) = 1
initialize!(::OneBandShortwave, ::PrimitiveEquation) = nothing


function shortwave_radiation!(
    column::ColumnVariables,
    scheme::OneBandShortwave,
    model::PrimitiveEquation,
)
    # ---- pull only fields that exist on ColumnVariables ----
    (; nlayers, cos_zenith, land_fraction,
       humid, sat_humid, albedo_ocean, albedo_land,
       optical_depth_shortwave, flux_temp_downward, flux_temp_upward,
       cloud_top, rain_rate_convection, rain_rate_large_scale) = column

    # ---- aliases / guards ----
    μ   = max(cos_zenith, 1f-3)                         # avoid divide-by-zero
    qv  = humid
    qsat = sat_humid
    T   = eltype(qv)
    τk = @view optical_depth_shortwave[1:nlayers, 1]

    # ---- derived diagnostics ----
    RH = @. clamp(qv / max(qsat, eps(T)), 0, 1) # clamp relative humidity between 0 and 1
    k_cl = (1 <= cloud_top <= nlayers) ? cloud_top : max(nlayers - 1, 1) # cloud top layer index
    RH_cl  = scheme.rh_cl 
    RH_clp = scheme.rh_clp
    r      = clamp((RH[min(k_cl, nlayers)] - RH_cl) / max(RH_clp - RH_cl, T(1e-6)), zero(T), one(T))

    P_mmday = T(86.4) * (rain_rate_convection + rain_rate_large_scale)
    P_mmday = min(scheme.pmax_cl, P_mmday)
    CLC     = min(one(T), scheme.wp_cl * sqrt(P_mmday) + r^2)

    # ---- transmissivity per layer: τ_layer = exp(- τk / μ) ----
    τ = similar(qv)
    @inbounds @simd for k in 1:nlayers
        τ[k] = exp(-τk[k] / μ)
    end

    # ---- 1) Downward beam ----
    S0 = model.planet.solar_constant
    D  = T(S0) * μ
    @inbounds flux_temp_downward[1] += D

    D_incident_ct = zero(T)
    # Clear sky portion until cloud top
    @inbounds for k in 1:(k_cl - 1)
        D *= τ[k]
        flux_temp_downward[k + 1] += D
    end
    
    # Cloudy portion from cloud top downward
    D_incident_ct = D
    D *= (one(T) - scheme.alb_cl * CLC)
    @inbounds for k in k_cl:nlayers
        D *= τ[k]
        flux_temp_downward[k + 1] += D
    end

    # ---- 2) Surface reflection ----
    column.surface_shortwave_down = D
    up_ocean = albedo_ocean * D
    up_land  = albedo_land  * D
    column.surface_shortwave_up_ocean = up_ocean
    column.surface_shortwave_up_land  = up_land

    # Computes grid-cell-average surface albedo and reflected shortwave flux, 
    albedo_sfc = (one(T) - land_fraction) * albedo_ocean + land_fraction * albedo_land
    column.surface_shortwave_up = (one(T) - land_fraction) * up_ocean + land_fraction * up_land

    U = albedo_sfc * D
    @inbounds flux_temp_upward[nlayers + 1] += U

    # ---- 3) Upward beam (+ extra cloud-top reflection) ----
    @inbounds for k in nlayers:-1:1
        U *= τ[k]
        flux_temp_upward[k] += U
        if k == k_cl
            dU = scheme.alb_cl * CLC * D_incident_ct
            U  += dU
            flux_temp_upward[k] += dU
        end
    end

    # ---- 4) TOA outgoing SW ----
    column.outgoing_shortwave_radiation += U

    return nothing
end