abstract type AbstractRadiation <: AbstractParameterization end
abstract type AbstractShortwave <: AbstractRadiation end
abstract type AbstractShortwaveRadiativeTransfer <: AbstractShortwave end

export TransparentShortwave
struct TransparentShortwave <: AbstractShortwave end
Adapt.@adapt_structure TransparentShortwave
TransparentShortwave(SG::SpectralGrid) = TransparentShortwave()

function variables(::AbstractShortwave)
    return (
        ParameterizationVariable(:surface_shortwave_down, Grid2D(), desc = "Surface shortwave radiation down", units = "W/m^2"),
        ParameterizationVariable(:surface_shortwave_down, Grid2D(), desc = "Surface shortwave radiation down over ocean", units = "W/m^2", namespace = :ocean),
        ParameterizationVariable(:surface_shortwave_down, Grid2D(), desc = "Surface shortwave radiation down over land", units = "W/m^2", namespace = :land),
        ParameterizationVariable(:surface_shortwave_up, Grid2D(), desc = "Surface shortwave radiation up", units = "W/m^2"),
        ParameterizationVariable(:surface_shortwave_up, Grid2D(), desc = "Surface shortwave radiation up over ocean", units = "W/m^2", namespace = :ocean),
        ParameterizationVariable(:surface_shortwave_up, Grid2D(), desc = "Surface shortwave radiation up over land", units = "W/m^2", namespace = :land),
        ParameterizationVariable(:outgoing_shortwave, Grid2D(), desc = "TOA Shortwave radiation up", units = "W/m^2"),
        ParameterizationVariable(:cos_zenith, Grid2D(), desc = "Cos zenith angle", units = "1"),
        ParameterizationVariable(:albedo, Grid2D(), desc = "Albedo", units = "1"),
        ParameterizationVariable(:albedo, Grid2D(), desc = "Albedo over ocean", units = "1", namespace = :ocean),
        ParameterizationVariable(:albedo, Grid2D(), desc = "Albedo over land", units = "1", namespace = :land),
    )
end

initialize!(::TransparentShortwave, ::PrimitiveEquation) = nothing

# function barrier
@propagate_inbounds function parameterization!(ij, vars, ::TransparentShortwave, model)

    planet = model.planet
    land_sea_mask = model.land_sea_mask.land_fraction

    cos_zenith = vars.parameterizations.cos_zenith[ij]
    land_fraction = land_sea_mask[ij]
    albedo_ocean = vars.parameterizations.ocean.albedo[ij]
    albedo_land = vars.parameterizations.land.albedo[ij]

    S₀ = planet.solar_constant
    D = S₀ * cos_zenith             # top of atmosphere downward radiation
    vars.parameterizations.surface_shortwave_down[ij] = D  # transparent atmosphere so same at surface (before albedo)
    vars.parameterizations.ocean.surface_shortwave_down[ij] = D
    vars.parameterizations.land.surface_shortwave_down[ij] = D

    # shortwave up is after albedo reflection, separated by ocean/land
    vars.parameterizations.ocean.surface_shortwave_up[ij] = albedo_ocean * D
    vars.parameterizations.land.surface_shortwave_up[ij] = albedo_land * D

    # land-sea mask-weighted
    albedo = (1 - land_fraction) * albedo_ocean + land_fraction * albedo_land
    vars.parameterizations.surface_shortwave_up[ij] = albedo * D
    vars.parameterizations.albedo[ij] = albedo   # store weighted albedo

    # transparent also for reflected shortwave radiation travelling up
    vars.parameterizations.outgoing_shortwave[ij] = vars.parameterizations.surface_shortwave_up[ij]
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

$(TYPEDFIELDS)
"""
@parameterized @kwdef struct OneBandShortwave{C, T, R} <: AbstractShortwave
    @component clouds::C
    @component transmissivity::T
    @component radiative_transfer::R
end
Adapt.@adapt_structure OneBandShortwave

# primitive wet model version
function OneBandShortwave(
        SG::SpectralGrid;
        clouds = DiagnosticClouds(SG),
        transmissivity = BackgroundShortwaveTransmissivity(SG),
        radiative_transfer = OneBandShortwaveRadiativeTransfer(SG),
    )
    return OneBandShortwave(clouds, transmissivity, radiative_transfer)
end

# primitive dry model version
export OneBandGreyShortwave
function OneBandGreyShortwave(
        SG::SpectralGrid;
        clouds = NoClouds(SG),
        transmissivity = ConstantShortwaveTransmissivity(SG),
        radiative_transfer = OneBandShortwaveRadiativeTransfer(SG),
    )
    return OneBandShortwave(clouds, transmissivity, radiative_transfer)
end

Base.show(io::IO, M::OneBandShortwave) = show(io, M, values = false)

# initialize one after another
function initialize!(radiation::OneBandShortwave, model::PrimitiveEquation)
    initialize!(radiation.clouds, model)
    initialize!(radiation.transmissivity, model)
    initialize!(radiation.radiative_transfer, model)
    return nothing
end

"""$(TYPEDSIGNATURES)
Calculate shortwave radiation using the one-band scheme with diagnostic clouds.
Computes cloud cover fraction from relative humidity and precipitation, then
integrates downward and upward radiative fluxes accounting for cloud albedo effects."""
@propagate_inbounds function parameterization!(
        ij,
        vars,
        radiation::OneBandShortwave,
        model,
    )
    clouds = clouds!(ij, vars, radiation.clouds, model)
    t = transmissivity!(ij, vars, clouds, radiation.transmissivity, model)
    shortwave_radiative_transfer!(ij, vars, t, clouds, radiation.radiative_transfer, model)
    return nothing
end

export OneBandShortwaveRadiativeTransfer
"""
    OneBandShortwaveRadiativeTransfer <: AbstractShortwaveRadiativeTransfer

$(TYPEDFIELDS)."""
@parameterized @kwdef struct OneBandShortwaveRadiativeTransfer{NF, F} <: AbstractShortwaveRadiativeTransfer
    "[OPTION] Total ozone absorption as fraction of incoming solar radiation (1)"
    @param ozone_absorption::NF = 0.01 (bounds = 0 .. 1,)

    "[OPTION] Ozone distribution above σ₀, has to be explicitly normalized to ∫dσ = 1 (1)"
    ozone_distribution::F
end
Adapt.@adapt_structure OneBandShortwaveRadiativeTransfer

# generator function
function OneBandShortwaveRadiativeTransfer(
        SG::SpectralGrid;
        ozone_distribution = (σ) -> 50 * max(0, 1 // 5 - σ),     # default distribution here
        kwargs...
    )
    return OneBandShortwaveRadiativeTransfer{SG.NF, typeof(ozone_distribution)}(;
        ozone_distribution = ozone_distribution, kwargs...
    )
end

initialize!(::OneBandShortwaveRadiativeTransfer, ::PrimitiveEquation) = nothing

"""$(TYPEDSIGNATURES)
One-band shortwave radiative transfer with cloud reflection and ozone absorption."""
@propagate_inbounds function shortwave_radiative_transfer!(
        ij,
        vars,
        t,          # Transmissivity array
        clouds,     # NamedTuple from clouds!
        radiation::OneBandShortwaveRadiativeTransfer,
        model,
    )

    O₃_absorption = radiation.ozone_absorption
    (; cloud_cover, cloud_top, stratocumulus_cover, cloud_albedo, stratocumulus_albedo) = clouds

    dTdt = get_tendency_step(vars.tendencies.grid.temperature, model.time_stepping, radiation)
    pₛ = vars.parameterizations.surface_pressure[ij]
    nlayers = size(dTdt, 2)
    σ = model.geometry.σ_levels_full
    Δσ = model.geometry.σ_levels_thick

    cos_zenith = vars.parameterizations.cos_zenith[ij]
    albedo_ocean = vars.parameterizations.ocean.albedo[ij]
    albedo_land = vars.parameterizations.land.albedo[ij]
    land_fraction = model.land_sea_mask.land_fraction[ij]
    cₚ = model.atmosphere.heat_capacity

    # Full TOA downward flux; ozone absorption is handled inside the layer loop below.
    D_toa = model.planet.solar_constant * cos_zenith
    D = D_toa

    # DOWNWARD BEAM
    U_reflected = zero(D)
    ozone_absorption = zero(D)

    for k in 1:nlayers
        # 1. cloud reflection?
        if k == cloud_top
            R = cloud_albedo * cloud_cover
            U_reflected = D * R
            D *= (1 - R)
        end

        # 2. ozone absorption in stratosphere layers above σ₀, distribution scaled by layer thickness
        O₃ = O₃_absorption * radiation.ozone_distribution(σ[k]) * Δσ[k]

        # 3. transmissivity of the layer
        D_out = (D - O₃ * D_toa) * t[ij, k]
        # Update temperature tendency due to absorbed shortwave radiation
        # from flux convergence D = D in from the top, D_out at the bottom of a layer
        dTdt[ij, k] += flux_to_tendency((D - D_out) / cₚ, pₛ, k, model)
        D = D_out
    end

    # Surface stratocumulus reflection
    U_stratocumulus = D * stratocumulus_albedo * stratocumulus_cover
    D_surface = D - U_stratocumulus
    vars.parameterizations.surface_shortwave_down[ij] = D_surface
    vars.parameterizations.ocean.surface_shortwave_down[ij] = D_surface
    vars.parameterizations.land.surface_shortwave_down[ij] = D_surface

    # Surface albedo reflections
    up_ocean = albedo_ocean * D_surface
    up_land = albedo_land * D_surface
    vars.parameterizations.ocean.surface_shortwave_up[ij] = up_ocean
    vars.parameterizations.land.surface_shortwave_up[ij] = up_land

    albedo = (1 - land_fraction) * albedo_ocean + land_fraction * albedo_land
    U_surface_albedo = albedo * D_surface
    vars.parameterizations.surface_shortwave_up[ij] = U_surface_albedo
    vars.parameterizations.albedo[ij] = albedo

    U = U_surface_albedo + U_stratocumulus
    for k in nlayers:-1:1
        U_out = U * t[ij, k]
        dTdt[ij, k] += flux_to_tendency((U - U_out) / cₚ, pₛ, k, model)
        U_out += ifelse(k == cloud_top, U_reflected, zero(U))
        U = U_out
    end

    vars.parameterizations.outgoing_shortwave[ij] = U
    return nothing
end

export TwoBandShortwave

"""
    TwoBandShortwave <: AbstractShortwave

A two-band shortwave radiation scheme with diagnostic clouds and ozone following
Fortran SPEEDY (Molteni, 2003; speedy.f90). A visible band (1 - `near_infrared_fraction`
of the incoming solar radiation) is absorbed by ozone in the stratosphere and by dry air,
aerosols, water vapor and clouds in the troposphere, reflected by clouds at the cloud top,
by stratocumulus clouds at the top of the boundary layer and by the surface albedo.
A near-infrared band is only absorbed by water vapor, bypassing cloud reflection,
and fully absorbed at the surface. Fields are

$(TYPEDFIELDS)
"""
@parameterized @kwdef struct TwoBandShortwave{C, T, O, R} <: AbstractShortwave
    @component clouds::C
    @component transmissivity::T
    @component ozone::O
    @component radiative_transfer::R
end
Adapt.@adapt_structure TwoBandShortwave

function TwoBandShortwave(
        SG::SpectralGrid;
        clouds = DiagnosticClouds(SG),
        transmissivity = TwoBandShortwaveTransmissivity(SG),
        ozone = SeasonalOzone(SG),
        radiative_transfer = TwoBandShortwaveRadiativeTransfer(SG),
    )
    return TwoBandShortwave(clouds, transmissivity, ozone, radiative_transfer)
end

Base.show(io::IO, M::TwoBandShortwave) = show(io, M, values = false)

# shortwave variables and those of the ozone component
variables(radiation::TwoBandShortwave) =
    (invoke(variables, Tuple{AbstractShortwave}, radiation)..., variables(radiation.ozone)...)

function initialize!(radiation::TwoBandShortwave, model::PrimitiveEquation)
    initialize!(radiation.clouds, model)
    initialize!(radiation.transmissivity, model)
    initialize!(radiation.ozone, model)
    initialize!(radiation.radiative_transfer, model)
    return nothing
end

# global (non-column) part: update the ozone distribution
function parameterization!(vars::Variables, radiation::TwoBandShortwave, model::PrimitiveEquation)
    ozone_absorption!(vars, radiation.ozone, model)
    return nothing
end

"""$(TYPEDSIGNATURES)
Calculate shortwave radiation using the two-band scheme with diagnostic clouds and ozone."""
@propagate_inbounds function parameterization!(
        ij,
        vars,
        radiation::TwoBandShortwave,
        model,
    )
    clouds = clouds!(ij, vars, radiation.clouds, model)
    t = transmissivity!(ij, vars, clouds, radiation.transmissivity, model)
    shortwave_radiative_transfer!(ij, vars, t, clouds, radiation.ozone, radiation.radiative_transfer, model)
    return nothing
end

export TwoBandShortwaveRadiativeTransfer
"""
    TwoBandShortwaveRadiativeTransfer <: AbstractShortwaveRadiativeTransfer

$(TYPEDFIELDS)"""
@parameterized @kwdef struct TwoBandShortwaveRadiativeTransfer{NF} <: AbstractShortwaveRadiativeTransfer
    "[OPTION] Fraction of incoming solar radiation in the near-infrared band (SPEEDY fband2) [1]"
    @param near_infrared_fraction::NF = 0.05 (bounds = 0 .. 1,)
end
Adapt.@adapt_structure TwoBandShortwaveRadiativeTransfer
TwoBandShortwaveRadiativeTransfer(SG::SpectralGrid; kwargs...) = TwoBandShortwaveRadiativeTransfer{SG.NF}(; kwargs...)
initialize!(::TwoBandShortwaveRadiativeTransfer, ::PrimitiveEquation) = nothing

"""$(TYPEDSIGNATURES)
Two-band shortwave radiative transfer. Visible band: ozone absorption, reflection by clouds
at the top of the cloud-top layer and by stratocumulus at the top of the boundary layer,
absorption with transmissivity `t.visible`, surface albedo reflection and absorption of the
upward beam. Near-infrared band: absorption with transmissivity `t.near_infrared` only,
no reflection by clouds or the surface."""
@propagate_inbounds function shortwave_radiative_transfer!(
        ij,
        vars,
        t,          # NamedTuple from transmissivity!
        clouds,     # NamedTuple from clouds!
        ozone,
        radiation::TwoBandShortwaveRadiativeTransfer,
        model,
    )
    (; cloud_cover, cloud_top, cloud_base, stratocumulus_cover, cloud_albedo, stratocumulus_albedo) = clouds
    (; visible, near_infrared, zenith_factor) = t

    dTdt = get_tendency_step(vars.tendencies.grid.temperature, model.time_stepping, radiation)
    pₛ = vars.parameterizations.surface_pressure[ij]
    nlayers = size(dTdt, 2)

    cos_zenith = vars.parameterizations.cos_zenith[ij]
    albedo_ocean = vars.parameterizations.ocean.albedo[ij]
    albedo_land = vars.parameterizations.land.albedo[ij]
    land_fraction = model.land_sea_mask.land_fraction[ij]
    cₚ = model.atmosphere.heat_capacity

    # top of atmosphere downward flux split into visible and near-infrared band
    D_toa = model.planet.solar_constant * cos_zenith
    D_nir = radiation.near_infrared_fraction * D_toa
    D = D_toa - D_nir

    # DOWNWARD BEAM
    U_cloud = zero(D)               # reflected by clouds at top of cloud-top layer
    U_stratocumulus = zero(D)       # reflected by stratocumulus at top of the boundary layer
    boundary_layer_top = cloud_base + 1     # first layer in the boundary layer

    for k in 1:nlayers
        # 1. reflection (visible only) at the top of layer k
        if k == cloud_top
            U_cloud = D * cloud_albedo * cloud_cover
            D -= U_cloud
        end
        if k == boundary_layer_top
            U_stratocumulus = D * stratocumulus_albedo * stratocumulus_cover
            D -= U_stratocumulus
        end

        # 2. ozone absorption (visible), fraction of TOA flux, corrected for slant path
        O₃ = min(D, ozone_absorption(ij, k, vars, ozone, model) * D_toa * zenith_factor)

        # 3. absorption in both bands
        D_out = (D - O₃) * visible[ij, k]
        D_nir_out = D_nir * near_infrared[ij, k]
        absorbed = (D - D_out) + (D_nir - D_nir_out)
        dTdt[ij, k] += flux_to_tendency(absorbed / cₚ, pₛ, k, model)
        D, D_nir = D_out, D_nir_out
    end

    # SURFACE, both bands reach the surface but only visible is reflected
    D_surface = D + D_nir
    vars.parameterizations.surface_shortwave_down[ij] = D_surface
    vars.parameterizations.ocean.surface_shortwave_down[ij] = D_surface
    vars.parameterizations.land.surface_shortwave_down[ij] = D_surface

    vars.parameterizations.ocean.surface_shortwave_up[ij] = albedo_ocean * D
    vars.parameterizations.land.surface_shortwave_up[ij] = albedo_land * D

    albedo = (1 - land_fraction) * albedo_ocean + land_fraction * albedo_land
    U = albedo * D
    vars.parameterizations.surface_shortwave_up[ij] = U
    vars.parameterizations.albedo[ij] = albedo

    # UPWARD BEAM (visible only), add reflected fluxes above the layer where they were reflected
    for k in nlayers:-1:1
        U_out = U * visible[ij, k]
        dTdt[ij, k] += flux_to_tendency((U - U_out) / cₚ, pₛ, k, model)
        U_out += ifelse(k == cloud_top, U_cloud, zero(U))
        U_out += ifelse(k == boundary_layer_top, U_stratocumulus, zero(U))
        U = U_out
    end

    vars.parameterizations.outgoing_shortwave[ij] = U
    return nothing
end
