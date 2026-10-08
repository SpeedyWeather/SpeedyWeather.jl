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

# variables of the shortwave scheme and those of its clouds component
# variables of the shortwave scheme and those of its clouds and radiative transfer components
variables(radiation::OneBandShortwave) = (
    invoke(variables, Tuple{AbstractShortwave}, radiation)...,
    variables(radiation.clouds)...,
    variables(radiation.radiative_transfer)...,
)

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

export CloudyShortwaveRadiativeTransfer

"""
    CloudyShortwaveRadiativeTransfer <: AbstractShortwaveRadiativeTransfer

One-band shortwave radiative transfer with clouds in every layer from the cloud state of a
prognostic cloud scheme, e.g. [`PrognosticCloudCondensation`](@ref), solved with the adding method
(multiple reflections between layers and the surface). In each layer the cloudy part, cloud
fraction `C`, has the reflectance `R_c` and transmittance `T_c` of a homogeneous, absorbing,
diffusely illuminated two-stream layer (quadrature coefficients) with the in-cloud optical depth

    τ = 3/2 (LWP/(ρ_w r_l) + IWP/(ρ_i r_i)) / C

and the layer's reflectance and transmittance are `R = C R_c` and `T = t ((1 - C) + C T_c)`
with the clear-sky transmissivity `t` (gas absorption). Clouds of different layers are
independent (random overlap). Ozone absorbs a fraction of the incoming beam above `σ = 0.2`
as in [`OneBandShortwaveRadiativeTransfer`](@ref). Use with clouds that do not reflect at a cloud
top themselves and a transmissivity without cloud absorption, see [`OneBandCloudyShortwave`](@ref).
Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct CloudyShortwaveRadiativeTransfer{NF, F} <: AbstractShortwaveRadiativeTransfer
    "[OPTION] Asymmetry factor of the scattering by cloud droplets and ice [1]"
    @param asymmetry_factor::NF = 0.85 (bounds = 0 .. 1,)

    "[OPTION] Single-scattering albedo of cloud droplets and ice, broadband [1]"
    @param single_scattering_albedo::NF = 0.999 (bounds = 0 .. 1,)

    "[OPTION] Density of ice [kg/m³] for the optical depth of cloud ice"
    @param ice_density::NF = 917 (bounds = Positive,)

    "[OPTION] Total ozone absorption as fraction of incoming solar radiation (1)"
    @param ozone_absorption::NF = 0.01 (bounds = 0 .. 1,)

    "[OPTION] Ozone distribution above σ₀, has to be explicitly normalized to ∫dσ = 1 (1)"
    ozone_distribution::F
end

Adapt.@adapt_structure CloudyShortwaveRadiativeTransfer

function CloudyShortwaveRadiativeTransfer(
        SG::SpectralGrid;
        ozone_distribution = (σ) -> 50 * max(0, 1 // 5 - σ),     # as OneBandShortwaveRadiativeTransfer
        kwargs...
    )
    return CloudyShortwaveRadiativeTransfer{SG.NF, typeof(ozone_distribution)}(; ozone_distribution, kwargs...)
end

initialize!(::CloudyShortwaveRadiativeTransfer, ::PrimitiveEquation) = nothing

# the shortwave diagnostics and the cloud state it reads
variables(radiative_transfer::CloudyShortwaveRadiativeTransfer) = (
    invoke(variables, Tuple{AbstractShortwave}, radiative_transfer)...,
    cloud_state_variables()...,
)

export OneBandCloudyShortwave

"""$(TYPEDSIGNATURES)
`OneBandShortwave` with clouds in every layer from a prognostic cloud scheme: [`PrognosticClouds`](@ref)
(column cloud cover and cloud top for output only), the background transmissivity without its
cloud absorption, and [`CloudyShortwaveRadiativeTransfer`](@ref)."""
function OneBandCloudyShortwave(
        SG::SpectralGrid;
        clouds = PrognosticClouds(SG),
        transmissivity = BackgroundShortwaveTransmissivity(SG; absorptivity_cloud_base = 0, absorptivity_cloud_limit = 0),
        radiative_transfer = CloudyShortwaveRadiativeTransfer(SG),
    )
    return OneBandShortwave(clouds, transmissivity, radiative_transfer)
end

"""$(TYPEDSIGNATURES)
Reflectance and transmittance of a homogeneous layer of optical depth `τ`, single-scattering albedo
`ω` and asymmetry factor `g` for diffuse illumination, from the two-stream equations with
quadrature coefficients `γ₁ = √3/2 (2 - ω(1 + g))`, `γ₂ = √3/2 ω(1 - g)` and `k = √(γ₁² - γ₂²)`

    R = γ₂ (1 - e^{-2kτ}) / (k (1 + e^{-2kτ}) + γ₁ (1 - e^{-2kτ}))
    T = 2k e^{-kτ} / (k (1 + e^{-2kτ}) + γ₁ (1 - e^{-2kτ}))

written with `expm1` so that the non-absorbing limit `R = γ₁τ/(1 + γ₁τ)`, `T = 1 - R` is
reached without cancellation."""
@inline function two_stream_diffuse_layer(τ, ω, g)
    NF = typeof(τ)
    γ₁ = sqrt(NF(3)) / 2 * (2 - ω * (1 + g))
    γ₂ = sqrt(NF(3)) / 2 * ω * (1 - g)
    k = max(sqrt(max(γ₁^2 - γ₂^2, zero(NF))), NF(1.0e-6))    # floor for ω = 1
    x = -expm1(-2k * τ)                                     # 1 - e^{-2kτ}
    denominator = k * (2 - x) + γ₁ * x
    R = γ₂ * x / denominator
    T = 2k * exp(-k * τ) / denominator
    return R, T
end

"""$(TYPEDSIGNATURES)
One-band shortwave radiative transfer of column `ij` with clouds in every layer, adding method,
see [`CloudyShortwaveRadiativeTransfer`](@ref). `t` is the clear-sky transmissivity (scratch
array, overwritten), `vars.scratch.grid.b` is used as work array."""
@propagate_inbounds function shortwave_radiative_transfer!(
        ij,
        vars,
        t,          # transmissivity array, overwritten
        clouds,     # NamedTuple from clouds!, not used for the transfer
        radiation::CloudyShortwaveRadiativeTransfer,
        model,
    )
    (; cloud_fraction, cloud_liquid_water, cloud_ice_water) = vars.parameterizations
    (; cloud_liquid_effective_radius, cloud_ice_effective_radius) = vars.parameterizations
    (; asymmetry_factor, single_scattering_albedo, ice_density, ozone_absorption) = radiation
    albedo_stack = vars.scratch.grid.b          # reflectance of layers, then albedo of the stack below

    dTdt = get_tendency_step(vars.tendencies.grid.temperature, model.time_stepping, radiation)
    NF = eltype(dTdt)
    pₛ = vars.parameterizations.surface_pressure[ij]
    nlayers = size(dTdt, 2)
    σ = model.geometry.σ_levels_full
    Δσ = model.geometry.σ_levels_thick
    coord = model.geometry.vertical_coordinates
    g = model.planet.gravity
    cₚ = model.atmosphere.heat_capacity
    ρ_water = model.atmosphere.water_density
    r_min = NF(1.0e-6)                          # [m], effective radius floor, also for absent cloud state

    cos_zenith = vars.parameterizations.cos_zenith[ij]
    albedo_ocean = vars.parameterizations.ocean.albedo[ij]
    albedo_land = vars.parameterizations.land.albedo[ij]
    land_fraction = model.land_sea_mask.land_fraction[ij]
    albedo = (1 - land_fraction) * albedo_ocean + land_fraction * albedo_land

    # 1. OZONE absorbs from the incoming beam in the stratosphere, LAYER OPTICS of the cloudy layers
    D_toa = model.planet.solar_constant * cos_zenith
    D = D_toa
    for k in 1:nlayers
        O₃ = ozone_absorption * radiation.ozone_distribution(σ[k]) * Δσ[k] * D_toa
        dTdt[ij, k] += flux_to_tendency(O₃ / cₚ, pₛ, k, model)
        D -= O₃

        C = cloud_fraction[ij, k]
        Δp_g = pressure_thickness(k, pₛ, coord) / g     # layer mass [kg/m²]
        r_liquid = max(cloud_liquid_effective_radius[ij, k], r_min)
        r_ice = max(cloud_ice_effective_radius[ij, k], r_min)
        τ_grid = 3 * Δp_g / 2 * (cloud_liquid_water[ij, k] / (ρ_water * r_liquid) + cloud_ice_water[ij, k] / (ice_density * r_ice))
        τ = τ_grid / max(C, eps(NF))                    # in-cloud optical depth
        R_cloud, T_cloud = two_stream_diffuse_layer(τ, single_scattering_albedo, asymmetry_factor)
        albedo_stack[ij, k] = C * R_cloud               # layer reflectance
        t[ij, k] *= (1 - C) + C * T_cloud               # layer transmittance
    end

    # 2. ADDING from the surface up: albedo of the stack below each layer top, and the transmission
    # of a downward flux through layer k including the multiple reflections with the stack below
    A = albedo
    for k in nlayers:-1:1
        R = albedo_stack[ij, k]
        T = t[ij, k]
        multiple_reflections = inv(1 - R * A)
        A = R + T^2 * A * multiple_reflections
        albedo_stack[ij, k] = A
        t[ij, k] = T * multiple_reflections
    end

    # 3. FLUXES from the top down, absorption in each layer from the net flux convergence
    D_top = D                                   # downward flux into the top layer after ozone
    for k in 1:nlayers
        U = albedo_stack[ij, k] * D             # upward at the top of layer k
        D_below = D * t[ij, k]
        A_below = ifelse(k < nlayers, albedo_stack[ij, min(k + 1, nlayers)], albedo)
        U_below = A_below * D_below
        absorbed = (D - U) - (D_below - U_below)
        dTdt[ij, k] += flux_to_tendency(absorbed / cₚ, pₛ, k, model)
        D = D_below
    end

    # surface fluxes, D is now the downward flux at the surface
    vars.parameterizations.surface_shortwave_down[ij] = D
    vars.parameterizations.ocean.surface_shortwave_down[ij] = D
    vars.parameterizations.land.surface_shortwave_down[ij] = D
    vars.parameterizations.ocean.surface_shortwave_up[ij] = albedo_ocean * D
    vars.parameterizations.land.surface_shortwave_up[ij] = albedo_land * D
    vars.parameterizations.surface_shortwave_up[ij] = albedo * D
    vars.parameterizations.albedo[ij] = albedo
    vars.parameterizations.outgoing_shortwave[ij] = albedo_stack[ij, 1] * D_top
    return nothing
end
