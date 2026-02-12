abstract type AbstractRadiation <: AbstractParameterization end
abstract type AbstractShortwave <: AbstractRadiation end
abstract type AbstractShortwaveRadiativeTransfer <: AbstractShortwave end

export TransparentShortwave
struct TransparentShortwave <: AbstractShortwave end
Adapt.@adapt_structure TransparentShortwave
TransparentShortwave(SG::SpectralGrid) = TransparentShortwave()

function variables(::AbstractShortwave)
    return (
        DiagnosticVariable(name = :surface_shortwave_down, dims = Grid2D(), desc = "Surface shortwave radiation down", units = "W/m^2"),
        DiagnosticVariable(name = :surface_shortwave_down, dims = Grid2D(), desc = "Surface shortwave radiation down over ocean", units = "W/m^2", namespace = :ocean),
        DiagnosticVariable(name = :surface_shortwave_down, dims = Grid2D(), desc = "Surface shortwave radiation down over land", units = "W/m^2", namespace = :land),
        DiagnosticVariable(name = :surface_shortwave_up, dims = Grid2D(), desc = "Surface shortwave radiation up", units = "W/m^2"),
        DiagnosticVariable(name = :surface_shortwave_up, dims = Grid2D(), desc = "Surface shortwave radiation up over ocean", units = "W/m^2", namespace = :ocean),
        DiagnosticVariable(name = :surface_shortwave_up, dims = Grid2D(), desc = "Surface shortwave radiation up over land", units = "W/m^2", namespace = :land),
        DiagnosticVariable(name = :outgoing_shortwave, dims = Grid2D(), desc = "TOA Shortwave radiation up", units = "W/m^2"),
        DiagnosticVariable(name = :cos_zenith, dims = Grid2D(), desc = "Cos zenith angle", units = "1"),
        DiagnosticVariable(name = :albedo, dims = Grid2D(), desc = "Albedo", units = "1"),
        DiagnosticVariable(name = :albedo, dims = Grid2D(), desc = "Albedo over ocean", units = "1", namespace = :ocean),
        DiagnosticVariable(name = :albedo, dims = Grid2D(), desc = "Albedo over land", units = "1", namespace = :land),
    )
end

initialize!(::TransparentShortwave, ::PrimitiveEquation) = nothing

# function barrier
@propagate_inbounds function parameterization!(ij, diagn, progn, ::TransparentShortwave, model)

    planet = model.planet
    land_sea_mask = model.land_sea_mask.mask

    (; surface_shortwave_down, surface_shortwave_up) = diagn.physics
    ssrd_ocean = diagn.physics.ocean.surface_shortwave_down
    ssrd_land = diagn.physics.land.surface_shortwave_down
    ssru_ocean = diagn.physics.ocean.surface_shortwave_up
    ssru_land = diagn.physics.land.surface_shortwave_up

    cos_zenith = diagn.physics.cos_zenith[ij]
    land_fraction = land_sea_mask[ij]
    albedo_ocean = diagn.physics.ocean.albedo[ij]
    albedo_land = diagn.physics.land.albedo[ij]

    S₀ = planet.solar_constant
    D = S₀ * cos_zenith             # top of atmosphere downward radiation
    surface_shortwave_down[ij] = D  # transparent atmosphere so same at surface (before albedo)
    ssrd_ocean[ij] = D
    ssrd_land[ij] = D

    # shortwave up is after albedo reflection, separated by ocean/land
    ssru_ocean[ij] = albedo_ocean * D
    ssru_land[ij] = albedo_land * D

    # land-sea mask-weighted
    albedo = (1 - land_fraction) * albedo_ocean + land_fraction * albedo_land
    surface_shortwave_up[ij] = albedo * D
    diagn.physics.albedo[ij] = albedo   # store weighted albedo

    # transparent also for reflected shortwave radiation travelling up
    diagn.physics.outgoing_shortwave[ij] = surface_shortwave_up[ij]
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
        transmissivity = TransparentShortwaveTransmissivity(SG),
        radiative_transfer = OneBandShortwaveRadiativeTransfer(SG),
    )
    return OneBandShortwave(clouds, transmissivity, radiative_transfer)
end

function Base.show(io::IO, M::OneBandShortwave)
    println(io, "OneBandShortwave <: AbstractShortwave")
    properties = propertynames(M)
    n = length(properties)
    for (i, key) in enumerate(properties)
        val = getfield(M, key)
        s = i == n ? "└" : "├"  # choose ending └ for last property
        p = i == n ? print : println
        p(io, "$s $key: $(typeof(val))")
    end
    return
end

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
        diagn,
        progn,
        radiation::OneBandShortwave,
        model,
    )
    clouds = clouds!(ij, diagn, progn, radiation.clouds, model)
    t = transmissivity!(ij, diagn, progn, clouds, radiation.transmissivity, model)
    shortwave_radiative_transfer!(ij, diagn, t, clouds, radiation.radiative_transfer, model)
    return nothing
end

export OneBandShortwaveRadiativeTransfer
"""
    OneBandShortwaveRadiativeTransfer <: AbstractShortwaveRadiativeTransfer

$(TYPEDFIELDS)."""
@parameterized @kwdef struct OneBandShortwaveRadiativeTransfer{NF} <: AbstractShortwaveRadiativeTransfer
    "[OPTION] Ozone absorption in upper stratosphere (layer 1) in fraction of incoming solar radiation (1)"
    @param ozone_absorption_upper::NF = 0.03 (bounds = 0..1,)

    "[OPTION] Ozone absorption in lower stratosphere (layer 2) in fraction of incoming solar radiation (1)"
    @param ozone_absorption_lower::NF = 0.024 (bounds = 0..1,)
end
Adapt.@adapt_structure OneBandShortwaveRadiativeTransfer

# generator function
OneBandShortwaveRadiativeTransfer(SG::SpectralGrid; kwargs...) = OneBandShortwaveRadiativeTransfer{SG.NF}(; kwargs...)
initialize!(::OneBandShortwaveRadiativeTransfer, ::PrimitiveEquation) = nothing

"""$(TYPEDSIGNATURES)
One-band shortwave radiative transfer with cloud reflection and ozone absorption."""
@propagate_inbounds function shortwave_radiative_transfer!(
        ij,
        diagn,
        t,          # Transmissivity array
        clouds,     # NamedTuple from clouds!
        radiation::OneBandShortwaveRadiativeTransfer,
        model,
    )
    (; ozone_absorption_upper, ozone_absorption_lower) = radiation
    (; cloud_cover, cloud_top, stratocumulus_cover, cloud_albedo, stratocumulus_albedo) = clouds

    dTdt = diagn.tendencies.temp_tend_grid
    pₛ = diagn.grid.pres_grid_prev[ij]
    nlayers = size(dTdt, 2)

    cos_zenith = diagn.physics.cos_zenith[ij]
    albedo_ocean = diagn.physics.ocean.albedo[ij]
    albedo_land = diagn.physics.land.albedo[ij]
    land_fraction = model.land_sea_mask.mask[ij]
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

        # 2. ozone absorption in stratosphere layers
        # TODO make this vertical coordinate and not layer index dependent
        if k == 1
            # ozone absorption is a fraction of incoming solar radiation at TOA
            # not dependent on latitude or year as done in Fortran SPEEDY
            ozone_absorption = ozone_absorption_upper * D_toa
        elseif k == 2
            # similar here but lower stratosphere absorption is typically weaker
            ozone_absorption = ozone_absorption_lower * D_toa
        else
            # no ozone absorption in troposphere layers
            ozone_absorption = zero(D)
        end

        # 3. transmissivity of the layer
        D_out = (D - ozone_absorption) * t[ij, k]
        # Update temperature tendency due to absorbed shortwave radiation
        # from flux convergence D = D in from the top, D_out at the bottom of a layer
        dTdt[ij, k] += flux_to_tendency((D - D_out) / cₚ, pₛ, k, model)
        D = D_out
    end

    # Surface stratocumulus reflection
    U_stratocumulus = D * stratocumulus_albedo * stratocumulus_cover
    D_surface = D - U_stratocumulus
    diagn.physics.surface_shortwave_down[ij] = D_surface
    diagn.physics.ocean.surface_shortwave_down[ij] = D_surface
    diagn.physics.land.surface_shortwave_down[ij] = D_surface

    # Surface albedo reflections
    up_ocean = albedo_ocean * D_surface
    up_land = albedo_land * D_surface
    diagn.physics.ocean.surface_shortwave_up[ij] = up_ocean
    diagn.physics.land.surface_shortwave_up[ij] = up_land

    albedo = (1 - land_fraction) * albedo_ocean + land_fraction * albedo_land
    U_surface_albedo = albedo * D_surface
    diagn.physics.surface_shortwave_up[ij] = U_surface_albedo
    diagn.physics.albedo[ij] = albedo

    U = U_surface_albedo + U_stratocumulus
    for k in nlayers:-1:1
        U_out = U * t[ij, k]
        dTdt[ij, k] += flux_to_tendency((U - U_out) / cₚ, pₛ, k, model)
        U_out += k == cloud_top ? U_reflected : zero(U)
        U = U_out
    end

    diagn.physics.outgoing_shortwave[ij] = U
    return nothing
end
