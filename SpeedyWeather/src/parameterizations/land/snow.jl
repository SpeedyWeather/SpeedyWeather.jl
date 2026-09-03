abstract type AbstractSnow <: AbstractLandComponent end

export SnowModel    # maybe change for a more concise name later

"""
$(TYPEDSIGNATURES)
Single-column snow bucket model in equivalent liquid water depth. Snow accumulates
from the diagnosed precipitation, melts once the top soil layer exceeds
`melting_threshold`, and is capped at `snow_depth_cap` to limit infinite snow/ice accumulation
over perennial ice caps and glaciers.

Optionally (`sea_ice_snow = true`) a second snow bucket accumulates over sea ice, where it
raises the surface albedo and insulates the surface heat and humidity fluxes. That bucket melts
with the surface air temperature above `melting_threshold` and is removed entirely once the sea
ice concentration drops below `sea_ice_threshold`, i.e. the snow disappears with the ice.
$(TYPEDFIELDS)"""
@parameterized @kwdef struct SnowModel{NF} <: AbstractSnow
    "[OPTION] Temperature threshold for snow melting [K]"
    @param melting_threshold::NF = 275 (bounds = Positive,)

    "[OPTION] Permanent snow/ice depth cap in equivalent liquid water depth [m]"
    snow_depth_cap::NF = 10

    "[OPTION] Accumulate snow on sea ice, affecting albedo and surface fluxes"
    sea_ice_snow::Bool = false

    "[OPTION] Melt rate of snow on sea ice per degree above the melting threshold [m/s/K]"
    @param sea_ice_melt_factor::NF = 1.0e-7 (bounds = Nonnegative,)

    "[OPTION] Sea ice concentration below which snow on sea ice is removed entirely [1]"
    @param sea_ice_threshold::NF = 0.01 (bounds = 0 .. 1,)
end

Adapt.@adapt_structure SnowModel

# generator function
SnowModel(SG::SpectralGrid, geometry::LandGeometryOrNothing = nothing; kwargs...) = SnowModel{SG.NF}(; kwargs...)

function variables(snow::SnowModel)
    vars = (
        PrognosticVariable(:snow_depth, Grid2D(), namespace = :land, units = "m", desc = "Snow depth in equivalent liquid water height"),
        PrognosticVariable(:soil_temperature, LandXYZ(), namespace = :land, units = "K", desc = "Soil temperature"),
        ParameterizationVariable(:snow_melt_rate, Grid2D(), namespace = :land, units = "kg/m²/s", desc = "Snow melt rate"),
    )

    # only allocate the snow on sea ice bucket if that effect is enabled
    snow.sea_ice_snow || return vars
    return (
        vars...,
        PrognosticVariable(:snow_depth, Grid2D(), namespace = :ocean, units = "m", desc = "Snow depth on sea ice in equivalent liquid water height"),
    )
end

# initialize component
initialize!(snow::SnowModel, model::PrimitiveEquation) = nothing

# set initial conditions for snow depth in initial conditions
function initialize!(vars::Variables, snow::SnowModel, model::PrimitiveEquation)
    set!(vars.prognostic.land, model.geometry, snow_depth = 0)
    snow.sea_ice_snow && haskey(vars.prognostic.ocean, :snow_depth) &&
        fill!(vars.prognostic.ocean.snow_depth, 0)
    return nothing
end

function timestep!(
        vars::Variables,
        snow::SnowModel,
        model::PrimitiveEquation,
    )

    (; Δt) = model.time_stepping                            # time step [s]
    (; snow_depth) = vars.prognostic.land                   # in equivalent liquid water height [m]
    (; soil_temperature) = vars.prognostic.land
    (; land_fraction) = model.land_sea_mask

    # Some thermodynamics needed by snow
    ρ_water = model.atmosphere.water_density                # water density [kg/m³]
    Lᵢ = model.atmosphere.latent_heat_fusion                # latent heat of fusion
    cₛ = model.land.thermodynamics.heat_capacity_dry_soil
    z₁ = model.land.geometry.layer_thickness[1]
    (; melting_threshold, snow_depth_cap) = snow

    # reset in any case
    vars.scratch.grid.a_2D .= 0

    # from precipitation schemes [m/s]
    snow_fall_rate = haskey(vars.parameterizations, :snow_rate) ?
        vars.parameterizations.snow_rate :
        vars.scratch.grid.a_2D

    snow_melt_rate = vars.parameterizations.land.snow_melt_rate     # for soil moisture model

    params = (; melting_threshold, cₛ, z₁, Δt, ρ_water, Lᵢ, snow_depth_cap)

    launch!(
        architecture(snow_depth), LinearWorkOrder, size(snow_depth), land_snow_kernel!,
        snow_depth, soil_temperature, snow_melt_rate, snow_fall_rate, land_fraction,
        params,
    )

    # snow on sea ice, only if enabled and there is sea ice to accumulate it on
    snow.sea_ice_snow && sea_ice_snow_timestep!(vars, snow, snow_fall_rate, model)
    return nothing
end

"""$(TYPEDSIGNATURES)
Time step the snow bucket on sea ice, a no-op unless both the snow on sea ice variable and a
sea ice concentration are available."""
function sea_ice_snow_timestep!(
        vars::Variables,
        snow::SnowModel,
        snow_fall_rate,
        model::PrimitiveEquation,
    )
    haskey(vars.prognostic.ocean, :snow_depth) || return nothing
    haskey(vars.prognostic.ocean, :sea_ice_concentration) || return nothing

    (; Δt) = model.time_stepping                            # time step [s]
    (; land_fraction) = model.land_sea_mask
    snow_depth = vars.prognostic.ocean.snow_depth           # per sea ice area [m]

    # sea ice concentration at the current prognostic step, hence lagged by one sea ice update
    ℵ = get_prognostic_step(vars.prognostic.ocean.sea_ice_concentration, model.time_stepping, model.sea_ice)

    # surface air temperature drives the melt, no surface energy budget over sea ice available
    (; surface_air_temperature) = vars.parameterizations
    (; melting_threshold, snow_depth_cap, sea_ice_melt_factor, sea_ice_threshold) = snow

    params = (; melting_threshold, snow_depth_cap, sea_ice_melt_factor, sea_ice_threshold, Δt)

    launch!(
        architecture(snow_depth), LinearWorkOrder, size(snow_depth), sea_ice_snow_kernel!,
        snow_depth, ℵ, surface_air_temperature, snow_fall_rate, land_fraction,
        params,
    )
    return nothing
end

@kernel inbounds = true function sea_ice_snow_kernel!(
        snow_depth, ℵ, surface_air_temperature, snow_fall_rate, land_fraction,
        params,
    )
    ij = @index(Global, Linear)             # every grid point ij

    if land_fraction[ij] < 1                # at least partially ocean

        (; melting_threshold, snow_depth_cap, sea_ice_melt_factor, sea_ice_threshold, Δt) = params

        # snow depth is a depth per sea ice area, so accumulation is not weighted by concentration
        # Term 1: snow fall rate from precipitation schemes [m/s]
        snow_fall = snow_fall_rate[ij]

        # Term 2: melt rate above the melting threshold [m/s], a degree-day-style melt as there is
        # no prognostic sea ice surface temperature to close a surface energy budget with
        T = surface_air_temperature[ij]
        δT_melt = isfinite(T) ? max(T - melting_threshold, 0) : zero(T)
        melt_rate = sea_ice_melt_factor * δT_melt

        # don't melt more snow than there is, so the bucket cannot go negative
        melt_rate = min(snow_depth[ij] / Δt, melt_rate)

        # Euler forward time step, floored at 0 depth and capped as over land
        snow_depth_forward = min(max(snow_depth[ij] + Δt * (snow_fall - melt_rate), 0), snow_depth_cap)

        # snow disappears with the sea ice it sits on. Mass is not conserved here, the melt water
        # is dropped rather than added to an ocean freshwater budget, which the model doesn't track
        snow_depth[ij] = ifelse(ℵ[ij] < sea_ice_threshold, zero(snow_depth_forward), snow_depth_forward)
    end
end

@kernel inbounds = true function land_snow_kernel!(
        snow_depth, soil_temperature, snow_melt_rate, snow_fall_rate, land_fraction,
        params,
    )
    ij = @index(Global, Linear)             # every grid point ij

    if land_fraction[ij] > 0               # at least partially land

        (; melting_threshold, cₛ, z₁, Δt, ρ_water, Lᵢ, snow_depth_cap) = params

        # check for melting of snow if temperature above melting threshold
        # check for NaNs here to prevent land temperatures read from NetCDF data
        # to cause an immediate blow up in case the land-sea mask doesn't align
        δT_melt = isfinite(soil_temperature[ij, 1]) ?
            max(soil_temperature[ij, 1] - melting_threshold, 0) : zero(soil_temperature[ij, 1])

        # energy available from soil warming above melting threshold [J/m²/s]
        # heat capacity per volume, so not *density needed
        E_avail = cₛ * δT_melt * z₁ / Δt  # [J/(m³ K)] * [K] * [m] / [s] = [J/m²/s]

        # Term 1: snow fall rate from precipitation schemes [m/s]
        snow_fall_rate_max = snow_fall_rate[ij]

        # Term 2: max melt rate allowed by available energy [m/s]
        melt_rate_max = E_avail / (ρ_water * Lᵢ)

        # Adding the terms change snow depth by falling snow minus melting and runoff [m/s]
        # maximum amount of snow change
        dsnow_max = snow_fall_rate_max - melt_rate_max

        # don't melt or runoff more than there is snow + snow that is falling down this time step
        # limited amount of snow change by how much there is
        dsnow = -min(snow_depth[ij] / Δt, -dsnow_max)

        # snow that we tried to melt/runoff that isn't available though
        dsnow_excess = dsnow_max - dsnow

        # store to pass to soil moisture [kg/m²/s], combined runoff with melt rate
        # limited to what's available to melt/runoff
        snow_melt_rate[ij] = (melt_rate_max + dsnow_excess) * ρ_water

        # Euler forward time step but cap at 0 depth to not melt more snow than available
        snow_depth_forward = max(snow_depth[ij] + Δt * dsnow, 0)

        # Conservation of mass is violated here by removing excess snow above the depth cap,
        # but we do it at the end here to avoid that excess adding to the melt rate
        snow_depth[ij] = min(snow_depth_forward, snow_depth_cap)
    end
end
