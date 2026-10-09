abstract type AbstractSnow <: AbstractLandComponent end
abstract type AbstractDynamicSnow <: AbstractSnow end
abstract type AbstractPrescribedSnow <: AbstractSnow end

export SnowModel    # maybe change for a more concise name later

"""
$(TYPEDSIGNATURES)
Single-column snow bucket model in equivalent liquid water depth. Snow accumulates
from the diagnosed precipitation, melts once the top soil layer exceeds
`melting_threshold`, and is capped at `snow_depth_cap` to limit infinite snow/ice accumulation
over perennial ice caps and glaciers. The snow depth is time stepped by the land time stepper
(`time_stepper(model.time_stepping, :land)`), the cap (and non-negativity) is applied in `filter!`.
$(TYPEDFIELDS)"""
@parameterized @kwdef struct SnowModel{NF} <: AbstractDynamicSnow
    "[OPTION] Temperature threshold for snow melting [K]"
    @param melting_threshold::NF = 275 (bounds = Positive,)

    "[OPTION] Permanent snow/ice depth cap in equivalent liquid water depth [m]"
    snow_depth_cap::NF = 10
end

Adapt.@adapt_structure SnowModel

# generator function
SnowModel(SG::SpectralGrid, geometry::LandGeometryOrNothing = nothing; kwargs...) = SnowModel{SG.NF}(; kwargs...)

# snow reads but doesn't define soil_temperature, whose dimensions are defined by the land temperature model
function variables(::SnowModel, model::AbstractModel)
    nsteps = get_nsteps(model.time_stepping, :land)
    pg = nsteps.prognostic
    tg = nsteps.tendency
    return (
        PrognosticVariable(:snow_depth, LandXYT(pg), namespace = :land, units = "m", desc = "Snow depth in equivalent liquid water height"),
        TendencyVariable(:snow_depth, LandXYT(tg), namespace = :land, units = "m/s", desc = "Tendency of snow depth"),
        ParameterizationVariable(:snow_melt_rate, Grid2D(), namespace = :land, units = "kg/m²/s", desc = "Snow melt rate"),
    )
end

# initialize component
initialize!(snow::SnowModel, model::PrimitiveEquation) = nothing

# set initial conditions for snow depth in initial conditions
initialize!(vars::Variables, snow::SnowModel, model::PrimitiveEquation) =
    set!(vars.prognostic.land, model.geometry, snow_depth = 0)

function timestep!(
        vars::Variables,
        snow::SnowModel,
        model::PrimitiveEquation,
    )

    land_time_stepper = time_stepper(model.time_stepping, :land)
    soil_temperature = get_prognostic_step(vars.prognostic.land.soil_temperature, land_time_stepper, snow)
    snow_depth = get_prognostic_step(vars.prognostic.land.snow_depth, land_time_stepper, snow)  # equivalent liquid water height [m]
    snow_depth_tendency = get_tendency_step(vars.tendencies.land.snow_depth, land_time_stepper, snow)

    Δt = time_step(model.time_stepping, :land, vars.prognostic.clock)    # time step [s] of the land time stepper
    (; land_fraction) = model.land_sea_mask

    # Some thermodynamics needed by snow
    ρ_water = model.atmosphere.water_density                # water density [kg/m³]
    Lᵢ = model.atmosphere.latent_heat_fusion                # latent heat of fusion
    cₛ = model.land.thermodynamics.heat_capacity_dry_soil
    z₁ = model.land.geometry.layer_thickness[1]
    (; melting_threshold) = snow

    # reset in any case
    vars.scratch.grid.a_2D .= 0

    # from precipitation schemes [m/s]
    snow_fall_rate = haskey(vars.parameterizations, :snow_rate) ?
        vars.parameterizations.snow_rate :
        vars.scratch.grid.a_2D

    snow_melt_rate = vars.parameterizations.land.snow_melt_rate     # for soil moisture model

    params = (; melting_threshold, cₛ, z₁, Δt, ρ_water, Lᵢ)

    launch!(
        architecture(snow_depth), LinearWorkOrder, size(snow_depth), land_snow_kernel!,
        snow_depth_tendency, snow_depth, soil_temperature, snow_melt_rate, snow_fall_rate, land_fraction,
        params,
    )
    return nothing
end

@kernel inbounds = true function land_snow_kernel!(
        snow_depth_tendency, snow_depth, soil_temperature, snow_melt_rate, snow_fall_rate, land_fraction,
        params,
    )
    ij = @index(Global, Linear)             # every grid point ij

    if land_fraction[ij] > 0               # at least partially land

        (; melting_threshold, cₛ, z₁, Δt, ρ_water, Lᵢ) = params

        # check for melting of snow if temperature above melting threshold
        δT_melt = max(soil_temperature[ij, 1] - melting_threshold, 0)

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
        # limited amount of snow change by how much there is (Δt is the land time step)
        dsnow = -min(snow_depth[ij] / Δt, -dsnow_max)

        # snow that we tried to melt/runoff that isn't available though
        dsnow_excess = dsnow_max - dsnow

        # store to pass to soil moisture [kg/m²/s], combined runoff with melt rate
        # limited to what's available to melt/runoff
        snow_melt_rate[ij] = (melt_rate_max + dsnow_excess) * ρ_water

        # the land time stepper steps snow depth forward with the uncapped tendency (the cap to the
        # available snow is not a tendency but an adjustment, done in filter!), with Euler forward
        # this is identical to the capped tendency as the melt rate above is capped
        snow_depth_tendency[ij] = dsnow_max
    end
end

"""$(TYPEDSIGNATURES)
Keep snow depth within [0, `snow_depth_cap`] after the time step. Conservation of mass is
violated here by removing excess snow above the depth cap, done after the time step so that
this excess does not add to the melt rate."""
function Base.filter!(vars::Variables, snow::SnowModel, model::PrimitiveEquation)
    snow_depth = vars.prognostic.land.snow_depth
    (; snow_depth_cap) = snow
    snow_depth .= clamp.(snow_depth, 0, snow_depth_cap)     # all steps
    return nothing
end
