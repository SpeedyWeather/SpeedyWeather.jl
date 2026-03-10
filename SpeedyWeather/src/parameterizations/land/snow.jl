abstract type AbstractSnow <: AbstractLandComponent end

export SnowModel    # maybe change for a more concise name later

"""
    SnowModel(; melting_threshold=275, runoff_time_scale=Year(1))

Single-column snow bucket model in equivalent liquid water depth. Snow accumulates
from the diagnosed precipitation, melts once the top soil layer exceeds
`melting_threshold`, and relaxes back to zero on the `runoff_time_scale` to prevent
unrealistic snow buildup as in reality the snow would turn into ice and a glacier / 
ice sheet would form instead.
$(TYPEDFIELDS)"""
@parameterized @kwdef struct SnowModel{NF} <: AbstractSnow
    "[OPTION] Temperature threshold for snow melting [K]"
    @param melting_threshold::NF = 275 (bounds = Positive,)

    "[OPTION] Time scale for snow runoff/leakage into soil moisture [s]"
    @param runoff_time_scale::Second = Year(4) (bounds = Positive,)
end

Adapt.@adapt_structure SnowModel

# generator function
SnowModel(SG::SpectralGrid, geometry::LandGeometryOrNothing = nothing; kwargs...) = SnowModel{SG.NF}(; kwargs...)

function variables(::SnowModel)
    return (
        PrognosticVariable(:snow_depth, Grid2D(), namespace = :land, units = "m", desc = "Snow depth in equivalent liquid water height"),
        PrognosticVariable(:soil_temperature, Land3D(), namespace = :land, units = "K", desc = "Soil temperature"),
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

    Δt = model.time_stepping.Δt_sec
    (; snow_depth) = vars.prognostic.land                   # in equivalent liquid water height [m]
    (; soil_temperature) = vars.prognostic.land
    (; mask) = model.land_sea_mask

    # Some thermodynamics needed by snow
    ρ_water = model.atmosphere.water_density                # water density [kg/m³]
    Lᵢ = model.atmosphere.latent_heat_fusion                  # latent heat of fusion
    cₛ = model.land.thermodynamics.heat_capacity_dry_soil
    z₁ = model.land.geometry.layer_thickness[1]
    (; melting_threshold) = snow
    r⁻¹ = inv(convert(eltype(snow_depth), Second(snow.runoff_time_scale).value))

    # reset in any case
    vars.scratch.grid.a_2D .= 0
    
    # from precipitation schemes [m/s]
    snow_fall_rate = haskey(vars.parameterizations, :snow_rate) ?
        vars.parameterizations.snow_rate :
        vars.scratch.grid.a_2D               
    
    snow_melt_rate = vars.parameterizations.land.snow_melt_rate     # for soil moisture model

    params = (; melting_threshold, cₛ, z₁, Δt, ρ_water, Lᵢ, r⁻¹)

    launch!(
        architecture(snow_depth), LinearWorkOrder, size(snow_depth), land_snow_kernel!,
        snow_depth, soil_temperature, snow_melt_rate, snow_fall_rate, mask,
        params,
    )
    return nothing
end

@kernel inbounds = true function land_snow_kernel!(
        snow_depth, soil_temperature, snow_melt_rate, snow_fall_rate, mask,
        params,
    )
    ij = @index(Global, Linear)             # every grid point ij

    if mask[ij] > 0                         # at least partially land

        (; melting_threshold, cₛ, z₁, Δt, ρ_water, Lᵢ, r⁻¹) = params

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

        # Term 3: runoff rate [m/s], a leakage term to prevent too much snow, leaks into soil moisture
        # maximum allowed by existing snow depth
        runoff_rate_max = r⁻¹ * snow_depth[ij]

        # Adding the terms change snow depth by falling snow minus melting and runoff [m/s]
        # maximum amount of snow change
        dsnow_max = snow_fall_rate_max - melt_rate_max - runoff_rate_max

        # don't melt or runoff more than there is snow + snow that is falling down this time step
        # limited amount of snow change by how much there is
        dsnow = -min(snow_depth[ij] / Δt, -dsnow_max)

        # snow that we tried to melt/runoff that isn't available though
        dsnow_excess = dsnow_max - dsnow

        # store to pass to soil moisture [kg/m²/s], combined runoff with melt rate
        # limited to what's available to melt/runoff
        snow_melt_rate[ij] = (melt_rate_max + runoff_rate_max + dsnow_excess) * ρ_water

        # Euler forward time step but cap at 0 depth to not melt more snow than available
        snow_depth[ij] = max(snow_depth[ij] + Δt * dsnow, 0)
    end
end
