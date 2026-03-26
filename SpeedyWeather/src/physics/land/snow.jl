abstract type AbstractSnow <: AbstractLandComponent end

export SnowModel    # maybe change for a more concise name later

"""
    SnowModel(; melting_threshold=275, snow_depth_cap=10)

Single-column snow bucket model in equivalent liquid water depth. Snow accumulates
from the diagnosed precipitation, melts once the top soil layer exceeds
`melting_threshold`, and is capped at `snow_depth_cap` to limit infinite snow/ice accumulation
over perennial ice caps and glaciers.
$(TYPEDFIELDS)"""
@parameterized @kwdef struct SnowModel{NF} <: AbstractSnow
    "[OPTION] Temperature threshold for snow melting [K]"
    @param melting_threshold::NF = 275 (bounds = Positive,)

    "[OPTION] Permanent snow/ice depth cap in equivalent liquid water depth [m]"
    snow_depth_cap::NF = 10
end

Adapt.@adapt_structure SnowModel

# generator function
SnowModel(SG::SpectralGrid, geometry::LandGeometryOrNothing = nothing; kwargs...) = SnowModel{SG.NF}(; kwargs...)

# initialize component
initialize!(snow::SnowModel, model::PrimitiveEquation) = nothing

# set initial conditions for snow depth in initial conditions
function initialize!(
        progn::PrognosticVariables,
        diagn::DiagnosticVariables,
        snow::SnowModel,
        model::PrimitiveEquation,
    )
    return set!(progn, model.geometry, snow_depth = 0)
end

function timestep!(
        progn::PrognosticVariables,
        diagn::DiagnosticVariables,
        snow::SnowModel,
        model::PrimitiveEquation,
    )
    Δt = model.time_stepping.Δt_sec
    (; snow_depth) = progn.land                             # in equivalent liquid water height [m]
    (; soil_temperature) = progn.land
    (; mask) = model.land_sea_mask

    # Some thermodynamics needed by snow
    ρ_water = model.atmosphere.water_density                # water density [kg/m³]
    Lᵢ = model.atmosphere.latent_heat_fusion                  # latent heat of fusion
    cₛ = model.land.thermodynamics.heat_capacity_dry_soil
    z₁ = model.land.geometry.layer_thickness[1]
    (; melting_threshold, snow_depth_cap) = snow

    snow_fall_rate = diagn.physics.snow_rate                # from precipitation schemes [m/s]
    snow_melt_rate = diagn.physics.land.snow_melt_rate      # for soil moisture model

    params = (; melting_threshold, cₛ, z₁, Δt, ρ_water, Lᵢ, snow_depth_cap)

    return launch!(
        architecture(snow_depth), LinearWorkOrder, size(snow_depth), land_snow_kernel!,
        snow_depth, soil_temperature, snow_melt_rate, snow_fall_rate, mask,
        params,
    )
end

@kernel inbounds = true function land_snow_kernel!(
        snow_depth, soil_temperature, snow_melt_rate, snow_fall_rate, mask,
        params,
    )
    ij = @index(Global, Linear)             # every grid point ij

    if mask[ij] > 0                         # at least partially land

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
        snow_depth[ij] = min(snow_depth_forward, snow_depth_cap)  # cap at permanent snow/ice depth
    end
end

function variables(::SnowModel)
    return (
        PrognosticVariable(name = :snow_depth, dims = Grid2D(), namespace = :land, units = "m", desc = "Snow depth in equivalent liquid water height"),
        PrognosticVariable(name = :soil_temperature, dims = Grid3D(), namespace = :land, units = "K", desc = "Soil temperature"),
        DiagnosticVariable(name = :snow_melt_rate, dims = Grid2D(), namespace = :land, units = "kg/m²/s", desc = "Snow melt rate"),
    )
end
