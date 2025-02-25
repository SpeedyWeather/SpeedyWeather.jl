export SurfaceEvaporation
@kwdef struct SurfaceEvaporation{Ocean, Land} <: AbstractSurfaceEvaporation
    ocean::Ocean = SurfaceOceanEvaporation()
    land::Land = SurfaceLandEvaporation()
end

function SurfaceEvaporation(
    SG::SpectralGrid; 
    ocean = SurfaceOceanEvaporation(SG),
    land = SurfaceLandEvaporation(SG))
    return SurfaceEvaporation(; ocean, land)
end

function initialize!(S::SurfaceEvaporation, model::PrimitiveWet)
    initialize!(S.ocean, model)
    initialize!(S.land, model)
end

function surface_evaporation!(   
    column::ColumnVariables,
    evaporation::SurfaceEvaporation,
    progn::PrognosticVariables,
    model::PrimitiveWet,
)   
    surface_evaporation!(column, evaporation.ocean, progn, model)
    surface_evaporation!(column, evaporation.land, progn, model)
end

## ----

export NoSurfaceEvaporation
struct NoSurfaceEvaporation <: AbstractSurfaceEvaporation end
NoSurfaceEvaporation(::SpectralGrid) = NoSurfaceEvaporation()
initialize!(::NoSurfaceEvaporation, ::PrimitiveWet) = nothing
surface_evaporation!(::ColumnVariables, ::NoSurfaceEvaporation, ::PrimitiveWet) = nothing

## ----

export SurfaceOceanEvaporation
@kwdef struct SurfaceOceanEvaporation{NF<:AbstractFloat} <: AbstractSurfaceEvaporation
    "Use column.boundary_layer_drag coefficient"
    use_boundary_layer_drag::Bool = true

    "Or drag coefficient for evaporation over ocean"
    moisture_exchange::NF = 0.9e-3
end

SurfaceOceanEvaporation(SG::SpectralGrid; kwargs...) = SurfaceOceanEvaporation{SG.NF}(; kwargs...)
initialize!(::SurfaceOceanEvaporation, ::PrimitiveWet) = nothing

function surface_evaporation!(
    column::ColumnVariables{NF},
    evaporation::SurfaceOceanEvaporation,
    model::PrimitiveWet,
) where NF

    (; skin_temperature_sea, pres) = column
    (; moisture_exchange) = evaporation

    # SATURATION HUMIDITY OVER OCEAN
    surface_pressure = pres[end]
    sat_humid_sea = saturation_humidity(skin_temperature_sea, surface_pressure, model.clausius_clapeyron)

    ρ = column.surface_air_density
    V₀ = column.surface_wind_speed
    land_fraction = column.land_fraction
    (; surface_humid) = column

    # drag coefficient either from SurfaceEvaporation or from a central drag coefficient
    drag_sea = evaporation.use_boundary_layer_drag ? column.boundary_layer_drag : moisture_exchange

    # SPEEDY documentation eq. 55/57, zero flux if sea surface temperature not available
    flux_sea = isfinite(skin_temperature_sea) ? ρ*drag_sea*V₀*max(sat_humid_sea  - surface_humid, zero(NF)) :
                                                    zero(skin_temperature_sea)
    column.evaporative_flux_ocean = flux_sea

    flux_sea *= (1 - land_fraction)             # weight by ocean fraction of land-sea mask
    column.flux_humid_upward[end] += flux_sea   # accumulate with += into end=lowermost layer total flux
    column.evaporative_flux = flux_sea          # output/diagnose: ocean sets flux (=), land accumulates (+=)
    return nothing
end

## ----

export SurfaceLandEvaporation
@kwdef struct SurfaceLandEvaporation{NF<:AbstractFloat} <: AbstractSurfaceEvaporation
    "Use column.boundary_layer_drag coefficient"
    use_boundary_layer_drag::Bool = true

    "Otherwise, use the following drag coefficient for evaporation over land"
    moisture_exchange::NF = 1.2e-3
end

SurfaceLandEvaporation(SG::SpectralGrid; kwargs...) = SurfaceLandEvaporation{SG.NF}(; kwargs...)
initialize!(::SurfaceLandEvaporation, ::PrimitiveWet) = nothing

function surface_evaporation!(
    column::ColumnVariables{NF},
    evaporation::SurfaceLandEvaporation,
    model::PrimitiveWet,
) where NF

    (; skin_temperature_land, pres) = column
    (; moisture_exchange) = evaporation
    α = column.soil_moisture_availability

    # SATURATION HUMIDITY OVER LAND
    surface_pressure = pres[end]
    sat_humid_land = saturation_humidity(skin_temperature_land, surface_pressure, model.clausius_clapeyron)

    ρ = column.surface_air_density
    V₀ = column.surface_wind_speed
    land_fraction = column.land_fraction
    (; surface_humid) = column

    # drag coefficient either from SurfaceLandEvaporation or from a central drag coefficient
    drag_land = evaporation.use_boundary_layer_drag ? column.boundary_layer_drag : moisture_exchange

    # SPEEDY documentation eq. 55/57, zero flux if land / soil moisture availability not available (=ocean)
    flux_land = isfinite(skin_temperature_land) && isfinite(α) ?
                ρ*drag_land*V₀*max(α*sat_humid_land  - surface_humid, zero(NF))*land_fraction :
                zero(NF)
    column.evaporative_flux_land = flux_land        # store flux separately for land
    flux_land *= land_fraction                      # weight by land fraction of land-sea mask
    column.flux_humid_upward[end] += flux_land      # end=lowermost layer, accumulate with (+=) to total flux
    column.evaporative_flux += flux_land            # ocean sets the flux (=), land accumulates (+=)
    return nothing
end

## ----

export PrescribedOceanEvaporation
struct PrescribedOceanEvaporation <: AbstractSurfaceEvaporation end
PrescribedOceanEvaporation(::SpectralGrid) = PrescribedOceanEvaporation()
initialize!(::PrescribedOceanEvaporation, ::PrimitiveWet) = nothing
function surface_evaporation!(
    column::ColumnVariables,
    ::PrescribedOceanEvaporation,
    progn::PrognosticVariables,
    model::PrimitiveWet)

    land_fraction = column.land_fraction

    # read in a prescribed flux
    flux = progn.ocean.evaporative_flux[column.ij]
    column.evaporative_flux_ocean = flux    # store ocean-only flux separately too

    flux *= (1-land_fraction)               # weight by ocean fraction of land-sea mask
    column.flux_humid_upward[end] += flux   # end=lowermost layer, accumulate with (+=) to total flux
    column.evaporative_flux = flux          # ocean sets the flux (=), land accumulates (+=)
end

## ----

export PrescribedLandEvaporation
struct PrescribedLandEvaporation <: AbstractSurfaceEvaporation end
PrescribedLandEvaporation(::SpectralGrid) = PrescribedLandEvaporation()
initialize!(::PrescribedLandEvaporation, ::PrimitiveWet) = nothing
function surface_evaporation!(
    column::ColumnVariables,
    ::PrescribedLandEvaporation,
    progn::PrognosticVariables,
    model::PrimitiveWet)

    land_fraction = column.land_fraction

    # read in a prescribed flux
    flux = progn.land.evaporative_flux[column.ij]
    column.evaporative_flux_land = flux     # store land-only flux separately

    flux *= land_fraction
    column.flux_humid_upward[end] += flux   # end=lowermost layer
    column.evaporative_flux += flux         # ocean sets the flux (=), land accumulates (+=)
end