export SurfaceHumidityFlux
@kwdef struct SurfaceHumidityFlux{Ocean, Land} <: AbstractSurfaceHumidityFlux
    ocean::Ocean = SurfaceOceanHumidityFlux()
    land::Land = SurfaceLandHumidityFlux()
end

function SurfaceHumidityFlux(
    SG::SpectralGrid; 
    ocean = SurfaceOceanHumidityFlux(SG),
    land = SurfaceLandHumidityFlux(SG))
    return SurfaceHumidityFlux(; ocean, land)
end

function initialize!(S::SurfaceHumidityFlux, model::PrimitiveWet)
    initialize!(S.ocean, model)
    initialize!(S.land, model)
end

function surface_humidity_flux!(   
    column::ColumnVariables,
    humidity_flux::SurfaceHumidityFlux,
    progn::PrognosticVariables,
    model::PrimitiveWet,
)   
    surface_humidity_flux!(column, humidity_flux.ocean, progn, model)
    surface_humidity_flux!(column, humidity_flux.land, progn, model)
end

## ----

export NoSurfaceHumidityFlux
struct NoSurfaceHumidityFlux <: AbstractSurfaceHumidityFlux end
NoSurfaceHumidityFlux(::SpectralGrid) = NoSurfaceHumidityFlux()
initialize!(::NoSurfaceHumidityFlux, ::PrimitiveWet) = nothing
surface_humidity_flux!(::ColumnVariables, ::NoSurfaceHumidityFlux, ::PrimitiveWet) = nothing

## ----

export SurfaceOceanHumidityFlux
@kwdef struct SurfaceOceanHumidityFlux{NF<:AbstractFloat} <: AbstractSurfaceHumidityFlux
    "[OPTION] Use column.boundary_layer_drag coefficient"
    use_boundary_layer_drag::Bool = true

    "[OPTION] Or fixed drag coefficient for humidity flux over ocean"
    drag::NF = 0.9e-3
end

SurfaceOceanHumidityFlux(SG::SpectralGrid; kwargs...) = SurfaceOceanHumidityFlux{SG.NF}(; kwargs...)
initialize!(::SurfaceOceanHumidityFlux, ::PrimitiveWet) = nothing

function surface_humidity_flux!(
    column::ColumnVariables{NF},
    humidity_flux::SurfaceOceanHumidityFlux,
    model::PrimitiveWet,
) where NF

    (; skin_temperature_sea, pres) = column
    (; drag) = humidity_flux

    # SATURATION HUMIDITY OVER OCEAN
    surface_pressure = pres[end]
    sat_humid_sea = saturation_humidity(skin_temperature_sea, surface_pressure, model.clausius_clapeyron)

    ρ = column.surface_air_density
    V₀ = column.surface_wind_speed
    land_fraction = column.land_fraction
    (; surface_humid) = column

    # drag coefficient either from SurfaceHumidityFlux or from a central drag coefficient
    drag_sea = humidity_flux.use_boundary_layer_drag ? column.boundary_layer_drag : drag

    # SPEEDY documentation eq. 55/57, zero flux if sea surface temperature not available
    # but remove the max( ,0) to allow for surface condensation
    flux_sea = isfinite(skin_temperature_sea) ? ρ*drag_sea*V₀*(sat_humid_sea  - surface_humid) :
                                                    zero(skin_temperature_sea)
    column.humidity_flux_ocean = flux_sea       # store without weighting by land fraction for coupling

    flux_sea *= (1 - land_fraction)             # weight by ocean fraction of land-sea mask
    column.flux_humid_upward[end] += flux_sea   # accumulate with += into end=lowermost layer total flux
    column.humidity_flux = flux_sea             # output/diagnose: ocean sets flux (=), land accumulates (+=)
    return nothing
end

## ----

export SurfaceLandHumidityFlux
@kwdef struct SurfaceLandHumidityFlux{NF<:AbstractFloat} <: AbstractSurfaceHumidityFlux
    "[OPTION] Use column.boundary_layer_drag coefficient"
    use_boundary_layer_drag::Bool = true

    "[OPTION] Otherwise, use the following drag coefficient for humidity flux (evaporation) over land"
    drag::NF = 1.2e-3
end
    
SurfaceLandHumidityFlux(SG::SpectralGrid; kwargs...) = SurfaceLandHumidityFlux{SG.NF}(; kwargs...)
initialize!(::SurfaceLandHumidityFlux, ::PrimitiveWet) = nothing

function surface_humidity_flux!(
    column::ColumnVariables{NF},
    humidity_flux::SurfaceLandHumidityFlux,
    model::PrimitiveWet,
) where NF

    (; skin_temperature_land, pres) = column
    (; drag) = humidity_flux
    α = column.soil_moisture_availability

    # SATURATION HUMIDITY OVER LAND
    surface_pressure = pres[end]
    sat_humid_land = saturation_humidity(skin_temperature_land, surface_pressure, model.clausius_clapeyron)

    ρ = column.surface_air_density
    V₀ = column.surface_wind_speed
    land_fraction = column.land_fraction
    (; surface_humid) = column

    # drag coefficient either from SurfaceLandHumidityFlux or from a central drag coefficient
    drag_land = humidity_flux.use_boundary_layer_drag ? column.boundary_layer_drag : drag

    # SPEEDY documentation eq. 55/57, zero flux if land / soil moisture availability not available (=ocean)
    # but remove the max( ,0) to allow for surface condensation
    flux_land = isfinite(skin_temperature_land) && isfinite(α) ?
                ρ*drag_land*V₀*(α*sat_humid_land  - surface_humid) : zero(NF)
    column.humidity_flux_land = flux_land           # store flux separately for land
    flux_land *= land_fraction                      # weight by land fraction of land-sea mask
    column.flux_humid_upward[end] += flux_land      # end=lowermost layer, accumulate with (+=) to total flux
    column.humidity_flux += flux_land               # ocean sets the flux (=), land accumulates (+=)
    return nothing
end

## ----

export PrescribedOceanHumidityFlux
struct PrescribedOceanHumidityFlux <: AbstractSurfaceHumidityFlux end
PrescribedOceanHumidityFlux(::SpectralGrid) = PrescribedOceanHumidityFlux()
initialize!(::PrescribedOceanHumidityFlux, ::PrimitiveWet) = nothing

function surface_humidity_flux!(
    column::ColumnVariables,
    ::PrescribedOceanHumidityFlux,
    progn::PrognosticVariables,
    model::PrimitiveWet)

    land_fraction = column.land_fraction

    # read in a prescribed flux
    flux = progn.ocean.humidity_flux[column.ij]
    column.humidity_flux_ocean = flux       # store ocean-only flux separately too

    flux *= (1-land_fraction)               # weight by ocean fraction of land-sea mask
    column.flux_humid_upward[end] += flux   # end=lowermost layer, accumulate with (+=) to total flux
    column.humidity_flux = flux             # ocean sets the flux (=), land accumulates (+=)
end

## ----

export PrescribedLandHumidityFlux
struct PrescribedLandHumidityFlux <: AbstractSurfaceHumidityFlux end
PrescribedLandHumidityFlux(::SpectralGrid) = PrescribedLandHumidityFlux()
initialize!(::PrescribedLandHumidityFlux, ::PrimitiveWet) = nothing

function surface_humidity_flux!(
    column::ColumnVariables,
    ::PrescribedLandHumidityFlux,
    progn::PrognosticVariables,
    model::PrimitiveWet)

    land_fraction = column.land_fraction

    # read in a prescribed flux
    flux = progn.land.humidity_flux[column.ij]
    column.humidity_flux_land = flux        # store land-only flux separately

    flux *= land_fraction
    column.flux_humid_upward[end] += flux   # end=lowermost layer
    column.humidity_flux += flux            # ocean sets the flux (=), land accumulates (+=)
end