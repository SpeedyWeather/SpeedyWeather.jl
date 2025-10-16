abstract type AbstractSurfaceHumidityFlux <: AbstractParameterization end

export SurfaceHumidityFlux

"""Composite surface humidity flux type that holds flux
for ocean and land (each of any type) separately. Fields are
$(TYPEDFIELDS)"""
@kwdef struct SurfaceHumidityFlux{Ocean, Land} <: AbstractSurfaceHumidityFlux
    ocean::Ocean = SurfaceOceanHumidityFlux()
    land::Land = SurfaceLandHumidityFlux()
end

# generator function and defaults
function SurfaceHumidityFlux(
    SG::SpectralGrid; 
    ocean = SurfaceOceanHumidityFlux(SG),
    land = SurfaceLandHumidityFlux(SG))
    return SurfaceHumidityFlux(; ocean, land)
end

# initialize both
function initialize!(S::SurfaceHumidityFlux, model::PrimitiveWet)
    initialize!(S.ocean, model)
    initialize!(S.land, model)
end

# call first ocean flux then land flux
function parameterization!(ij, diagn, progn, humidity_flux::SurfaceHumidityFlux, model)
    surface_humidity_flux!(ij, diagn, progn, humidity_flux.ocean, model)
    surface_humidity_flux!(ij, diagn, progn, humidity_flux.land,  model)
end

export SurfaceOceanHumidityFlux
@kwdef struct SurfaceOceanHumidityFlux{NF<:AbstractFloat} <: AbstractSurfaceHumidityFlux
    "[OPTION] Use drag coefficient from calculated following model.boundary_layer_drag"
    use_boundary_layer_drag::Bool = true

    "[OPTION] Or fixed drag coefficient for humidity flux over ocean"
    drag::NF = 0.9e-3
end

SurfaceOceanHumidityFlux(SG::SpectralGrid; kwargs...) = SurfaceOceanHumidityFlux{SG.NF}(; kwargs...)
initialize!(::SurfaceOceanHumidityFlux, ::PrimitiveWet) = nothing

function surface_humidity_flux!(ij, diagn, progn, humidity_flux::SurfaceOceanHumidityFlux, model)

    # TODO use a skin temperature?
    T = progn.ocean.sea_surface_temperature[ij]

    # SATURATION HUMIDITY OVER OCEAN
    pₛ = diagn.grid.pres_grid_prev[ij]          # surface pressure [Pa]
    sat_humid_ocean = saturation_humidity(T, pₛ, model.clausius_clapeyron)

    ρ = diagn.physics.surface_air_density[ij]
    V₀ = diagn.physics.surface_wind_speed[ij]
    land_fraction = model.land_sea_mask[ij]
    surface_humid = diagn.grid.humid_grid[ij, end]

    # drag coefficient either from SurfaceHumidityFlux or from a central drag coefficient
    d = diagn.physics.boundary_layer_drag[ij]
    drag_ocean = humidity_flux.use_boundary_layer_drag ? d : humidity_flux.drag

    # SPEEDY documentation eq. 55/57, zero flux if sea surface temperature not available
    # but remove the max( ,0) to allow for surface condensation
    flux_ocean = isfinite(T) ? ρ*drag_ocean*V₀*(sat_humid_ocean  - surface_humid) : zero(T)

    # store without weighting by land fraction for coupling
    diagn.physics.ocean.surface_humidity_flux[ij] = flux_ocean         
    flux_ocean *= (1 - land_fraction)             # weight by ocean fraction of land-sea mask

    # output/diagnose: ocean sets flux (=), land accumulates (+=)
    diagn.physics.surface_humidity_flux[ij] = flux_ocean
    
    # accumulate with += into end=lowermost layer total flux
    diagn.tendencies.humid_tend_grid[ij, end] += surface_flux_to_tendency(flux_ocean, pₛ, model)
    return nothing
end

export SurfaceLandHumidityFlux
@kwdef struct SurfaceLandHumidityFlux{NF<:AbstractFloat} <: AbstractSurfaceHumidityFlux
    "[OPTION] Use column.boundary_layer_drag coefficient"
    use_boundary_layer_drag::Bool = true

    "[OPTION] Or fixed drag coefficient for humidity flux over land"
    drag::NF = 1.2e-3
end
    
SurfaceLandHumidityFlux(SG::SpectralGrid; kwargs...) = SurfaceLandHumidityFlux{SG.NF}(; kwargs...)
initialize!(::SurfaceLandHumidityFlux, ::PrimitiveWet) = nothing

function surface_humidity_flux!(ij, diagn, progn, humidity_flux::SurfaceLandHumidityFlux, model)

    # TODO use a skin temperature?
    T = progn.land.soil_temperature[ij, 1]  # uppermost land layer with index 1
    α = diagn.physics.land.soil_moisture_availability[ij]

    # SATURATION HUMIDITY OVER LAND
    pₛ = diagn.grid.pres_grid_prev[ij]          # surface pressure [Pa]
    sat_humid_land = saturation_humidity(T, pₛ, model.clausius_clapeyron)

    ρ = diagn.physics.surface_air_density[ij]
    V₀ = diagn.physics.surface_wind_speed[ij]
    land_fraction = model.land_sea_mask[ij]
    surface_humid = diagn.grid.humid_grid[ij, end]

    # drag coefficient either from SurfaceLandHumidityFlux or from a central drag coefficient
    d = diagn.physics.boundary_layer_drag[ij]
    drag_land = humidity_flux.use_boundary_layer_drag ? d : humidity_flux.drag

    # SPEEDY documentation eq. 55/57, zero flux if land / soil moisture availability not available (=ocean)
    # but remove the max( ,0) to allow for surface condensation
    flux_land = isfinite(T) && isfinite(α) ? ρ*drag_land*V₀*(α*sat_humid_land  - surface_humid) : zero(T)
    
    # store without weighting by land fraction for coupling
    diagn.physics.land.surface_humidity_flux[ij] = flux_land         
    flux_land *= land_fraction             # weight by land fraction of land-sea mask

    # output/diagnose: ocean sets flux (=), land accumulates (+=)
    diagn.physics.surface_humidity_flux[ij] += flux_land
    
    # accumulate with += into end=lowermost layer total flux
    diagn.tendencies.humid_tend_grid[ij, end] += surface_flux_to_tendency(flux_land, pₛ, model)
    return nothing
end

## ----

export PrescribedOceanHumidityFlux
struct PrescribedOceanHumidityFlux <: AbstractSurfaceHumidityFlux end
PrescribedOceanHumidityFlux(::SpectralGrid) = PrescribedOceanHumidityFlux()
initialize!(::PrescribedOceanHumidityFlux, ::PrimitiveWet) = nothing

function surface_humidity_flux!(ij, diagn, progn, ::PrescribedOceanHumidityFlux, model)
    land_fraction = model.land_sea_mask[ij]
    pₛ = diagn.grid.pres_grid_prev[ij]          # surface pressure [Pa]

    # read in a prescribed flux
    flux_ocean = progn.ocean.surface_humidity_flux[ij]

    # store without weighting by land fraction for coupling
    diagn.physics.ocean.surface_humidity_flux[ij] = flux_ocean         
    flux_ocean *= (1 - land_fraction)             # weight by ocean fraction of land-sea mask

    # output/diagnose: ocean sets flux (=), land accumulates (+=)
    diagn.physics.surface_humidity_flux[ij] = flux_ocean
    
    # accumulate with += into end=lowermost layer total flux
    diagn.tendencies.humid_tend_grid[ij, end] += surface_flux_to_tendency(flux_ocean, pₛ, model)
    return nothing
end

## ----

export PrescribedLandHumidityFlux
struct PrescribedLandHumidityFlux <: AbstractSurfaceHumidityFlux end
PrescribedLandHumidityFlux(::SpectralGrid) = PrescribedLandHumidityFlux()
initialize!(::PrescribedLandHumidityFlux, ::PrimitiveWet) = nothing

function surface_humidity_flux!(ij, diagn, progn, ::PrescribedLandHumidityFlux, model)
    land_fraction = model.land_sea_mask[ij]
    pₛ = diagn.grid.pres_grid_prev[ij]          # surface pressure [Pa]

    # read in a prescribed flux
    flux_land = progn.land.surface_humidity_flux[ij]

    # store without weighting by land fraction for coupling
    diagn.physics.land.surface_humidity_flux[ij] = flux_land       
    flux_land *= land_fraction                  # weight by ocean fraction of land-sea mask

    # output/diagnose: ocean sets flux (=), land accumulates (+=)
    diagn.physics.surface_humidity_flux[ij] += flux_land
    
    # accumulate with += into end=lowermost layer total flux
    diagn.tendencies.humid_tend_grid[ij, end] += surface_flux_to_tendency(flux_land, pₛ, model)
    return nothing
end