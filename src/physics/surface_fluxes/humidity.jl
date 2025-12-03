abstract type AbstractSurfaceHumidityFlux <: AbstractParameterization end

export SurfaceHumidityFlux

"""Composite surface humidity flux type that holds flux
for ocean and land (each of any type) separately. Fields are
$(TYPEDFIELDS)"""
@kwdef struct SurfaceHumidityFlux{Ocean, Land} <: AbstractSurfaceHumidityFlux
    "[OPTION] Surface humidity flux parameterization for ocean surfaces"
    ocean::Ocean

    "[OPTION] Surface humidity flux parameterization for land surfaces"
    land::Land
end

Adapt.@adapt_structure SurfaceHumidityFlux

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
# order important as ocean sets the flux (=) and land accumulates (+=)
@propagate_inbounds function parameterization!(ij, diagn, progn, humidity_flux::SurfaceHumidityFlux, model)
    surface_humidity_flux!(ij, diagn, progn, humidity_flux.ocean, model)
    surface_humidity_flux!(ij, diagn, progn, humidity_flux.land,  model)
end

export SurfaceOceanHumidityFlux
"""Humidity flux parameterization over ocean surfaces. Fields are $(TYPEDFIELDS)"""
@kwdef struct SurfaceOceanHumidityFlux{NF} <: AbstractSurfaceHumidityFlux
    "[OPTION] Use drag coefficient from calculated following model.boundary_layer_drag"
    use_boundary_layer_drag::Bool = true

    "[OPTION] Or fixed drag coefficient for humidity flux over ocean"
    drag::NF = 0.9e-3
end

Adapt.@adapt_structure SurfaceOceanHumidityFlux
SurfaceOceanHumidityFlux(SG::SpectralGrid; kwargs...) = SurfaceOceanHumidityFlux{SG.NF}(; kwargs...)
initialize!(::SurfaceOceanHumidityFlux, ::PrimitiveWet) = nothing

@propagate_inbounds function surface_humidity_flux!(ij, diagn, progn, humidity_flux::SurfaceOceanHumidityFlux, model)

    surface = model.geometry.nlayers
    T = progn.ocean.sea_surface_temperature[ij]

    # SATURATION HUMIDITY OVER OCEAN
    pₛ = diagn.grid.pres_grid_prev[ij]          # surface pressure [Pa]
    sat_humid_ocean = saturation_humidity(T, pₛ, model.clausius_clapeyron)

    ρ = diagn.physics.surface_air_density[ij]
    V₀ = diagn.physics.surface_wind_speed[ij]
    land_fraction = model.land_sea_mask.mask[ij]
    surface_humid = diagn.grid.humid_grid_prev[ij, surface]

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
    diagn.tendencies.humid_tend_grid[ij, surface] += surface_flux_to_tendency(flux_ocean, pₛ, model)
    return nothing
end

export SurfaceLandHumidityFlux
"""Humidity flux parameterization over land surfaces. Fields are $(TYPEDFIELDS)"""
@kwdef struct SurfaceLandHumidityFlux{NF} <: AbstractSurfaceHumidityFlux
    "[OPTION] Use column.boundary_layer_drag coefficient"
    use_boundary_layer_drag::Bool = true

    "[OPTION] Or fixed drag coefficient for humidity flux over land"
    drag::NF = 1.2e-3
end

Adapt.@adapt_structure SurfaceLandHumidityFlux
SurfaceLandHumidityFlux(SG::SpectralGrid; kwargs...) = SurfaceLandHumidityFlux{SG.NF}(; kwargs...)
initialize!(::SurfaceLandHumidityFlux, ::PrimitiveWet) = nothing

@propagate_inbounds function surface_humidity_flux!(ij, diagn, progn, humidity_flux::SurfaceLandHumidityFlux, model)

    # TODO use a skin temperature?
    T = progn.land.soil_temperature[ij, 1]  # uppermost land layer with index 1
    α = diagn.physics.land.soil_moisture_availability[ij]

    # SATURATION HUMIDITY OVER LAND
    pₛ = diagn.grid.pres_grid_prev[ij]          # surface pressure [Pa]
    sat_humid_land = saturation_humidity(T, pₛ, model.clausius_clapeyron)

    ρ = diagn.physics.surface_air_density[ij]
    V₀ = diagn.physics.surface_wind_speed[ij]
    land_fraction = model.land_sea_mask.mask[ij]
    surface = model.geometry.nlayers            # indexing top to bottom
    surface_humid = diagn.grid.humid_grid_prev[ij, surface]

    # drag coefficient either from SurfaceLandHumidityFlux or from a central drag coefficient
    d = diagn.physics.boundary_layer_drag[ij]
    drag_land = humidity_flux.use_boundary_layer_drag ? d : humidity_flux.drag

    # SPEEDY documentation eq. 55/57, zero flux if land / soil moisture availability not available (=ocean)
    # but remove the max( ,0) to allow for surface condensation
    flux_land = isfinite(T) && isfinite(α) ? ρ*drag_land*V₀*(α*sat_humid_land - surface_humid) : zero(T)
    
    # store without weighting by land fraction for coupling
    # TODO multiply by Lᵥ for W/m^2 here?
    diagn.physics.land.surface_humidity_flux[ij] = flux_land    # [kg/m²/s]
    flux_land *= land_fraction             # weight by land fraction of land-sea mask

    # output/diagnose: ocean sets flux (=), land accumulates (+=)
    diagn.physics.surface_humidity_flux[ij] += flux_land
    
    # accumulate with += into end=lowermost layer total flux
    diagn.tendencies.humid_tend_grid[ij, surface] += surface_flux_to_tendency(flux_land, pₛ, model)
    return nothing
end

## ----

export PrescribedOceanHumidityFlux
"""Prescribed humidity flux over ocean surfaces. Applies surface humidity flux from
`progn.ocean.surface_humidity_flux`. $(TYPEDFIELDS)"""
struct PrescribedOceanHumidityFlux <: AbstractSurfaceHumidityFlux end
Adapt.@adapt_structure PrescribedOceanHumidityFlux
PrescribedOceanHumidityFlux(::SpectralGrid) = PrescribedOceanHumidityFlux()
initialize!(::PrescribedOceanHumidityFlux, ::PrimitiveWet) = nothing

@propagate_inbounds function surface_humidity_flux!(ij, diagn, progn, ::PrescribedOceanHumidityFlux, model)
    land_fraction = model.land_sea_mask.mask[ij]
    pₛ = diagn.grid.pres_grid_prev[ij]          # surface pressure [Pa]
    surface = model.geometry.nlayers

    # read in a prescribed flux
    flux_ocean = progn.ocean.surface_humidity_flux[ij]

    # store without weighting by land fraction for coupling
    diagn.physics.ocean.surface_humidity_flux[ij] = flux_ocean         
    flux_ocean *= (1 - land_fraction)             # weight by ocean fraction of land-sea mask

    # output/diagnose: ocean sets flux (=), land accumulates (+=)
    diagn.physics.surface_humidity_flux[ij] = flux_ocean
    
    # accumulate with += into end=lowermost layer total flux
    diagn.tendencies.humid_tend_grid[ij, surface] += surface_flux_to_tendency(flux_ocean, pₛ, model)
    return nothing
end


## ----

export PrescribedLandHumidityFlux
"""Prescribed humidity flux over land surfaces. Applies surface humidity flux from
`progn.land.surface_humidity_flux`. $(TYPEDFIELDS)"""
struct PrescribedLandHumidityFlux <: AbstractSurfaceHumidityFlux end
Adapt.@adapt_structure PrescribedLandHumidityFlux
PrescribedLandHumidityFlux(::SpectralGrid) = PrescribedLandHumidityFlux()
initialize!(::PrescribedLandHumidityFlux, ::PrimitiveWet) = nothing

@propagate_inbounds function surface_humidity_flux!(ij, diagn, progn, ::PrescribedLandHumidityFlux, model)
    land_fraction = model.land_sea_mask.mask[ij]
    pₛ = diagn.grid.pres_grid_prev[ij]          # surface pressure [Pa]
    surface = model.geometry.nlayers

    # read in a prescribed flux
    flux_land = progn.land.surface_humidity_flux[ij]

    # store without weighting by land fraction for coupling
    diagn.physics.land.surface_humidity_flux[ij] = flux_land       
    flux_land *= land_fraction                  # weight by ocean fraction of land-sea mask

    # output/diagnose: ocean sets flux (=), land accumulates (+=)
    diagn.physics.surface_humidity_flux[ij] += flux_land
    
    # accumulate with += into end=lowermost layer total flux
    diagn.tendencies.humid_tend_grid[ij, surface] += surface_flux_to_tendency(flux_land, pₛ, model)
    return nothing
end

function variables(::AbstractSurfaceHumidityFlux)
    return (
        DiagnosticVariable(name=:surface_humidity_flux, dims=Grid2D(), desc="Total surface humidity flux", units="kg/m²/s"),
        DiagnosticVariable(name=:surface_humidity_flux, dims=Grid2D(), desc="Ocean surface humidity flux", units="kg/m²/s", namespace=:ocean),
        DiagnosticVariable(name=:surface_humidity_flux, dims=Grid2D(), desc="Land surface humidity flux", units="kg/m²/s", namespace=:land),
        PrognosticVariable(name=:surface_humidity_flux, dims=Grid2D(), desc="Prescribed Ocean surface humidity flux", units="kg/m²/s", namespace=:ocean),
        PrognosticVariable(name=:surface_humidity_flux, dims=Grid2D(), desc="Prescribed Land surface humidity flux", units="kg/m²/s", namespace=:land),
    )
end
