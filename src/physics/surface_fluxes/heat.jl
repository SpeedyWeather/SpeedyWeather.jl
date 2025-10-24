abstract type AbstractSurfaceHeatFlux <: AbstractParameterization end

export SurfaceHeatFlux

"""Composite surface (sensible) heat flux type that holds flux
for ocean and land (each of any type) separately. Fields are
$(TYPEDFIELDS)"""
@kwdef struct SurfaceHeatFlux{Ocean, Land} <: AbstractSurfaceHeatFlux
    ocean::Ocean = SurfaceOceanHeatFlux()
    land::Land = SurfaceLandHeatFlux()
end

# generator function
function SurfaceHeatFlux(
    SG::SpectralGrid; 
    ocean = SurfaceOceanHeatFlux(SG),
    land = SurfaceLandHeatFlux(SG))
    return SurfaceHeatFlux(; ocean, land)
end

function initialize!(S::SurfaceHeatFlux, model::PrimitiveEquation)
    initialize!(S.ocean, model)
    initialize!(S.land, model)
end

# call first ocean flux then land flux
function parameterization!(ij, diagn, progn, heat_flux::SurfaceHeatFlux, model)
    surface_heat_flux!(ij, diagn, progn, heat_flux.ocean, model)
    surface_heat_flux!(ij, diagn, progn, heat_flux.land,  model)
end

Adapt.@adapt_structure SurfaceHeatFlux

export SurfaceOceanHeatFlux

"""Surface sensible heat flux parameterization over ocean. Calculates the 
turbulent exchange of sensible heat between the ocean surface and the atmosphere
based on temperature differences and wind speed. Uses bulk aerodynamic formulas
with drag coefficients. Fields are $(TYPEDFIELDS)"""
@kwdef struct SurfaceOceanHeatFlux{NF} <: AbstractSurfaceHeatFlux
    "[OPTION] Use drag coefficient from calculated following model.boundary_layer_drag"
    use_boundary_layer_drag::Bool = true

    "[OPTION] Or fixed drag coefficient for heat fluxes over ocean"
    drag::NF = 0.9e-3
end

SurfaceOceanHeatFlux(SG::SpectralGrid; kwargs...) = SurfaceOceanHeatFlux{SG.NF}(; kwargs...)
initialize!(::SurfaceOceanHeatFlux, ::PrimitiveEquation) = nothing

function surface_heat_flux!(ij, diagn, progn, heat_flux::SurfaceOceanHeatFlux, model)

    nlayers = model.geometry.nlayers
    cₚ = model.atmosphere.heat_capacity
    ρ = diagn.physics.surface_air_density[ij]
    V₀ = diagn.physics.surface_wind_speed[ij]

    # TODO actually implement skin temperature?
    T_skin_ocean = progn.ocean.sea_surface_temperature[ij]
    T = diagn.physics.surface_air_temperature[ij]
    land_fraction = model.land_sea_mask[ij]
    pₛ = diagn.grid.pres_grid_prev[ij]

    # drag coefficient
    d = diagn.physics.boundary_layer_drag[ij]
    drag_ocean = heat_flux.use_boundary_layer_drag ? d : heat_flux.drag

    # SPEEDY documentation Eq. 54/56, land/sea fraction included
    # Only flux from ocean if available (not NaN) otherwise zero flux
    # leave out *cₚ here but include below to avoid division
    flux_ocean  = isfinite(T_skin_ocean) ? ρ*drag_ocean*V₀*(T_skin_ocean  - T) : zero(T_skin_ocean)
    
    # store without weighting by land fraction for coupling
    diagn.physics.ocean.sensible_heat_flux[ij] = flux_ocean*cₚ  # to store ocean flux separately too
    
    flux_ocean *= (1-land_fraction)                             # weight by ocean fraction of land-sea mask
    
    # output/diagnostics: ocean sets the flux (=), land accumulates (+=)
    diagn.physics.sensible_heat_flux[ij] = flux_ocean*cₚ
    
    # accumulate with += into end=lowermost layer total flux
    diagn.tendencies.temp_tend_grid[ij, nlayers] += surface_flux_to_tendency(flux_ocean, pₛ, model)
    return nothing
end

Adapt.@adapt_structure SurfaceOceanHeatFlux

export SurfaceLandHeatFlux

"""Surface sensible heat flux parameterization over land. Calculates the 
turbulent exchange of sensible heat between the land surface and the atmosphere
based on soil temperature differences and wind speed. Uses bulk aerodynamic 
formulas with drag coefficients that can account for surface roughness. 
Fields are $(TYPEDFIELDS)"""
@kwdef struct SurfaceLandHeatFlux{NF} <: AbstractSurfaceHeatFlux
    "[OPTION] Use drag coefficient from calculated following model.boundary_layer_drag"
    use_boundary_layer_drag::Bool = true

    "[OPTION] Or fixed drag coefficient for heat fluxes over land"
    drag::NF = 1.2e-3       # for neutral stability
end

SurfaceLandHeatFlux(SG::SpectralGrid; kwargs...) = SurfaceLandHeatFlux{SG.NF}(; kwargs...)
initialize!(::SurfaceLandHeatFlux, ::PrimitiveEquation) = nothing

function surface_heat_flux!(ij, diagn, progn, heat_flux::SurfaceLandHeatFlux, model)
    
    cₚ = model.atmosphere.heat_capacity
    pₛ = diagn.grid.pres_grid_prev[ij]                  # surface pressure [Pa]
    ρ = diagn.physics.surface_air_density[ij]
    V₀ = diagn.physics.surface_wind_speed[ij]

    # TODO actually implement skin temperature?
    T_skin_land = progn.land.soil_temperature[ij, 1]    # uppermost land layer with index 1
    T = diagn.physics.surface_air_temperature[ij]
    land_fraction = model.land_sea_mask[ij]

    # drag coefficient
    d = diagn.physics.boundary_layer_drag[ij]
    drag_land = heat_flux.use_boundary_layer_drag ? d : heat_flux.drag

    # SPEEDY documentation Eq. 54/56, land/sea fraction included
    # Only flux from land if available (not NaN) otherwise zero flux
    # leave out *cₚ here but include below to avoid division
    flux_land  = isfinite(T_skin_land) ? ρ*drag_land*V₀*(T_skin_land  - T) : zero(T_skin_land)
    
    # store without weighting by land fraction for coupling
    diagn.physics.land.sensible_heat_flux[ij] = flux_land*cₚ  # store land flux separately too
    flux_land *= land_fraction                  # weight by land fraction of land-sea mask
    
    # output/diagnose: ocean sets flux (=), land accumulates (+=)
    diagn.physics.sensible_heat_flux[ij] += flux_land*cₚ
    
    # accumulate with += into end=lowermost layer total flux
    nlayers = model.geometry.nlayers
    diagn.tendencies.temp_tend_grid[ij, nlayers] += surface_flux_to_tendency(flux_land, pₛ, model)
    return nothing
end

Adapt.@adapt_structure SurfaceLandHeatFlux

## ----

export PrescribedOceanHeatFlux

"""Prescribed surface sensible heat flux over ocean. Reads heat flux values
from external data sources instead of calculating them dynamically. Useful for
coupled model runs where ocean heat fluxes are provided by an ocean model
or for sensitivity experiments with fixed flux patterns."""
struct PrescribedOceanHeatFlux <: AbstractSurfaceHeatFlux end
PrescribedOceanHeatFlux(::SpectralGrid) = PrescribedOceanHeatFlux()
initialize!(::PrescribedOceanHeatFlux, ::PrimitiveEquation) = nothing

function surface_heat_flux!(ij, diagn, progn, ::PrescribedOceanHeatFlux, model)
    land_fraction = model.land_sea_mask[ij]
    pₛ = diagn.grid.pres_grid_prev[ij]          # surface pressure [Pa]

    # read in a prescribed flux
    flux_ocean = progn.ocean.sensible_heat_flux[ij]
    
    # store without weighting by land fraction for coupling
    diagn.physics.ocean.sensible_heat_flux[ij] = flux_ocean    # store ocean flux separately too
    
    flux_ocean *= (1-land_fraction)             # weight by ocean fraction of land-sea mask
    
    # output/diagnose: ocean sets flux (=), land accumulates (+=)
    diagn.physics.sensible_heat_flux[ij] = flux_ocean
    
    # accumulate with += into end=lowermost layer total flux
    nlayers = model.geometry.nlayers
    diagn.tendencies.temp_tend_grid[ij, nlayers] += surface_flux_to_tendency(flux_ocean, pₛ, model)
    return nothing
end

Adapt.@adapt_structure PrescribedOceanHeatFlux

## ----

export PrescribedLandHeatFlux

"""Prescribed surface sensible heat flux over land. Reads heat flux values
from external data sources instead of calculating them dynamically. Useful for
coupled model runs where land heat fluxes are provided by a land model
or for sensitivity experiments with fixed flux patterns."""
struct PrescribedLandHeatFlux <: AbstractSurfaceHeatFlux end
PrescribedLandHeatFlux(::SpectralGrid) = PrescribedLandHeatFlux()
initialize!(::PrescribedLandHeatFlux, ::PrimitiveEquation) = nothing

function surface_heat_flux!(ij, diagn, progn, ::PrescribedLandHeatFlux, model)
    land_fraction = model.land_sea_mask[ij]
    pₛ = diagn.grid.pres_grid_prev[ij]          # surface pressure [Pa]

    # read in a prescribed flux
    flux_land = progn.land.sensible_heat_flux[ij]
    
    # store without weighting by land fraction for coupling
    diagn.physics.land.sensible_heat_flux[ij] = flux_land  # store land flux separately too

    flux_land *= land_fraction                  # weight by land fraction of land-sea mask
    
    # output/diagnose: ocean sets flux (=), land accumulates (+=)
    diagn.physics.sensible_heat_flux[ij] += flux_land
    
    # accumulate with += into end=lowermost layer total flux
    nlayers = model.geometry.nlayers
    diagn.tendencies.temp_tend_grid[ij, nlayers] += surface_flux_to_tendency(flux_land, pₛ, model)
    return nothing
end 

Adapt.@adapt_structure PrescribedLandHeatFlux
