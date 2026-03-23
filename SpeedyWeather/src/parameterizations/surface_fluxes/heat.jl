abstract type AbstractSurfaceHeatFlux <: AbstractParameterization end

export SurfaceHeatFlux

"""Composite surface (sensible) heat flux type that holds flux
for ocean and land (each of any type) separately. Fields are
$(TYPEDFIELDS)"""
@parameterized @kwdef struct SurfaceHeatFlux{Ocean, Land} <: AbstractSurfaceHeatFlux
    @component ocean::Ocean = SurfaceOceanHeatFlux()
    @component land::Land = SurfaceLandHeatFlux()
end

Adapt.@adapt_structure SurfaceHeatFlux

# abbreviate show for SurfaceHeatFlux to avoid printing all fields of ocean and land albedo in the main show
Base.show(io::IO, SHF::SurfaceHeatFlux) = show(io, SHF, values = false)

# generator function
function SurfaceHeatFlux(
        SG::SpectralGrid;
        ocean = SurfaceOceanHeatFlux(SG),
        land = SurfaceLandHeatFlux(SG)
    )
    return SurfaceHeatFlux(; ocean, land)
end

# determine variables by combining ocean and land variables
variables(SHF::SurfaceHeatFlux) = (
    ParameterizationVariable(:sensible_heat_flux, Grid2D(), desc = "Total surface sensible heat flux", units = "W/m^2"),
    variables(SHF.ocean)...,
    variables(SHF.land)...,
)

function initialize!(S::SurfaceHeatFlux, model::PrimitiveEquation)
    initialize!(S.ocean, model)
    initialize!(S.land, model)
    return nothing
end

# call first ocean flux then land flux
@propagate_inbounds function parameterization!(ij, vars, heat_flux::SurfaceHeatFlux, model)
    surface_heat_flux!(ij, vars, heat_flux.ocean, model)
    surface_heat_flux!(ij, vars, heat_flux.land, model)
    return nothing
end

export SurfaceOceanHeatFlux

"""Surface sensible heat flux parameterization over ocean. Calculates the 
turbulent exchange of sensible heat between the ocean surface and the atmosphere
based on temperature differences and wind speed. Uses bulk aerodynamic formulas
with drag coefficients. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct SurfaceOceanHeatFlux{NF} <: AbstractSurfaceHeatFlux
    "[OPTION] Use drag coefficient from calculated following model.boundary_layer_drag"
    use_boundary_layer_drag::Bool = true

    "[OPTION] Or fixed drag coefficient for heat fluxes over ocean"
    @param drag::NF = 0.9e-3 (bounds = Nonnegative,)

    "[OPTION] Sea ice insulating surface heat fluxes [1]"
    @param sea_ice_insulation::NF = 0.01 (bounds = Positive,)
end

Adapt.@adapt_structure SurfaceOceanHeatFlux
SurfaceOceanHeatFlux(SG::SpectralGrid; kwargs...) = SurfaceOceanHeatFlux{SG.NF}(; kwargs...)
initialize!(::SurfaceOceanHeatFlux, ::PrimitiveEquation) = nothing

variables(::SurfaceOceanHeatFlux) = (
    ParameterizationVariable(:sensible_heat_flux, Grid2D(), desc = "Ocean sensible heat flux", units = "W/m^2", namespace = :ocean),
)

@propagate_inbounds function surface_heat_flux!(ij, vars, heat_flux::SurfaceOceanHeatFlux, model)

    nlayers = model.geometry.nlayers
    cₚ = model.atmosphere.heat_capacity
    ρ = vars.parameterizations.surface_air_density[ij]
    V₀ = vars.parameterizations.surface_wind_speed[ij]
    sea_ice_concentration = haskey(vars.prognostic.ocean, :sea_ice_concentration) ?
        vars.prognostic.ocean.sea_ice_concentration[ij] : zero(ρ)

    # TODO actually implement skin temperature?
    SST = vars.prognostic.ocean.sea_surface_temperature[ij]
    T = vars.parameterizations.surface_air_temperature[ij]
    land_fraction = model.land_sea_mask.mask[ij]
    pₛ = vars.grid.pres_prev[ij]

    # drag coefficient
    d = vars.parameterizations.boundary_layer_drag[ij]
    drag_ocean = ifelse(heat_flux.use_boundary_layer_drag, d, heat_flux.drag)

    # SPEEDY documentation Eq. 54/56, land/sea fraction included
    # Only flux from ocean if available (not NaN) otherwise zero flux
    # leave out *cₚ here but include below to avoid division
    flux_ocean = ifelse(isfinite(SST), ρ * drag_ocean * V₀ * (SST - T), zero(SST))

    # sea ice insulation: more sea ice ⇒ smaller flux (ℵ / ℵ₀ scaling)
    flux_ocean /= 1 + sea_ice_concentration / heat_flux.sea_ice_insulation

    # store without weighting by land fraction for coupling [W/m²]
    vars.parameterizations.ocean.sensible_heat_flux[ij] = flux_ocean * cₚ  # to store ocean flux separately too

    flux_ocean *= (1 - land_fraction)                             # weight by ocean fraction of land-sea mask

    # output/diagnostics: ocean sets the flux (=), land accumulates (+=)
    vars.parameterizations.sensible_heat_flux[ij] = flux_ocean * cₚ        # [W/m²]

    # accumulate with += into end=lowermost layer total flux
    vars.tendencies.grid.temp[ij, nlayers] += surface_flux_to_tendency(flux_ocean, pₛ, model)
    return nothing
end

export SurfaceLandHeatFlux

"""Surface sensible heat flux parameterization over land. Calculates the 
turbulent exchange of sensible heat between the land surface and the atmosphere
based on soil temperature differences and wind speed. Uses bulk aerodynamic 
formulas with drag coefficients that can account for surface roughness. 
Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct SurfaceLandHeatFlux{NF} <: AbstractSurfaceHeatFlux
    "[OPTION] Use drag coefficient from calculated following model.boundary_layer_drag"
    use_boundary_layer_drag::Bool = true

    "[OPTION] Or fixed drag coefficient for heat fluxes over land"
    @param drag::NF = 1.2e-3 (bounds = Nonnegative,)       # for neutral stability

    "[OPTION] Snow insulation depth [m], e-folding depth reducing surface heat fluxes"
    @param snow_insulation_depth::NF = 0.05 (bounds = Positive,)
end

Adapt.@adapt_structure SurfaceLandHeatFlux
SurfaceLandHeatFlux(SG::SpectralGrid; kwargs...) = SurfaceLandHeatFlux{SG.NF}(; kwargs...)
initialize!(::SurfaceLandHeatFlux, ::PrimitiveEquation) = nothing

variables(::SurfaceLandHeatFlux) = (
    ParameterizationVariable(:sensible_heat_flux, Grid2D(), desc = "Land sensible heat flux", units = "W/m^2", namespace = :land),
)

@propagate_inbounds function surface_heat_flux!(ij, vars, heat_flux::SurfaceLandHeatFlux, model)

    surface = model.geometry.nlayers
    cₚ = model.atmosphere.heat_capacity
    pₛ = vars.grid.pres_prev[ij]                  # surface pressure [Pa]
    ρ = vars.parameterizations.surface_air_density[ij]
    V₀ = vars.parameterizations.surface_wind_speed[ij]

    # TODO actually implement skin temperature?
    T_skin_land = vars.prognostic.land.soil_temperature[ij, 1]    # uppermost land layer with index 1
    T = vars.parameterizations.surface_air_temperature[ij]
    land_fraction = model.land_sea_mask.mask[ij]
    snow_depth = haskey(vars.prognostic.land, :snow_depth) ? vars.prognostic.land.snow_depth[ij] : zero(T)

    # drag coefficient
    d = vars.parameterizations.boundary_layer_drag[ij]
    drag_land = ifelse(heat_flux.use_boundary_layer_drag, d, heat_flux.drag)

    # SPEEDY documentation Eq. 54/56, land/sea fraction included
    # Only flux from land if available (not NaN) otherwise zero flux
    # leave out *cₚ here but include below to avoid division
    flux_land = ifelse(isfinite(T_skin_land), ρ * drag_land * V₀ * (T_skin_land - T), zero(T))

    # snow insulation: deeper snow ⇒ smaller flux (S / S₀ depth scaling)
    flux_land /= 1 + snow_depth / heat_flux.snow_insulation_depth

    # store without weighting by land fraction for coupling [W/m²]
    vars.parameterizations.land.sensible_heat_flux[ij] = flux_land * cₚ  # store land flux separately too
    flux_land *= land_fraction                  # weight by land fraction of land-sea mask

    # output/diagnose: ocean sets flux (=), land accumulates (+=)
    vars.parameterizations.sensible_heat_flux[ij] += flux_land * cₚ    # [W/m²]

    # accumulate with += into end=lowermost layer total flux
    vars.tendencies.grid.temp[ij, surface] += surface_flux_to_tendency(flux_land, pₛ, model)
    return nothing
end

export PrescribedOceanHeatFlux

"""Prescribed surface sensible heat flux over ocean. Reads heat flux values
from external data sources instead of calculating them dynamically. Useful for
coupled model runs where ocean heat fluxes are provided by an ocean model
or for sensitivity experiments with fixed flux patterns."""
struct PrescribedOceanHeatFlux <: AbstractSurfaceHeatFlux end
Adapt.@adapt_structure PrescribedOceanHeatFlux
PrescribedOceanHeatFlux(::SpectralGrid) = PrescribedOceanHeatFlux()
initialize!(::PrescribedOceanHeatFlux, ::PrimitiveEquation) = nothing

variables(::PrescribedOceanHeatFlux) = (
    PrognosticVariable(:sensible_heat_flux, Grid2D(), desc = "Prescribed ocean sensible heat flux", units = "W/m^2", namespace = :ocean),
    ParameterizationVariable(:sensible_heat_flux, Grid2D(), desc = "Ocean sensible heat flux", units = "W/m^2", namespace = :ocean),
)

@propagate_inbounds function surface_heat_flux!(ij, vars, ::PrescribedOceanHeatFlux, model)
    land_fraction = model.land_sea_mask.mask[ij]
    pₛ = vars.grid.pres_prev[ij]          # surface pressure [Pa]
    cₚ = model.atmosphere.heat_capacity
    surface = model.geometry.nlayers

    # read in a prescribed flux
    flux_ocean = vars.prognostic.ocean.sensible_heat_flux[ij]

    # store without weighting by land fraction for coupling
    vars.parameterizations.ocean.sensible_heat_flux[ij] = flux_ocean    # store ocean flux separately too

    flux_ocean *= (1 - land_fraction)             # weight by ocean fraction of land-sea mask

    # output/diagnose: ocean sets flux (=), land accumulates (+=)
    vars.parameterizations.sensible_heat_flux[ij] = flux_ocean

    # accumulate with += into end=lowermost layer total flux
    flux_ocean /= cₚ                            # [W/m²] -> [K/s]
    vars.tendencies.grid.temp[ij, surface] += surface_flux_to_tendency(flux_ocean, pₛ, model)
    return nothing
end

export PrescribedLandHeatFlux

"""Prescribed surface sensible heat flux over land. Reads heat flux values
from external data sources instead of calculating them dynamically. Useful for
coupled model runs where land heat fluxes are provided by a land model
or for sensitivity experiments with fixed flux patterns."""
struct PrescribedLandHeatFlux <: AbstractSurfaceHeatFlux end
Adapt.@adapt_structure PrescribedLandHeatFlux
PrescribedLandHeatFlux(::SpectralGrid) = PrescribedLandHeatFlux()
initialize!(::PrescribedLandHeatFlux, ::PrimitiveEquation) = nothing

variables(::PrescribedLandHeatFlux) = (
    PrognosticVariable(:sensible_heat_flux, Grid2D(), desc = "Prescribed land sensible heat flux", units = "W/m^2", namespace = :land),
    ParameterizationVariable(:sensible_heat_flux, Grid2D(), desc = "Land sensible heat flux", units = "W/m^2", namespace = :land),
)

@propagate_inbounds function surface_heat_flux!(ij, vars, ::PrescribedLandHeatFlux, model)
    land_fraction = model.land_sea_mask.mask[ij]
    pₛ = vars.grid.pres_prev[ij]                # surface pressure [Pa]
    cₚ = model.atmosphere.heat_capacity
    surface = vars.nlayers                     # indexing top to bottom

    # read in a prescribed flux
    flux_land = vars.prognostic.land.sensible_heat_flux[ij]

    # store without weighting by land fraction for coupling
    vars.physics.land.sensible_heat_flux[ij] = flux_land  # store land flux separately too

    flux_land *= land_fraction                  # weight by land fraction of land-sea mask

    # output/diagnose: ocean sets flux (=), land accumulates (+=)
    vars.physics.sensible_heat_flux[ij] += flux_land

    # accumulate with += into end=lowermost layer total flux
    flux_land /= cₚ                             # [W/m²] -> [K/s]
    vars.tendencies.temp_tend_grid[ij, surface] += surface_flux_to_tendency(flux_land, pₛ, model)
    return nothing
end