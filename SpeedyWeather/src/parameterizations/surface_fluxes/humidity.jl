abstract type AbstractSurfaceHumidityFlux <: AbstractParameterization end

export SurfaceHumidityFlux

"""Composite surface humidity flux type that holds flux
for ocean and land (each of any type) separately. Fields are
$(TYPEDFIELDS)"""
@parameterized @kwdef struct SurfaceHumidityFlux{Ocean, Land} <: AbstractSurfaceHumidityFlux
    "[OPTION] Surface humidity flux parameterization for ocean surfaces"
    @component ocean::Ocean

    "[OPTION] Surface humidity flux parameterization for land surfaces"
    @component land::Land
end

Adapt.@adapt_structure SurfaceHumidityFlux

# abbreviate show for SurfaceHumidityFlux to avoid printing all fields of ocean and land albedo in the main show
Base.show(io::IO, SHF::SurfaceHumidityFlux) = show(io, SHF, values = false)

function SurfaceHumidityFlux(
        SG::SpectralGrid;
        ocean = SurfaceOceanHumidityFlux(SG),
        land = SurfaceLandHumidityFlux(SG)
    )
    return SurfaceHumidityFlux(; ocean, land)
end

# determine variables by combining ocean and land variables
variables(SHF::SurfaceHumidityFlux) = (
    ParameterizationVariable(:surface_humidity_flux, Grid2D(), desc = "Total surface humidity flux", units = "kg/m²/s"),
    variables(SHF.ocean)...,
    variables(SHF.land)...,
)

# initialize both
function initialize!(S::SurfaceHumidityFlux, model::PrimitiveWet)
    initialize!(S.ocean, model)
    initialize!(S.land, model)
    return nothing
end

# call first ocean flux then land flux
# order important as ocean sets the flux (=) and land accumulates (+=)
@propagate_inbounds function parameterization!(ij, vars, humidity_flux::SurfaceHumidityFlux, model)
    surface_humidity_flux!(ij, vars, humidity_flux.ocean, model)
    surface_humidity_flux!(ij, vars, humidity_flux.land, model)
    return nothing
end

export SurfaceOceanHumidityFlux
"""Humidity flux parameterization over ocean surfaces. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct SurfaceOceanHumidityFlux{NF} <: AbstractSurfaceHumidityFlux
    "[OPTION] Use drag coefficient from calculated following model.boundary_layer_drag"
    use_boundary_layer_drag::Bool = true

    "[OPTION] Or fixed drag coefficient for humidity flux over ocean"
    @param drag::NF = 0.9e-3 (bounds = 0 .. 1,)

    "[OPTION] Sea ice insulating surface humidity fluxes [1]"
    @param sea_ice_insulation::NF = 0.01 (bounds = 0 .. 1,)
end

Adapt.@adapt_structure SurfaceOceanHumidityFlux
SurfaceOceanHumidityFlux(SG::SpectralGrid; kwargs...) = SurfaceOceanHumidityFlux{SG.NF}(; kwargs...)
initialize!(::SurfaceOceanHumidityFlux, ::PrimitiveWet) = nothing

variables(::SurfaceOceanHumidityFlux) = (
    ParameterizationVariable(:surface_humidity_flux, Grid2D(), desc = "Ocean surface humidity flux", units = "kg/m²/s", namespace = :ocean),
)

@propagate_inbounds function surface_humidity_flux!(
        ij, vars, humidity_flux::SurfaceOceanHumidityFlux, model
    )

    surface = model.geometry.nlayers
    SST = vars.prognostic.ocean.sea_surface_temperature[ij]

    # SATURATION HUMIDITY OVER OCEAN
    pₛ = vars.grid.pres_prev[ij]          # surface pressure [Pa]
    sat_humid_ocean = saturation_humidity(SST, pₛ, model.atmosphere)

    ρ = vars.parameterizations.surface_air_density[ij]
    V₀ = vars.parameterizations.surface_wind_speed[ij]
    land_fraction = model.land_sea_mask.mask[ij]
    surface_humid = vars.grid.humid_grid_prev[ij, surface]
    sea_ice_concentration = haskey(vars.prognostic.ocean, :sea_ice_concentration) ?
        vars.prognostic.ocean.sea_ice_concentration[ij] : zero(SST)

    # drag coefficient either from SurfaceHumidityFlux or from a central drag coefficient
    d = vars.parameterizations.boundary_layer_drag[ij]
    drag_ocean = humidity_flux.use_boundary_layer_drag ? d : humidity_flux.drag

    # SPEEDY documentation eq. 55/57, zero flux if sea surface temperature not available
    # but remove the max( ,0) to allow for surface condensation
    flux_ocean = isfinite(SST) ? ρ * drag_ocean * V₀ * (sat_humid_ocean - surface_humid) : zero(SST)

    # sea ice insulation: more sea ice ⇒ smaller flux (ℵ / ℵ₀ scaling)
    flux_ocean /= 1 + sea_ice_concentration / humidity_flux.sea_ice_insulation

    # store without weighting by land fraction for coupling
    vars.parameterizations.ocean.surface_humidity_flux[ij] = flux_ocean
    flux_ocean *= (1 - land_fraction)             # weight by ocean fraction of land-sea mask

    # output/diagnose: ocean sets flux (=), land accumulates (+=)
    vars.parameterizations.surface_humidity_flux[ij] = flux_ocean

    # accumulate with += into end=lowermost layer total flux
    vars.tendencies.grid.humid[ij, surface] += surface_flux_to_tendency(flux_ocean, pₛ, model)
    return nothing
end

export SurfaceLandHumidityFlux
"""Humidity flux parameterization over land surfaces. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct SurfaceLandHumidityFlux{NF} <: AbstractSurfaceHumidityFlux
    "[OPTION] Use column.boundary_layer_drag coefficient"
    use_boundary_layer_drag::Bool = true

    "[OPTION] Or fixed drag coefficient for humidity flux over land"
    @param drag::NF = 1.2e-3 (bounds = 0 .. 1,)

    "[OPTION] Snow insulation depth [m], e-folding depth controlling how snow insulates surface humidity fluxes"
    @param snow_insulation_depth::NF = 0.05 (bounds = 0 .. 1,)
end

Adapt.@adapt_structure SurfaceLandHumidityFlux
SurfaceLandHumidityFlux(SG::SpectralGrid; kwargs...) = SurfaceLandHumidityFlux{SG.NF}(; kwargs...)
initialize!(::SurfaceLandHumidityFlux, ::PrimitiveWet) = nothing

variables(::SurfaceLandHumidityFlux) = (
    ParameterizationVariable(:surface_humidity_flux, Grid2D(), desc = "Land surface humidity flux", units = "kg/m²/s", namespace = :land),
)

@propagate_inbounds function surface_humidity_flux!(ij, vars, humidity_flux::SurfaceLandHumidityFlux, model)

    # TODO use a skin temperature?
    T = vars.prognostic.land.soil_temperature[ij, 1]  # uppermost land layer with index 1
    snow_depth = vars.prognostic.land.snow_depth[ij]
    α = vars.parameterizations.land.soil_moisture_availability[ij]

    # SATURATION HUMIDITY OVER LAND
    pₛ = vars.grid.pres_prev[ij]          # surface pressure [Pa]
    sat_humid_land = saturation_humidity(T, pₛ, model.atmosphere)

    ρ = vars.parameterizations.surface_air_density[ij]
    V₀ = vars.parameterizations.surface_wind_speed[ij]
    land_fraction = model.land_sea_mask.mask[ij]
    surface = model.geometry.nlayers            # indexing top to bottom
    surface_humid = vars.grid.humid_grid_prev[ij, surface]

    # drag coefficient either from SurfaceLandHumidityFlux or from a central drag coefficient
    d = vars.parameterizations.boundary_layer_drag[ij]
    drag_land = humidity_flux.use_boundary_layer_drag ? d : humidity_flux.drag

    # SPEEDY documentation eq. 55/57, zero flux if land / soil moisture availability not available (=ocean)
    # but remove the max( ,0) to allow for surface condensation
    flux_land = ifelse(isfinite(T) && isfinite(α), ρ * drag_land * V₀ * (α * sat_humid_land - surface_humid), zero(T))

    # snow insulation: deeper snow ⇒ smaller flux (S / S₀ depth scaling)
    flux_land /= 1 + snow_depth / humidity_flux.snow_insulation_depth

    # store without weighting by land fraction for coupling
    vars.parameterizations.land.surface_humidity_flux[ij] = flux_land    # [kg/m²/s]
    flux_land *= land_fraction             # weight by land fraction of land-sea mask

    # output/diagnose: ocean sets flux (=), land accumulates (+=)
    vars.parameterizations.surface_humidity_flux[ij] += flux_land

    # accumulate with += into end=lowermost layer total flux
    vars.tendencies.grid.humid[ij, surface] += surface_flux_to_tendency(flux_land, pₛ, model)
    return nothing
end

export PrescribedOceanHumidityFlux
"""Prescribed humidity flux over ocean surfaces. Applies surface humidity flux from
`progn.ocean.surface_humidity_flux`. $(TYPEDFIELDS)"""
struct PrescribedOceanHumidityFlux <: AbstractSurfaceHumidityFlux end
Adapt.@adapt_structure PrescribedOceanHumidityFlux
PrescribedOceanHumidityFlux(::SpectralGrid) = PrescribedOceanHumidityFlux()
initialize!(::PrescribedOceanHumidityFlux, ::PrimitiveWet) = nothing

variables(::PrescribedOceanHumidityFlux) = (
    PrognosticVariable(:surface_humidity_flux, Grid2D(), desc = "Prescribed ocean surface humidity flux", units = "kg/m²/s", namespace = :ocean),
    ParameterizationVariable(:surface_humidity_flux, Grid2D(), desc = "Ocean surface humidity flux", units = "kg/m²/s", namespace = :ocean),
)


@propagate_inbounds function surface_humidity_flux!(ij, vars, ::PrescribedOceanHumidityFlux, model)
    land_fraction = model.land_sea_mask.mask[ij]
    pₛ = vars.grid.pres_prev[ij]          # surface pressure [Pa]
    surface = model.geometry.nlayers

    # read in a prescribed flux
    flux_ocean = vars.prognostic.ocean.surface_humidity_flux[ij]

    # store without weighting by land fraction for coupling
    vars.parameterizations.ocean.surface_humidity_flux[ij] = flux_ocean
    flux_ocean *= (1 - land_fraction)             # weight by ocean fraction of land-sea mask

    # output/diagnose: ocean sets flux (=), land accumulates (+=)
    vars.parameterizations.surface_humidity_flux[ij] = flux_ocean

    # accumulate with += into end=lowermost layer total flux
    vars.tendencies.grid.humid[ij, surface] += surface_flux_to_tendency(flux_ocean, pₛ, model)
    return nothing
end

export PrescribedLandHumidityFlux
"""Prescribed humidity flux over land surfaces. Applies surface humidity flux from
`progn.land.surface_humidity_flux`. $(TYPEDFIELDS)"""
struct PrescribedLandHumidityFlux <: AbstractSurfaceHumidityFlux end
Adapt.@adapt_structure PrescribedLandHumidityFlux
PrescribedLandHumidityFlux(::SpectralGrid) = PrescribedLandHumidityFlux()
initialize!(::PrescribedLandHumidityFlux, ::PrimitiveWet) = nothing

variables(::PrescribedLandHumidityFlux) = (
    PrognosticVariable(:surface_humidity_flux, Grid2D(), desc = "Prescribed land surface humidity flux", units = "kg/m²/s", namespace = :land),
    ParameterizationVariable(:surface_humidity_flux, Grid2D(), desc = "Land surface humidity flux", units = "kg/m²/s", namespace = :land),
)

@propagate_inbounds function surface_humidity_flux!(ij, vars, ::PrescribedLandHumidityFlux, model)
    land_fraction = model.land_sea_mask.mask[ij]
    pₛ = vars.grid.pres_prev[ij]          # surface pressure [Pa]
    surface = model.geometry.nlayers

    # read in a prescribed flux
    flux_land = vars.prognostic.land.surface_humidity_flux[ij]

    # store without weighting by land fraction for coupling
    vars.parameterizations.land.surface_humidity_flux[ij] = flux_land
    flux_land *= land_fraction                  # weight by ocean fraction of land-sea mask

    # output/diagnose: ocean sets flux (=), land accumulates (+=)
    vars.parameterizations.surface_humidity_flux[ij] += flux_land

    # accumulate with += into end=lowermost layer total flux
    vars.tendencies.grid.humid[ij, surface] += surface_flux_to_tendency(flux_land, pₛ, model)
    return nothing
end