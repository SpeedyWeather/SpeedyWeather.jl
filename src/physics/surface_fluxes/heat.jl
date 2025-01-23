export SurfaceHeatFlux
@kwdef struct SurfaceHeatFlux{Ocean, Land} <: AbstractSurfaceHeatFlux
    ocean::Ocean = SurfaceOceanHeatFlux()
    land::Land = SurfaceLandHeatFlux()
end

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

function surface_heat_flux!(   
    column::ColumnVariables,
    heat_flux::SurfaceHeatFlux,
    diagn::DiagnosticVariables,
    model::PrimitiveEquation,
)   
    surface_heat_flux!(column, heat_flux.ocean, diagn, model)
    surface_heat_flux!(column, heat_flux.land, diagn, model)
end

## ----

export NoSurfaceHeatFlux
struct NoSurfaceHeatFlux <: AbstractSurfaceHeatFlux end
NoSurfaceHeatFlux(::SpectralGrid) = NoSurfaceHeatFlux()
initialize!(::NoSurfaceHeatFlux, ::PrimitiveEquation) = nothing
surface_heat_flux!(::ColumnVariables, ::NoSurfaceHeatFlux, ::PrimitiveEquation) = nothing

## ----

export SurfaceOceanHeatFlux
@kwdef struct SurfaceOceanHeatFlux{NF} <: AbstractSurfaceHeatFlux
    "Use (possibly) flow-dependent column.boundary_layer_drag coefficient"
    use_boundary_layer_drag::Bool = true

    "Otherwise, use the following drag coefficient for heat fluxes over ocean"
    heat_exchange::NF = 0.9e-3
end

SurfaceOceanHeatFlux(SG::SpectralGrid; kwargs...) = SurfaceOceanHeatFlux{SG.NF}(; kwargs...)
initialize!(::SurfaceOceanHeatFlux, ::PrimitiveEquation) = nothing

function surface_heat_flux!(
    column::ColumnVariables,
    heat_flux::SurfaceOceanHeatFlux,
    diagn::DiagnosticVariables,
    model::PrimitiveEquation,
)
    cₚ = model.atmosphere.heat_capacity
    (; heat_exchange) = heat_flux

    ρ = column.surface_air_density
    V₀ = column.surface_wind_speed
    T_skin_ocean = column.skin_temperature_sea
    T = column.surface_temp
    land_fraction = column.land_fraction

    # drag coefficient
    drag_ocean = heat_flux.use_boundary_layer_drag ? column.boundary_layer_drag : heat_exchange

    # SPEEDY documentation Eq. 54/56, land/sea fraction included
    # Only flux from sea if available (not NaN) otherwise zero flux
    flux_ocean  = isfinite(T_skin_ocean) ? ρ*drag_ocean*V₀*cₚ*(T_skin_ocean  - T)*(1-land_fraction) : 0
    column.flux_temp_upward[end] += flux_ocean
    column.sensible_heat_flux = flux_ocean      # ocean sets the flux (=), land accumulates (+=)

    return nothing
end

## ----

export SurfaceLandHeatFlux
@kwdef struct SurfaceLandHeatFlux{NF} <: AbstractSurfaceHeatFlux
    "Use (possibly) flow-dependent column.boundary_layer_drag coefficient"
    use_boundary_layer_drag::Bool = true

    "Otherwise, use the following drag coefficient for heat fluxes over land"
    heat_exchange::NF = 1.2e-3    # for neutral stability
end

SurfaceLandHeatFlux(SG::SpectralGrid; kwargs...) = SurfaceLandHeatFlux{SG.NF}(; kwargs...)
initialize!(::SurfaceLandHeatFlux, ::PrimitiveEquation) = nothing

function surface_heat_flux!(
    column::ColumnVariables,
    heat_flux::SurfaceLandHeatFlux,
    diagn::DiagnosticVariables,
    model::PrimitiveEquation,
)
    cₚ = model.atmosphere.heat_capacity
    (; heat_exchange) = heat_flux

    ρ = column.surface_air_density
    V₀ = column.surface_wind_speed
    T_skin_land = column.skin_temperature_land
    T = column.surface_temp
    land_fraction = column.land_fraction

    # drag coefficient
    drag_land = heat_flux.use_boundary_layer_drag ? column.boundary_layer_drag : heat_exchange

    # SPEEDY documentation Eq. 54/56, land/sea fraction included
    # Only flux from sea if available (not NaN) otherwise zero flux
    flux_land  = isfinite(T_skin_land) ? ρ*drag_land*V₀*cₚ*(T_skin_land  - T)*land_fraction : 0
    column.flux_temp_upward[end] += flux_land
    column.sensible_heat_flux += flux_land      # ocean sets the flux (=), land accumulates (+=)

    return nothing
end

## ----

export PrescribedOceanHeatFlux
struct PrescribedOceanHeatFlux <: AbstractSurfaceHeatFlux end
PrescribedOceanHeatFlux(::SpectralGrid) = PrescribedOceanHeatFlux()
initialize!(::PrescribedOceanHeatFlux, ::PrimitiveEquation) = nothing

function surface_heat_flux!(
    column::ColumnVariables,
    fluxes::PrescribedOceanHeatFlux,
    diagn::DiagnosticVariables,
    model::PrimitiveEquation,
)
    land_fraction = column.land_fraction

    # read in a prescribed flux
    flux = diagn.physics.sensible_heat_flux[column.ij]*(1-land_fraction)
    column.flux_temp_upward[end] += flux    # end=lowermost layer
    column.sensible_heat_flux = flux        # ocean sets the flux (=), land accumulates (+=)
end

## ----

export PrescribedLandHeatFlux
struct PrescribedLandHeatFlux <: AbstractSurfaceHeatFlux end
PrescribedLandHeatFlux(::SpectralGrid) = PrescribedLandHeatFlux()
initialize!(::PrescribedLandHeatFlux, ::PrimitiveEquation) = nothing

function surface_heat_flux!(
    column::ColumnVariables,
    fluxes::PrescribedLandHeatFlux,
    diagn::DiagnosticVariables,
    model::PrimitiveEquation,
)
    land_fraction = column.land_fraction

    # read in a prescribed flux
    flux = diagn.physics.sensible_heat_flux[column.ij]*land_fraction
    column.flux_temp_upward[end] += flux    # end=lowermost layer
    column.sensible_heat_flux += flux        # ocean sets the flux (=), land accumulates (+=)
end