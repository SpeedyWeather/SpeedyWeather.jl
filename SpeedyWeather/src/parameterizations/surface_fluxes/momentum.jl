abstract type AbstractSurfaceMomentumFlux <: AbstractParameterization end

export SurfaceMomentumFlux

"""Surface momentum flux parameterization for wind stress. Calculates the 
turbulent exchange of momentum between the surface and atmosphere, representing
surface drag that slows down near-surface winds. Uses bulk aerodynamic formulas
with separate drag coefficients for ocean and land surfaces to account for
different surface roughness characteristics. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct SurfaceMomentumFlux{NF, FI <: AbstractFluxIntegration} <: AbstractSurfaceMomentumFlux
    "[OPTION] Use drag coefficient from calculated following model.boundary_layer.drag"
    use_boundary_layer_drag::Bool = true

    "[OPTION] Or fixed drag coefficient for momentum fluxes over ocean"
    @param drag_ocean::NF = 1.8e-3 (bounds = 0 .. 1,)

    "[OPTION] Or fixed drag coefficient for momentum fluxes over land"
    @param drag_land::NF = 2.4e-3 (bounds = 0 .. 1,)

    "[OPTION] Time integration of the flux into the lowermost layer"
    time_integration::FI = DEFAULT_FLUX_INTEGRATION()
end

Adapt.@adapt_structure SurfaceMomentumFlux
SurfaceMomentumFlux(SG::SpectralGrid; time_integration = DEFAULT_FLUX_INTEGRATION(), kwargs...) =
    SurfaceMomentumFlux{SG.NF, typeof(time_integration)}(; time_integration, kwargs...)
initialize!(::SurfaceMomentumFlux, ::PrimitiveEquation) = nothing

# function barrier
@propagate_inbounds function parameterization!(ij, vars, momentum_flux::SurfaceMomentumFlux, model)
    surface_wind_stress!(ij, vars, momentum_flux, model)
    return nothing
end

@propagate_inbounds function surface_wind_stress!(ij, vars, momentum_flux::SurfaceMomentumFlux, model)

    (; drag_land, drag_ocean) = momentum_flux
    land_fraction = model.land_sea_mask.land_fraction[ij]
    surface = model.geometry.nlayers

    d_land = vars.parameterizations.land.boundary_layer_drag_momentum[ij]
    d_ocean = vars.parameterizations.ocean.boundary_layer_drag_momentum[ij]
    d_combined = land_fraction * d_land + (1 - land_fraction) * d_ocean
    
    # drag coefficient either from SurfaceMomentumFlux or from a central drag coefficient
    drag = ifelse(momentum_flux.use_boundary_layer_drag, d_combined, land_fraction * drag_land + (1 - land_fraction) * drag_ocean)

    # Fortran SPEEDY documentation eq. 52, 53, accumulate fluxes with +=
    V₀ = vars.parameterizations.surface_wind_speed[ij]
    ρ = vars.parameterizations.surface_air_density[ij]
    uₛ = vars.parameterizations.surface_wind_u[ij]
    vₛ = vars.parameterizations.surface_wind_v[ij]

    # flux into lowermost layer [kg/m³ * m/s * m/s = kg/(m·s²) = Pa]
    conductance = ρ * drag * V₀ * f             # = -∂flux/∂u [kg/m²/s]
    flux_u_upward = -conductance * uₛ
    flux_v_upward = -conductance * vₛ

    # convert fluxes to tendencies
    pₛ = vars.parameterizations.surface_pressure[ij]       # surface pressure [Pa]
    u_tend = get_tendency_step(vars.tendencies.grid.u, model.time_stepping, momentum_flux)
    v_tend = get_tendency_step(vars.tendencies.grid.v, model.time_stepping, momentum_flux)

    scheme = momentum_flux.time_integration
    u_tend[ij, surface] += surface_flux_to_tendency(scheme, flux_u_upward, conductance, pₛ, model)
    v_tend[ij, surface] += surface_flux_to_tendency(scheme, flux_v_upward, conductance, pₛ, model)
    return nothing
end
