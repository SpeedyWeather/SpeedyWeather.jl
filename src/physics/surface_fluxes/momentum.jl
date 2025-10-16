abstract type AbstractSurfaceMomentumFlux <: AbstractParameterization end

export SurfaceMomentumFlux

"""Surface momentum flux parameterization for wind stress. Calculates the 
turbulent exchange of momentum between the surface and atmosphere, representing
surface drag that slows down near-surface winds. Uses bulk aerodynamic formulas
with separate drag coefficients for ocean and land surfaces to account for
different surface roughness characteristics. Fields are $(TYPEDFIELDS)"""
@kwdef struct SurfaceMomentumFlux{NF} <: AbstractSurfaceMomentumFlux
    "[OPTION] Near-surface wind slowdown"
    wind_slowdown::NF = 0.95

    "[OPTION] Use drag coefficient from calculated following model.boundary_layer_drag"
    use_boundary_layer_drag::Bool = true

    "[OPTION] Or fixed drag coefficient for momentum fluxes over ocean"
    drag_ocean::NF = 1.8e-3

    "[OPTION] Or fixed drag coefficient for momentum fluxes over land"
    drag_land::NF = 2.4e-3
end

SurfaceMomentumFlux(SG::SpectralGrid; kwargs...) = SurfaceMomentumFlux{SG.NF}(; kwargs...)
initialize!(::SurfaceMomentumFlux, ::PrimitiveEquation) = nothing

# function barrier
function parameterization!(ij, diagn, progn, momentum_flux::SurfaceMomentumFlux, model)
    surface_wind_stress!(ij, diagn, momentum_flux, model)
end

function surface_wind_stress!(ij, diagn, momentum_flux::SurfaceMomentumFlux, model)

    (; drag_land, drag_ocean) = momentum_flux
    f = momentum_flux.wind_slowdown

    # drag coefficient either from SurfaceMomentumFlux or from a central drag coefficient
    d = diagn.physics.boundary_layer_drag[ij]
    drag = momentum_flux.use_boundary_layer_drag ? d : land_fraction*drag_land + (1-land_fraction)*drag_ocean
    
    # Fortran SPEEDY documentation eq. 52, 53, accumulate fluxes with +=
    V₀ = diagn.physics.surface_wind_speed[ij]
    ρ = diagn.physics.surface_air_density[ij]
    surface_u = diagn.grid.u_grid_prev[ij, end]
    surface_v = diagn.grid.v_grid_prev[ij, end]

    flux_u_upward = -ρ*drag*V₀*f*surface_u
    flux_v_upward = -ρ*drag*V₀*f*surface_v
    
    # convert fluxes to tendencies
    pₛ = diagn.grid.pres_grid_prev[ij]          # surface pressure [Pa]
    diagn.tendencies.u_tend_grid[ij, end] += surface_flux_to_tendency(flux_u_upward, pₛ, model)
    diagn.tendencies.v_tend_grid[ij, end] += surface_flux_to_tendency(flux_v_upward, pₛ, model)

    return nothing
end