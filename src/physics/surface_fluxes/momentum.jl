abstract type AbstractSurfaceMomentumFlux <: AbstractParameterization end

export SurfaceMomentumFlux
@kwdef struct SurfaceMomentumFlux{NF} <: AbstractSurfaceMomentumFlux
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

    # drag coefficient either from SurfaceMomentumFlux or from a central drag coefficient
    d = diagn.physics.boundary_layer_drag[ij]
    drag = surface_wind.use_boundary_layer_drag ? d : land_fraction*drag_land + (1-land_fraction)*drag_ocean
    
    # Fortran SPEEDY documentation eq. 52, 53, accumulate fluxes with +=
    V₀ = diagn.physics.surface_wind_speed
    flux_u_upward = -ρ*drag*V₀*surface_u
    flux_v_upward = -ρ*drag*V₀*surface_v
    
    # convert fluxes to tendencies
    pₛ = diagn.grid.pres_grid_prev[ij]          # surface pressure [Pa]
    diagn.tendencies.u_tend_grid[ij, end] += flux_to_tendency(flux_u_upward, pₛ, model)
    diagn.tendencies.v_tend_grid[ij, end] += flux_to_tendency(flux_v_upward, pₛ, model)

    return nothing
end