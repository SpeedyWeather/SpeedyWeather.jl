abstract type AbstractSurfaceMomentumFlux <: AbstractParameterization end

export SurfaceMomentumFlux

"""Surface momentum flux parameterization for wind stress. Calculates the 
turbulent exchange of momentum between the surface and atmosphere, representing
surface drag that slows down near-surface winds. Uses bulk aerodynamic formulas
with separate drag coefficients for ocean and land surfaces to account for
different surface roughness characteristics. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct SurfaceMomentumFlux{NF} <: AbstractSurfaceMomentumFlux
    "[OPTION] Near-surface wind slowdown"
    @param wind_slowdown::NF = 0.95 (bounds=0..1,)

    "[OPTION] Use drag coefficient from calculated following model.boundary_layer_drag"
    use_boundary_layer_drag::Bool = true

    "[OPTION] Or fixed drag coefficient for momentum fluxes over ocean"
    @param drag_ocean::NF = 1.8e-3 (bounds=0..1,)

    "[OPTION] Or fixed drag coefficient for momentum fluxes over land"
    @param drag_land::NF = 2.4e-3 (bounds=0..1,)
end

Adapt.@adapt_structure SurfaceMomentumFlux
SurfaceMomentumFlux(SG::SpectralGrid; kwargs...) = SurfaceMomentumFlux{SG.NF}(; kwargs...)
initialize!(::SurfaceMomentumFlux, ::PrimitiveEquation) = nothing

function variables(::SurfaceMomentumFlux)
    return (
        DiagnosticVariable(name = :boundary_layer_drag, dims = Grid2D(), desc = "Boundary layer drag coefficient", units = "1"),
        DiagnosticVariable(name = :surface_wind_speed, dims = Grid2D(), desc = "Surface wind speed", units = "m/s"),
        DiagnosticVariable(name = :surface_air_density, dims = Grid2D(), desc = "Surface air density", units = "kg/m³"),
        DiagnosticVariable(name = :surface_air_temperature, dims = Grid2D(), desc = "Surface air temperature", units = "K"),
    )
end

# function barrier
@propagate_inbounds function parameterization!(ij, diagn, progn, momentum_flux::SurfaceMomentumFlux, model)
    return surface_wind_stress!(ij, diagn, momentum_flux, model)
end

@propagate_inbounds function surface_wind_stress!(ij, diagn, momentum_flux::SurfaceMomentumFlux, model)

    (; drag_land, drag_ocean) = momentum_flux
    land_fraction = model.land_sea_mask.mask[ij]
    surface = model.geometry.nlayers

    # drag coefficient either from SurfaceMomentumFlux or from a central drag coefficient
    d = diagn.physics.boundary_layer_drag[ij]
    drag = momentum_flux.use_boundary_layer_drag ?
        d : land_fraction * drag_land + (1 - land_fraction) * drag_ocean

    # Fortran SPEEDY documentation eq. 52, 53, accumulate fluxes with +=
    V₀ = diagn.physics.surface_wind_speed[ij]
    ρ = diagn.physics.surface_air_density[ij]
    u = diagn.grid.u_grid_prev[ij, surface]
    v = diagn.grid.v_grid_prev[ij, surface]

    # fraction to slow down the lowermost layer wind u,v to surface wind
    f = momentum_flux.wind_slowdown

    # flux into lowermost layer [kg/m³ * m/s * m/s = kg/(m·s²) = Pa]
    flux_u_upward = -ρ * drag * V₀ * f * u
    flux_v_upward = -ρ * drag * V₀ * f * v

    # convert fluxes to tendencies
    pₛ = diagn.grid.pres_grid_prev[ij]          # surface pressure [Pa]
    diagn.tendencies.u_tend_grid[ij, surface] += surface_flux_to_tendency(flux_u_upward, pₛ, model)
    diagn.tendencies.v_tend_grid[ij, surface] += surface_flux_to_tendency(flux_v_upward, pₛ, model)
    return nothing
end
