abstract type AbstractSurfaceMomentumFlux <: AbstractParameterization end

export SurfaceMomentumFlux

"""Surface momentum flux parameterization for wind stress. Calculates the 
turbulent exchange of momentum between the surface and atmosphere, representing
surface drag that slows down near-surface winds. Uses bulk aerodynamic formulas
with separate drag coefficients for ocean and land surfaces to account for
different surface roughness characteristics. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct SurfaceMomentumFlux{NF} <: AbstractSurfaceMomentumFlux
    "[OPTION] Near-surface wind slowdown"
    @param wind_slowdown::NF = 0.95 (bounds = 0 .. 1,)

    "[OPTION] Use drag coefficient from calculated following model.boundary_layer_drag"
    use_boundary_layer_drag::Bool = true

    "[OPTION] Or fixed drag coefficient for momentum fluxes over ocean"
    @param drag_ocean::NF = 1.8e-3 (bounds = 0 .. 1,)

    "[OPTION] Or fixed drag coefficient for momentum fluxes over land"
    @param drag_land::NF = 2.4e-3 (bounds = 0 .. 1,)
end

Adapt.@adapt_structure SurfaceMomentumFlux
SurfaceMomentumFlux(SG::SpectralGrid; kwargs...) = SurfaceMomentumFlux{SG.NF}(; kwargs...)
initialize!(::SurfaceMomentumFlux, ::PrimitiveEquation) = nothing

function variables(::SurfaceMomentumFlux)
    return (
        ParameterizationVariable(:boundary_layer_drag, Grid2D(), desc = "Boundary layer drag coefficient", units = "1"),
        ParameterizationVariable(:surface_wind_speed, Grid2D(), desc = "Surface wind speed", units = "m/s"),
        ParameterizationVariable(:surface_air_density, Grid2D(), desc = "Surface air density", units = "kg/m^2"),
        ParameterizationVariable(:surface_air_temperature, Grid2D(), desc = "Surface air temperature", units = "K"),
    )
end

# function barrier
@propagate_inbounds function parameterization!(ij, vars, momentum_flux::SurfaceMomentumFlux, model)
    surface_wind_stress!(ij, vars, momentum_flux, model)
    return nothing
end

@propagate_inbounds function surface_wind_stress!(ij, vars, momentum_flux::SurfaceMomentumFlux, model)

    (; drag_land, drag_ocean) = momentum_flux
    land_fraction = model.land_sea_mask.mask[ij]
    surface = model.geometry.nlayers

    # drag coefficient either from SurfaceMomentumFlux or from a central drag coefficient
    d = vars.parameterizations.boundary_layer_drag[ij]
    drag = ifelse(momentum_flux.use_boundary_layer_drag, d, land_fraction * drag_land + (1 - land_fraction) * drag_ocean)

    # Fortran SPEEDY documentation eq. 52, 53, accumulate fluxes with +=
    V₀ = vars.parameterizations.surface_wind_speed[ij]
    ρ = vars.parameterizations.surface_air_density[ij]
    u = vars.grid.u_prev[ij, surface]
    v = vars.grid.v_prev[ij, surface]

    # fraction to slow down the lowermost layer wind u,v to surface wind
    f = momentum_flux.wind_slowdown

    # flux into lowermost layer [kg/m³ * m/s * m/s = kg/(m·s²) = Pa]
    flux_u_upward = -ρ * drag * V₀ * f * u
    flux_v_upward = -ρ * drag * V₀ * f * v

    # convert fluxes to tendencies
    pₛ = vars.grid.pres_prev[ij]          # surface pressure [Pa]
    vars.tendencies.grid.u[ij, surface] += surface_flux_to_tendency(flux_u_upward, pₛ, model)
    vars.tendencies.grid.v[ij, surface] += surface_flux_to_tendency(flux_v_upward, pₛ, model)
    return nothing
end
