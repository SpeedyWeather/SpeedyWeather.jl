abstract type AbstractSurfaceWind <: AbstractParameterization end

export SurfaceWind
@kwdef struct SurfaceWind{NF} <: AbstractSurfaceWind
    "Ratio of near-surface wind to lowest-level wind [1]"
    f_wind::NF = 0.95

    "Wind speed of sub-grid scale gusts [m/s]"
    V_gust::NF = 5

    "Use (possibly) flow-dependent column.boundary_layer_drag coefficient"
    use_boundary_layer_drag::Bool = true

    "Otherwise, drag coefficient over land (orography = 0) [1]"
    drag_land::NF = 2.4e-3
    
    "Otherwise, Drag coefficient over sea [1]"
    drag_sea::NF = 1.8e-3
end

SurfaceWind(SG::SpectralGrid; kwargs...) = SurfaceWind{SG.NF}(; kwargs...)
initialize!(::SurfaceWind, ::PrimitiveEquation) = nothing

# function barrier
function parameterization!(ij, diagn, progn, sw::SurfaceWind, model)
    surface_wind_stress!(ij, diagn, sw, model.land_sea_mask, model.atmosphere, model.planet, model.geometry)
end

function surface_wind_stress!(ij, diagn, surface_wind::SurfaceWind, land_sea_mask, atmosphere, planet, geometry)

    land_fraction = land_sea_mask[ij]
    (; f_wind, V_gust, drag_land, drag_sea) = surface_wind

    # Fortran SPEEDY documentation eq. 49 but use previous time step for numerical stability
    surface_u = f_wind*diagn.grid.u_grid_prev[ij, end] 
    surface_v = f_wind*diagn.grid.v_grid_prev[ij, end]

    # Fortran SPEEDY documentation eq. 50
    surface_wind_speed = sqrt(surface_u^2 + surface_v^2 + V_gust^2)
    diagn.physics.surface_wind_speed[ij] = surface_wind_speed

    # drag coefficient either from SurfaceWind or from a central drag coefficient
    (; boundary_layer_drag) = diagn.physics
    drag = surface_wind.use_boundary_layer_drag ? boundary_layer_drag :
                                land_fraction*drag_land + (1-land_fraction)*drag_sea

    (; R_dry, κ) = atmosphere
    pₛ = diagn.grid.pres_grid_prev[ij]          # surface pressure [Pa]
    Tᵥ = diagn.physics.temp_virt[ij, end]       # virtual temperature at lowest model level [K]
    Tᵥ *= σ^(-κ)                                # lower to surface assuming dry adiabatic lapse rate
    ρ = pₛ/(R_dry*Tᵥ)                           # surface air density [kg/m³] from ideal gas law
    diagn.physics.surface_air_density[ij] = ρ   # store for surface temp/humidity fluxes

    # Fortran SPEEDY documentation eq. 52, 53, accumulate fluxes with +=
    V₀ = surface_wind_speed
    flux_u_upward = -ρ*drag*V₀*surface_u
    flux_v_upward = -ρ*drag*V₀*surface_v
    
    # convert fluxes to tendencies
    k = size(u_tend_grid, 2)                    # lowest model level index
    m = (planet=planet, geometry=geometry)
    diagn.tendencies.u_tend_grid[ij, end] += flux_to_tendency(k, flux_u_upward, pₛ, m)
    diagn.tendencies.v_tend_grid[ij, end] += flux_to_tendency(k, flux_v_upward, pₛ, m)

    return nothing
end