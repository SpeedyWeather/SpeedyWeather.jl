abstract type AbstractSurfaceCondition <: AbstractParameterization end

export SurfaceCondition

"""Surface condition parameterization that calculates near-surface atmospheric
variables needed for surface flux calculations. Computes surface wind speed
including sub-grid scale gusts, surface air density, and surface air temperature
by extrapolating from the lowest model level to the surface using standard
atmospheric relationships. Fields are $(TYPEDFIELDS)"""
@kwdef struct SurfaceCondition{NF} <: AbstractSurfaceCondition
    "[OPTION] Ratio of near-surface wind to lowest-level wind [1]"
    wind_slowdown::NF = 0.95

    "[OPTION] Wind speed of sub-grid scale gusts [m/s]"
    gust_speed::NF = 5
end

SurfaceCondition(SG::SpectralGrid; kwargs...) = SurfaceCondition{SG.NF}(; kwargs...)
initialize!(::SurfaceCondition, ::PrimitiveEquation) = nothing

# function barrier
function parameterization!(ij, diagn, progn, sc::SurfaceCondition, model)
    surface_condition!(ij, diagn, sc, model)
end

function surface_condition!(ij, diagn, surface_condition::SurfaceCondition, model)

    (; wind_slowdown, gust_speed) = surface_condition

    # Fortran SPEEDY documentation eq. 49 but use previous time step for numerical stability
    surface_u = wind_slowdown*diagn.grid.u_grid_prev[ij, end] 
    surface_v = wind_slowdown*diagn.grid.v_grid_prev[ij, end]

    # Fortran SPEEDY documentation eq. 50
    surface_wind_speed = sqrt(surface_u^2 + surface_v^2 + gust_speed^2)
    diagn.physics.surface_wind_speed[ij] = surface_wind_speed

    # Surface air density
    (; R_dry, κ) = model.atmosphere
    σ = model.geometry.σ_levels_full[end]           # σ vertical coordinate at lowest model level
    pₛ = diagn.grid.pres_grid_prev[ij]              # surface pressure [Pa]
    Tᵥ = diagn.grid.temp_virt_grid[ij, end]         # virtual temperature at lowest model level [K]
    σ⁻ᵏ = σ^(-κ)                                    # precalculate
    Tᵥ *= σ⁻ᵏ                                       # lower to surface assuming dry adiabatic lapse rate
    ρ = pₛ/(R_dry*Tᵥ)                               # surface air density [kg/m³] from ideal gas law
    diagn.physics.surface_air_density[ij] = ρ       # store for surface temp/humidity fluxes

    # Surface air temperature
    # TODO use _prev and add implicit.temperature_profile!
    T = diagn.grid.temp_grid[ij, end]               # virtual temperature at lowest model level [K]
    T *= σ⁻ᵏ                                        # lower to surface assuming dry adiabatic lapse rate
    diagn.physics.surface_air_temperature[ij] = T   # store for surface temp/humidity fluxes
    return nothing
end

Adapt.@adapt_structure SurfaceCondition
