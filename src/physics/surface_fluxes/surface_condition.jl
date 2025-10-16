abstract type AbstractSurfaceCondition <: AbstractParameterization end

export SurfaceCondition

@kwdef struct SurfaceCondition{NF} <: AbstractSurfaceCondition
    "Ratio of near-surface wind to lowest-level wind [1]"
    f_wind::NF = 0.95

    "Wind speed of sub-grid scale gusts [m/s]"
    gust_speed::NF = 5
end

SurfaceCondition(SG::SpectralGrid; kwargs...) = SurfaceCondition{SG.NF}(; kwargs...)
initialize!(::SurfaceCondition, ::PrimitiveEquation) = nothing

# function barrier
function parameterization!(ij, diagn, progn, sc::SurfaceCondition, model)
    surface_condition!(ij, diagn, sc, model)
end

function surface_condition!(ij, diagn, surface_condition::SurfaceCondition, model)

    (; f_wind, gust_speed) = surface_condition

    # Fortran SPEEDY documentation eq. 49 but use previous time step for numerical stability
    surface_u = f_wind*diagn.grid.u_grid_prev[ij, end] 
    surface_v = f_wind*diagn.grid.v_grid_prev[ij, end]

    # Fortran SPEEDY documentation eq. 50
    surface_wind_speed = sqrt(surface_u^2 + surface_v^2 + gust_speed^2)
    diagn.physics.surface_wind_speed[ij] = surface_wind_speed

    # Surface air density
    (; R_dry, κ) = model.atmosphere
    σ = model.geometry.σ_levels_full[k]
    pₛ = diagn.grid.pres_grid_prev[ij]          # surface pressure [Pa]
    Tᵥ = diagn.grid.temp_virt_grid[ij, end]     # virtual temperature at lowest model level [K]
    σ⁻ᵏ = σ^(-κ)                                # precalculate
    Tᵥ *= σ⁻ᵏ                                   # lower to surface assuming dry adiabatic lapse rate
    ρ = pₛ/(R_dry*Tᵥ)                           # surface air density [kg/m³] from ideal gas law
    diagn.physics.surface_air_density[ij] = ρ   # store for surface temp/humidity fluxes

    # Surface air temperature
    T = diagn.grid.temp_grid[ij, end]               # virtual temperature at lowest model level [K]
    T *= σ⁻ᵏ                                        # lower to surface assuming dry adiabatic lapse rate
    diagn.physics.surface_air_temperature[ij] = T   # store for surface temp/humidity fluxes
    return nothing
end