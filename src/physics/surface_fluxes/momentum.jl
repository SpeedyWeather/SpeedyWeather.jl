export NoSurfaceWind
struct NoSurfaceWind <: AbstractSurfaceWind end
NoSurfaceWind(::SpectralGrid) = NoSurfaceWind()
initialize!(::NoSurfaceWind, ::PrimitiveEquation) = nothing
surface_wind_stress!(::ColumnVariables, ::NoSurfaceWind, ::PrimitiveEquation) = nothing

export SurfaceWind
Base.@kwdef struct SurfaceWind{NF<:AbstractFloat} <: AbstractSurfaceWind
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

function surface_wind_stress!(  column::ColumnVariables,
                                surface_wind::SurfaceWind,
                                model::PrimitiveEquation)

    (; land_fraction) = column
    (; f_wind, V_gust, drag_land, drag_sea) = surface_wind

    # SPEEDY documentation eq. 49, but use previous time step for numerical stability
    column.surface_u = f_wind*column.u[end] 
    column.surface_v = f_wind*column.v[end]
    (; surface_u, surface_v) = column

    # SPEEDY documentation eq. 50
    column.surface_wind_speed = sqrt(surface_u^2 + surface_v^2 + V_gust^2)

    # drag coefficient either from SurfaceEvaporation or from a central drag coefficient
    drag_sea, drag_land = surface_wind.use_boundary_layer_drag ?
                                (column.boundary_layer_drag, column.boundary_layer_drag) : 
                                (drag_sea, drag_land)
    
    # surface wind stress: quadratic drag, fractional land-sea mask
    ρ = column.surface_air_density
    V₀ = column.surface_wind_speed
    drag = land_fraction*drag_land + (1-land_fraction)*drag_sea

    # SPEEDY documentation eq. 52, 53, accumulate fluxes with +=
    column.flux_u_upward[end] -= ρ*drag*V₀*surface_u
    column.flux_v_upward[end] -= ρ*drag*V₀*surface_v
    
    return nothing
end