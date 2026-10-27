abstract type AbstractBoundaryLayerDrag <: AbstractParameterization end

# every boundary layer is expected to write these on every time step
variables(::AbstractBoundaryLayerDrag) = (
    ParameterizationVariable(:boundary_layer_drag_momentum, GridXY(), desc = "Boundary layer coefficient for momentum", units = "1"),
    ParameterizationVariable(:boundary_layer_drag_heat, GridXY(), desc = "Boundary layer drag coefficient for heat", units = "1"),
    ParameterizationVariable(:boundary_layer_drag_humidity, GridXY(), desc = "Boundary layer drag coefficient for humidity", units = "1"),
    ParameterizationVariable(:boundary_layer_drag_momentum, GridXY(), desc = "Boundary layer land drag coefficient for momentum", units = "1", namespace = :land),
    ParameterizationVariable(:boundary_layer_drag_heat, GridXY(), desc = "Boundary layer land drag coefficient for heat", units = "1", namespace = :land),
    ParameterizationVariable(:boundary_layer_drag_humidity, GridXY(), desc = "Boundary layer land drag coefficient for humidity", units = "1", namespace = :land),
    ParameterizationVariable(:boundary_layer_drag_momentum, GridXY(), desc = "Boundary layer ocean drag coefficient for momentum", units = "1", namespace = :ocean),
    ParameterizationVariable(:boundary_layer_drag_heat, GridXY(), desc = "Boundary layer ocean drag coefficient for heat", units = "1", namespace = :ocean),
    ParameterizationVariable(:boundary_layer_drag_humidity, GridXY(), desc = "Boundary layer ocean drag coefficient for humidity", units = "1", namespace = :ocean),
)

export ConstantDrag
"""Constant boundary layer drag coefficient. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct ConstantDrag{NF} <: AbstractBoundaryLayerDrag
    "[OPTION] Constant land drag coefficient [1]"
    @param land_drag::NF = 1.0e-3 (bounds = Nonnegative,)

    "[OPTION] Constant ocean drag coefficient [1]"
    @param ocean_drag::NF = 1.0e-3 (bounds = Nonnegative,)
end

Adapt.@adapt_structure ConstantDrag
ConstantDrag(SG::SpectralGrid; kwargs...) = ConstantDrag{SG.NF}(; kwargs...)
initialize!(::ConstantDrag, ::PrimitiveEquation) = nothing
@propagate_inbounds parameterization!(ij, vars, drag::ConstantDrag, boundary_layer, model) =
    boundary_layer_drag!(ij, vars, drag, model)

@propagate_inbounds function boundary_layer_drag!(ij, vars, scheme::ConstantDrag, model)
    land_fraction = model.land_sea_mask.land_fraction[ij]

    vars.parameterizations.boundary_layer_drag_momentum[ij] = land_fraction * scheme.land_drag + (1 - land_fraction) * scheme.ocean_drag
    vars.parameterizations.boundary_layer_drag_heat[ij] = land_fraction * scheme.land_drag + (1 - land_fraction) * scheme.ocean_drag
    vars.parameterizations.boundary_layer_drag_humidity[ij] = land_fraction * scheme.land_drag + (1 - land_fraction) * scheme.ocean_drag

    vars.parameterizations.land.boundary_layer_drag_momentum[ij] = scheme.land_drag
    vars.parameterizations.land.boundary_layer_drag_heat[ij] = scheme.land_drag
    vars.parameterizations.land.boundary_layer_drag_humidity[ij] = scheme.land_drag

    vars.parameterizations.ocean.boundary_layer_drag_momentum[ij] = scheme.ocean_drag
    vars.parameterizations.ocean.boundary_layer_drag_heat[ij] = scheme.ocean_drag
    vars.parameterizations.ocean.boundary_layer_drag_humidity[ij] = scheme.ocean_drag
    return nothing
end

export BulkRichardsonDrag

"""Boundary layer drag coefficient from the bulk Richardson number,
following Frierson, 2006, Journal of the Atmospheric Sciences.
$(TYPEDFIELDS)"""
@parameterized @kwdef struct BulkRichardsonDrag{NF} <: AbstractBoundaryLayerDrag
    "[OPTION] von Kármán constant [1]"
    @param von_Karman::NF = 0.4 (bounds = 0 .. 1,)

    "[OPTION] Critical Richardson number for stable mixing cutoff [1]"
    @param critical_Richardson::NF = 10 (bounds = Positive,)

    "[OPTION] Drag minimum to avoid zero surface fluxes in stable conditions [1]"
    @param drag_min::NF = 1.0e-5 (bounds = Nonnegative,)
end

Adapt.@adapt_structure BulkRichardsonDrag
BulkRichardsonDrag(SG::SpectralGrid; kwargs...) = BulkRichardsonDrag{SG.NF}(; kwargs...)
initialize!(::BulkRichardsonDrag, ::PrimitiveEquation) = nothing

# function barrier
@propagate_inbounds parameterization!(ij, vars, drag::BulkRichardsonDrag, boundary_layer, model) =
    boundary_layer_drag!(ij, vars, drag, model)

@propagate_inbounds function boundary_layer_drag!(
        ij,
        vars,
        drag::BulkRichardsonDrag,
        model
    )
    (; atmosphere, planet, orography, time_stepping, land_sea_mask) = model
    # Height z [m] of lowermost layer above ground
    surface = size(vars.dynamics.geopotential, 2)    # surface index = nlayers
    (; gravity) = planet
    land_fraction = land_sea_mask.land_fraction[ij]
    z = vars.dynamics.geopotential[ij, surface] / gravity - orography.orography[ij]

    # maximum drag Cmax from that height, stable conditions would decrease Cmax towards 0
    # Frierson 2006, eq (12)
    κ = drag.von_Karman

    # Get surface roughness length (computed by the surface_roughness parameterization)
    z₀M_land = vars.parameterizations.land.momentum_roughness[ij]
    z₀H_land = vars.parameterizations.land.heat_roughness[ij]
    z₀Q_land = vars.parameterizations.land.humidity_roughness[ij]

    z₀M_ocean = vars.parameterizations.ocean.momentum_roughness[ij]
    z₀H_ocean = vars.parameterizations.ocean.heat_roughness[ij]
    z₀Q_ocean = vars.parameterizations.ocean.humidity_roughness[ij]

    C_M_land = neutral_transfer_coefficient(z, z₀M_land, z₀M_land, κ)
    C_H_land = neutral_transfer_coefficient(z, z₀M_land, z₀H_land, κ)
    C_Q_land = neutral_transfer_coefficient(z, z₀M_land, z₀Q_land, κ)

    C_M_ocean = neutral_transfer_coefficient(z, z₀M_ocean, z₀M_ocean, κ)
    C_H_ocean = neutral_transfer_coefficient(z, z₀M_ocean, z₀H_ocean, κ)
    C_Q_ocean = neutral_transfer_coefficient(z, z₀M_ocean, z₀Q_ocean, κ)

    # bulk Richardson number at lowermost layer from Frierson, 2006, eq. (15)
    # they call it Ri_a = Ri here
    ΔΦ₀ = gravity * z     # geopotential height relative to surface
    Ri = bulk_richardson_surface(ij, ΔΦ₀, vars, atmosphere, drag, time_stepping)
    Ri_c = drag.critical_Richardson
    (; drag_min) = drag

    # clamp to get the cases, eq (12-14)
    # if Ri > Ri_c then C = 0
    # if Ri_c > Ri > 0 then = κ^2/log(z/z₀)^2 * (1-Ri/Ri_c)^2
    # if Ri_c < 0 then κ^2/log(z/z₀)^2
    Ri = clamp(Ri, 0, Ri_c)

    # unpack for readability
    (; land, ocean) = vars.parameterizations

    land.boundary_layer_drag_momentum[ij] = max(drag_min, C_M_land * (1 - Ri / Ri_c)^2)
    land.boundary_layer_drag_heat[ij] = max(drag_min, C_H_land * (1 - Ri / Ri_c)^2)
    land.boundary_layer_drag_humidity[ij] = max(drag_min, C_Q_land * (1 - Ri / Ri_c)^2)

    ocean.boundary_layer_drag_momentum[ij] = max(drag_min, C_M_ocean * (1 - Ri / Ri_c)^2)
    ocean.boundary_layer_drag_heat[ij] = max(drag_min, C_H_ocean * (1 - Ri / Ri_c)^2)
    ocean.boundary_layer_drag_humidity[ij] = max(drag_min, C_Q_ocean * (1 - Ri / Ri_c)^2)

    # For output reading
    vars.parameterizations.boundary_layer_drag_momentum[ij] =
        land_fraction * land.boundary_layer_drag_momentum[ij] + (1 - land_fraction) * ocean.boundary_layer_drag_momentum[ij]
    vars.parameterizations.boundary_layer_drag_heat[ij] =
        land_fraction * land.boundary_layer_drag_heat[ij] + (1 - land_fraction) * ocean.boundary_layer_drag_heat[ij]
    vars.parameterizations.boundary_layer_drag_humidity[ij] =
        land_fraction * land.boundary_layer_drag_humidity[ij] + (1 - land_fraction) * ocean.boundary_layer_drag_humidity[ij]

    return nothing
end

"""
$(TYPEDSIGNATURES)
Bulk transfer coefficient in the neutral limit, κ²/(log(z/z₀M) log(z/z₀X)), for the height `z` [m] of
the lowermost layer, the momentum roughness length `z₀M` [m] and the roughness length `z₀X` [m] of
momentum, heat or humidity. Should be z > z₀ (z = z₀ means infinitely high drag), so capped one order
of magnitude above the roughness length. A roughness length of exactly zero marks a surface type that
is absent at a grid point (the coefficient is then weighted by a zero land or ocean fraction). 1e-5 m
is used instead, as `log(z/0) = Inf` gives a zero coefficient but a NaN derivative (0 · ∞) in reverse
mode, e.g. with Reactant. Its value is irrelevant as long as it's positive, but it is only a replacement
for zero, not a lower bound: real roughness lengths below 1e-5 m are kept."""
@inline function neutral_transfer_coefficient(z, z₀M, z₀X, κ)
    z₀M = ifelse(z₀M > 0, z₀M, oftype(z₀M, 1.0e-5))
    z₀X = ifelse(z₀X > 0, z₀X, oftype(z₀X, 1.0e-5))
    zM = max(z, 10z₀M)
    zX = max(z, 10z₀X)
    return κ^2 / (log(zM / z₀M) * log(zX / z₀X))
end

"""
$(TYPEDSIGNATURES)
Calculate the bulk Richardson number following Frierson, 2006.
For vertical stability in the boundary layer."""
@propagate_inbounds function bulk_richardson_surface(ij, ΔΦ₀, vars, atmosphere, drag, time_stepping)
    cₚ = atmosphere.heat_capacity
    temp = get_prognostic_step(vars.grid.temperature, time_stepping, drag)
    NF = eltype(temp)
    surface = size(temp, 2)     # surface index = nlayers

    Vₛ = vars.parameterizations.surface_wind_speed[ij]
    T = temp[ij, surface]
    q = haskey(vars.grid, :humidity) ? get_prognostic_step(vars.grid.humidity, time_stepping, drag)[ij, surface] : zero(NF)
    Tᵥ = virtual_temperature(T, q, atmosphere)

    # bulk Richardson number at lowermost layer N from Frierson, 2006, eq. (15)
    Θ₀ = cₚ * Tᵥ          # virtual dry static energy at surface (z=0)
    Θ₁ = Θ₀ + ΔΦ₀       # virtual dry static energy at first model level (z=z)
    bulk_richardson = ΔΦ₀ * (Θ₁ - Θ₀) / (Θ₀ * Vₛ^2)
    return bulk_richardson
end
