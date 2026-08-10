abstract type AbstractBoundaryLayer <: AbstractParameterization end

variables(::AbstractBoundaryLayer) = (
    ParameterizationVariable(:surface_wind_speed, Grid2D(), desc = "Surface wind speed", units = "m/s"),
    ParameterizationVariable(:surface_air_density, Grid2D(), desc = "Surface air density", units = "kg/m³"),
    ParameterizationVariable(:surface_air_temperature, Grid2D(), desc = "Surface air temperature", units = "K"),
)

export ConstantDrag
"""Constant boundary layer drag coefficient. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct ConstantDrag{NF} <: AbstractParameterization
    "[OPTION] Constant drag coefficient [1]"
    @param drag::NF = 1.0e-3
end

variables(::ConstantDrag) = (
    ParameterizationVariable(:boundary_layer_drag_momentum, Grid2D(), desc = "Boundary layer drag coefficient for momentum", units = "1"),
    ParameterizationVariable(:boundary_layer_drag_heat, Grid2D(), desc = "Boundary layer drag coefficient for heat", units = "1"),
    ParameterizationVariable(:boundary_layer_drag_humidity, Grid2D(), desc = "Boundary layer drag coefficient for humidity", units = "1"),
)

Adapt.@adapt_structure ConstantDrag
ConstantDrag(SG::SpectralGrid; kwargs...) = ConstantDrag{SG.NF}(; kwargs...)
initialize!(::ConstantDrag, ::PrimitiveEquation) = nothing
@propagate_inbounds parameterization!(ij, vars, drag::ConstantDrag, model) =
    boundary_layer_drag!(ij, vars, drag)

@propagate_inbounds function boundary_layer_drag!(ij, vars, scheme::ConstantDrag)
    vars.parameterizations.boundary_layer_drag[ij] = scheme.drag
    return nothing
end

export NeutralWindSpeed
"""Ocean-based neutral wind speed calculation from actual wind speed, 
derived from ERA5 data via symbolic regression. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct NeutralWindSpeed{NF} <: AbstractParameterization
    # Parameters for neutral wind calculation
    "[OPTION] symbolic regression parameter 1 [1]"
    @param c1::NF = -0.039317116

    "[OPTION] symbolic regression parameter 2 [1]"
    @param c2::NF = -2.9858496

    "[OPTION] symbolic regression parameter 3 [1]"
    @param c3::NF = 2.0046231e-10

    "[OPTION] symbolic regression parameter 4 [1]"
    @param c4::NF = 1.0768474

    "[OPTION] symbolic regression parameter 5 [1]"
    @param c5::NF = 0.20268184

    "[OPTION] symbolic regression parameter 6 [1]"
    @param c6::NF = 1.2684147

    "[OPTION] symbolic regression parameter 7 [1]"
    @param c7::NF = -0.94933933

    "[OPTION] symbolic regression parameter 8 [1]"
    @param c8::NF = 0.041551278

    "[OPTION] symbolic regression parameter 9 [1]"
    @param c9::NF = 5.8649142

    # Safety parameters
    "[OPTION] Minimum wind speed to avoid zero surface fluxes [m/s]"
    minimum_wind::NF = 1.0e-6

    "[OPTION] Minimum argument for logarithm to avoid log(0) [1]"
    log_safety::NF = 1.0e-8
end

variables(::NeutralWindSpeed) = (
    ParameterizationVariable(:neutral_wind_speed, Grid2D(), desc = "Neutral surface wind speed", units = "m/s"),
)


Adapt.@adapt_structure NeutralWindSpeed
NeutralWindSpeed(SG::SpectralGrid; kwargs...) = NeutralWindSpeed{SG.NF}(; kwargs...)
initialize!(::NeutralWindSpeed, ::PrimitiveEquation) = nothing

Base.show(io::IO, M::AbstractBoundaryLayer) = Base.show(io, M, values = false)

export BoundaryLayer
"""Composite type, containing surface roughness computation
and drag coefficient computation. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct BoundaryLayer{SC, NW, SR, D} <: AbstractBoundaryLayer
    @component surface_condition::SC
    @component neutral_wind_speed::NW
    @component surface_roughness::SR
    @component drag::D
end

Adapt.@adapt_structure BoundaryLayer
function BoundaryLayer(
        SG::SpectralGrid;
        neutral_wind_speed = NeutralWindSpeed(SG),
        surface_roughness = LearnedSurfaceRoughness(SG),
        surface_condition = SurfaceCondition(SG),
        drag = BulkRichardsonDrag(SG),
    )
    return BoundaryLayer(surface_condition, neutral_wind_speed, surface_roughness, drag)
end

function initialize!(BL::BoundaryLayer)
    initialize!(BL.neutral_wind_speed)
    initialize!(BL.surface_roughness)
    initialize!(BL.surface_condition)
    initialize!(BL.drag)
    return nothing
end

# variables of boundary layer are the union of surface roughness and drag variables
variables(BL::BoundaryLayer) = (
    variables(BL.neutral_wind_speed)...,
    variables(BL.surface_roughness)...,
    variables(BL.surface_condition)...,
    variables(BL.drag)...,
)

# just call the sub-paramterizations one after another
@propagate_inbounds function parameterization!(ij, vars, BL::BoundaryLayer, model)
    parameterization!(ij, vars, BL.neutral_wind_speed, model)
    parameterization!(ij, vars, BL.surface_roughness, model)
    parameterization!(ij, vars, BL.surface_condition, model)
    parameterization!(ij, vars, BL.drag, model)
    return nothing
end

export BulkRichardsonDrag

"""Boundary layer drag coefficient from the bulk Richardson number,
following Frierson, 2006, Journal of the Atmospheric Sciences.
$(TYPEDFIELDS)"""
@parameterized @kwdef struct BulkRichardsonDrag{NF} <: AbstractParameterization
    "[OPTION] von Kármán constant [1]"
    @param von_Karman::NF = 0.4 (bounds = 0 .. 1,)

    "[OPTION] Critical Richardson number for stable mixing cutoff [1]"
    @param critical_Richardson::NF = 10 (bounds = Positive,)

    "[OPTION] Drag minimum to avoid zero surface fluxes in stable conditions [1]"
    @param drag_min::NF = 1.0e-5 (bounds = Nonnegative,)
end

variables(::BulkRichardsonDrag) = (
    ParameterizationVariable(:boundary_layer_drag_momentum, Grid2D(), desc = "Boundary layer drag coefficient for momentum", units = "1"),
    ParameterizationVariable(:boundary_layer_drag_heat, Grid2D(), desc = "Boundary layer drag coefficient for heat", units = "1"),
    ParameterizationVariable(:boundary_layer_drag_humidity, Grid2D(), desc = "Boundary layer drag coefficient for humidity", units = "1"),
)

Adapt.@adapt_structure BulkRichardsonDrag
BulkRichardsonDrag(SG::SpectralGrid, kwargs...) = BulkRichardsonDrag{SG.NF}(; kwargs...)
initialize!(::BulkRichardsonDrag, ::PrimitiveEquation) = nothing

# function barrier
@propagate_inbounds parameterization!(ij, vars, drag::BulkRichardsonDrag, model) =
    boundary_layer_drag!(ij, vars, drag, model.atmosphere, model.planet, model.orography, model.time_stepping)

@propagate_inbounds function boundary_layer_drag!(
        ij,
        vars,
        drag::BulkRichardsonDrag,
        atmosphere,
        planet,
        orography,
        time_stepping,
    )

    # Height z [m] of lowermost layer above ground
    surface = size(vars.dynamics.geopotential, 2)    # surface index = nlayers
    (; gravity) = planet
    z = vars.dynamics.geopotential[ij, surface] / gravity - orography.orography[ij]

    # maximum drag Cmax from that height, stable conditions would decrease Cmax towards 0
    # Frierson 2006, eq (12)
    κ = drag.von_Karman

    # Get surface roughness length (computed by the surface_roughness parameterization)
    z₀M = vars.parameterizations.momentum_roughness[ij]
    z₀H = vars.parameterizations.heat_roughness[ij]
    z₀Q = vars.parameterizations.moisture_roughness[ij]

    function calc_drag_max(z₀)
        # should be z > z₀, z=z₀ means an infinitely high drag, choose one order higher than roughness length at least
        # 0 < z < z₀ doesn't make sense so cap here
        z_calc = max(z, 10z₀)
        return (κ / log(z_calc / z₀))^2
    end

    drag_max_momentum = calc_drag_max(z₀M)
    drag_max_heat = calc_drag_max(z₀H)
    drag_max_humidity = calc_drag_max(z₀Q)

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
    vars.parameterizations.boundary_layer_drag_momentum[ij] = max(drag_min, drag_max_momentum * (1 - Ri / Ri_c)^2)
    vars.parameterizations.boundary_layer_drag_heat[ij] = max(drag_min, drag_max_heat * (1 - Ri / Ri_c)^2)
    vars.parameterizations.boundary_layer_drag_humidity[ij] = max(drag_min, drag_max_humidity * (1 - Ri / Ri_c)^2)
    return nothing
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

@propagate_inbounds function parameterization!(ij, vars, nw::NeutralWindSpeed{NF}, model) where {NF}
    (; land_fraction) = model.land_sea_mask
    (land_fraction[ij] < 1) || return nothing # TODO train a land-based neutral wind speed parameterization
    return neutral_wind_speed(ij, vars, nw, model)
end

@propagate_inbounds function neutral_wind_speed(ij, vars, nw::NeutralWindSpeed{NF}, model) where {NF}
    (; surface_wind_speed) = vars.parameterizations
    (; surface_air_temperature) = vars.parameterizations

    sst = vars.prognostic.ocean.sea_surface_temperature[ij]
    t_diff = surface_air_temperature[ij] - sst # TODO: replace SST with ocean skin temperature
    ws_safe = max(surface_wind_speed[ij], nw.minimum_wind)
    log_arg = max(nw.c1 * ws_safe * t_diff + exp(t_diff), nw.log_safety)

    numerator = 2 * t_diff + nw.c8 * exp(t_diff) - nw.c3 * (nw.c4^surface_air_temperature[ij])
    denominator = t_diff * (log(log_arg) + nw.c2) + nw.c5 * (nw.c6^ws_safe) + nw.c9 * (ws_safe^nw.c7) + ws_safe

    vars.parameterizations.neutral_wind_speed[ij] = max(surface_wind_speed[ij] - (numerator / denominator), 0)

    return nothing
end
