abstract type AbstractHoleFilling <: AbstractModelComponent end

# HOLE FILLING
# `hole_filling!(vars, hole_filling, model)` is called in `transform!` after the spectral prognostic
# humidity of the current step has been transformed to the grid and changes that grid humidity only.
# Everything that reads humidity afterwards (dynamics, parameterizations) uses the filled grid field,
# but the time stepping integrates the spectral humidity. A hole filling that changes the grid field
# therefore also adds that change as a tendency in the next time step (`hole_filling_tendency!`), so
# that the spectral humidity receives the moisture that the tendencies assume.

# no hole filling
hole_filling!(::Variables, ::Nothing, ::AbstractModel) = nothing
hole_filling_tendency!(::Variables, ::Nothing, ::AbstractModel) = nothing

# dispatch over model, only PrimitiveWet has humidity to fill
hole_filling_tendency!(vars::Variables, model::PrimitiveWet) =
    hole_filling_tendency!(vars, model.hole_filling, model)
hole_filling_tendency!(::Variables, ::AbstractModel) = nothing

# Fallbacks for custom hole fillings that only define `hole_filling!(humid_grid, hole_filling, model)`:
# fill the grid humidity of the current step, without a tendency the spectral humidity doesn't see the fill
function hole_filling!(vars::Variables, hole_filling::AbstractHoleFilling, model::AbstractModel)
    humid_grid = get_prognostic_step(vars.grid.humidity, model.time_stepping, model.spectral_transform)
    hole_filling!(humid_grid, hole_filling, model)
    return nothing
end

hole_filling_tendency!(::Variables, ::AbstractHoleFilling, ::AbstractModel) = nothing

export ClipNegatives

"""Hole filling that sets negative grid-point humidity to zero, `q = max(q, 0)`.
The humidity this adds at every grid point, `max(-q, 0)`, is stored in
`vars.dynamics.humidity_hole_fill` and added as a tendency in the next time step,
so that the spectral prognostic humidity is stepped forward from the clipped humidity
that the tendencies were computed from (with `Leapfrog` from the clipped previous step).
Without that tendency the spectral humidity would keep (and accumulate) negative values
where all tendencies see zero humidity. Not water conserving, the water to fill the
holes is added."""
struct ClipNegatives <: AbstractHoleFilling end
ClipNegatives(::SpectralGrid) = ClipNegatives()

function variables(::ClipNegatives, model::AbstractModel)
    # as many steps as the grid humidity: with Leapfrog the fill of the current step is
    # computed in transform! but the fill of the previous step is added as tendency
    pg = get_nsteps(model.time_stepping, model).prognostic_grid
    return (
        DynamicsVariable(
            :humidity_hole_fill, GridXYZT(pg), units = "kg/kg",
            desc = "Humidity added by the hole filling, clipped minus unclipped humidity"
        ),
    )
end

"""$(TYPEDSIGNATURES)
Clip negative grid-point humidity of the current step and store the humidity added
in `vars.dynamics.humidity_hole_fill`."""
function hole_filling!(vars::Variables, ::ClipNegatives, model::AbstractModel)
    (; time_stepping, spectral_transform) = model
    humid_grid = get_prognostic_step(vars.grid.humidity, time_stepping, spectral_transform)
    humid_fill = get_prognostic_step(vars.dynamics.humidity_hole_fill, time_stepping, spectral_transform)
    clip_negatives!(humid_grid, humid_fill)
    return nothing
end

"""$(TYPEDSIGNATURES)
Clip negative values of `humid_grid` in place and write what was added into `humid_fill`."""
function clip_negatives!(humid_grid::AbstractField, humid_fill::AbstractField)
    # broadcast over .data, a Field broadcast fails Enzyme's type analysis in reverse mode
    q = humid_grid.data
    zero_humidity = zero(eltype(q))
    humid_fill.data .= max.(.-q, zero_humidity)     # fill = clipped - unclipped ≥ 0
    q .= max.(q, zero_humidity)
    return nothing
end

# clip a grid field without storing the fill
hole_filling!(humid_grid::AbstractField, hole_filling::ClipNegatives, ::AbstractModel) =
    hole_filling!(humid_grid, hole_filling)

function hole_filling!(humid_grid::AbstractField, ::ClipNegatives)
    zero_humidity = zero(eltype(humid_grid))
    return @. humid_grid = max(humid_grid, zero_humidity)
end

"""$(TYPEDSIGNATURES)
Add the humidity of the hole filling, `vars.dynamics.humidity_hole_fill`, divided by the time step
to the grid humidity tendency. With `Leapfrog`, which steps forward from the previous step, the fill
of the previous step is used, so with Δt of this step (Δt/2 and Δt in the leapfrog start, then 2Δt)
the spectral humidity is stepped forward from the clipped previous humidity. Called after the
parameterizations (e.g. stochastic physics should not perturb the fill) and before the dynamics,
which scales all tendencies."""
function hole_filling_tendency!(vars::Variables, hole_filling::ClipNegatives, model::AbstractModel)
    (; time_stepping) = model
    Δt = time_step(time_stepping, vars.prognostic.clock)    # time step of this step [s]
    humid_tend_grid = get_tendency_step(vars.tendencies.grid.humidity, time_stepping, hole_filling)
    humid_fill = get_prognostic_step(vars.dynamics.humidity_hole_fill, time_stepping, hole_filling)
    humid_tend_grid.data .+= humid_fill.data ./ Δt    # .data: a Field broadcast fails Enzyme in reverse mode
    return nothing
end
