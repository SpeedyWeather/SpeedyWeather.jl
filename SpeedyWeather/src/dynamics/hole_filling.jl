abstract type AbstractHoleFilling <: AbstractModelComponent end

# Hole filling for the specific humidity 
# Note that this is a relatively hacky solution to the problem of negative humidity values in the spectral 
# model. We might revist how humidity is handled in the model in the future, but for now this is a simple 
# way to ensure that the model does not produce unphysical negative humidity values.

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

"""Hole filling that sets negative grid-point humidity to zero, `q = max(q, 0)`, after every
transform to grid-point space, so that dynamics and parameterizations see non-negative humidity.
The spectral prognostic humidity keeps its negative values, where all tendencies see zero humidity.
With `fill_tendency``, the humidity the clipping adds at every grid point, `max(-q, 0)`, is stored
in `vars.dynamics.humidity_hole_fill` and added as a tendency in the next time step, so that the spectral
humidity is stepped forward from the clipped humidity that the tendencies were computed from (with
`Leapfrog` from the clipped previous step). Without it holes in the spectral humidity can deepen and
spread, as with the additive noise of `StochasticLeapfrog` in the dry upper layers. Not water conserving,
the water to fill the holes is added. Construct with

    ClipNegatives(spectral_grid; fill_tendency = nothing)

`fill_tendency = true` or `false` adds the fill as tendency or not, `nothing` (default) adds it
only with a `StochasticLeapfrog` time stepper."""
struct ClipNegatives{F} <: AbstractHoleFilling end     # F = fill_tendency: true, false or nothing (automatic)

ClipNegatives(; fill_tendency::Union{Nothing, Bool} = nothing) = ClipNegatives{fill_tendency}()
ClipNegatives(::SpectralGrid; kwargs...) = ClipNegatives(; kwargs...)

"""$(TYPEDSIGNATURES)
Whether the hole filling adds its fill as tendency with `time_stepping`: as chosen by `fill_tendency`,
with `fill_tendency = nothing` only for `StochasticLeapfrog` (method defined there). Known at compile time."""
adds_fill_tendency(::ClipNegatives{true}, time_stepping) = true
adds_fill_tendency(::ClipNegatives{false}, time_stepping) = false
adds_fill_tendency(::ClipNegatives{nothing}, time_stepping) = false

function variables(hole_filling::ClipNegatives, model::AbstractModel)
    adds_fill_tendency(hole_filling, model.time_stepping) || return ()
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
Clip negative grid-point humidity of the current step, with the fill tendency also store the
humidity added in `vars.dynamics.humidity_hole_fill`."""
function hole_filling!(vars::Variables, hole_filling::ClipNegatives, model::AbstractModel)
    (; time_stepping, spectral_transform) = model
    humid_grid = get_prognostic_step(vars.grid.humidity, time_stepping, spectral_transform)
    if adds_fill_tendency(hole_filling, time_stepping)      # resolved at compile time
        humid_fill = get_prognostic_step(vars.dynamics.humidity_hole_fill, time_stepping, spectral_transform)
        clip_negatives!(humid_grid, humid_fill)
    else
        hole_filling!(humid_grid, hole_filling)
    end
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
With the fill tendency (see `ClipNegatives`), add the humidity of the hole filling,
`vars.dynamics.humidity_hole_fill`, divided by the time step to the grid humidity tendency. With `Leapfrog`, which steps forward from the previous step, the fill
of the previous step is used, so with Δt of this step (Δt/2 and Δt in the leapfrog start, then 2Δt)
the spectral humidity is stepped forward from the clipped previous humidity. Called after the
parameterizations (e.g. stochastic physics should not perturb the fill) and before the dynamics,
which scales all tendencies."""
function hole_filling_tendency!(vars::Variables, hole_filling::ClipNegatives, model::AbstractModel)
    (; time_stepping) = model
    adds_fill_tendency(hole_filling, time_stepping) || return nothing
    Δt = time_step(time_stepping, vars.prognostic.clock)    # time step of this step [s]
    humid_tend_grid = get_tendency_step(vars.tendencies.grid.humidity, time_stepping, hole_filling)
    humid_fill = get_prognostic_step(vars.dynamics.humidity_hole_fill, time_stepping, hole_filling)
    humid_tend_grid.data .+= humid_fill.data ./ Δt    # .data: a Field broadcast fails Enzyme in reverse mode
    return nothing
end
