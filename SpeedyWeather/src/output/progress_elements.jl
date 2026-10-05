"""
    SpeedyWeather.ProgressElements

Elements of the progress line printed by [`Feedback`](@ref SpeedyWeather.Feedback) while a simulation is running. Contains
the generic elements of ProgressMeter.jl (`Description`, `Percentage`, `Bar`, `ETA`, `Speed`,
`ElapsedTime`, `Counter`, `Colored`) and the ones specific to SpeedyWeather (`SimulationTime`,
`SimulationSpeed`, `MaximumWindSpeed`, `TemperatureRange`, `VerticalCourantNumber`). Pass a tuple of
them as `Feedback(elements = ...)`, [`default_elements`](@ref)
returns the default ones, e.g.

```julia
feedback = Feedback(elements = (ProgressElements.default_elements()..., ProgressElements.VerticalCourantNumber()))
```"""
module ProgressElements

using DocStringExtensions: TYPEDSIGNATURES
import Dates
import Printf: @sprintf
import ProgressMeter
import ProgressMeter.Elements: AbstractProgressElement, print_element, Colored, Description,
    Percentage, Bar, ETA, Speed, ElapsedTime, Counter

export AbstractProgressElement, print_element, Colored, Description, Percentage, Bar, ETA, Speed,
    ElapsedTime, Counter, SimulationTime, SimulationSpeed, MaximumWindSpeed, TemperatureRange,
    VerticalCourantNumber, default_elements

# Elements of the progress line. They are created unbound, e.g. `VerticalCourantNumber()`, and
# `bind_element` returns a copy that holds what it needs from the simulation (`Variables`, time
# step, ...) in `initialize!(::Feedback, ...)`. `print_element` is only called when the progress
# meter is redrawn, so diagnostics are only computed when they are displayed.

abstract type AbstractSimulationElement <: AbstractProgressElement end

# elements that do not need anything from the simulation are their own bound version
bind_element(element, vars, model) = element

# bound elements hold the `Variables`, don't print them
Base.show(io::IO, element::AbstractSimulationElement) = print(io, nameof(typeof(element)), "()")

"""$(TYPEDSIGNATURES)
Progress line element that shows the current simulation date."""
struct SimulationTime{C} <: AbstractSimulationElement
    clock::C
end
SimulationTime() = SimulationTime(nothing)
bind_element(::SimulationTime, vars, model) = SimulationTime(vars.prognostic.clock)
print_element(element::SimulationTime, p) = string(Dates.Date(element.clock.time))
print_element(::SimulationTime{Nothing}, p) = ""

"""$(TYPEDSIGNATURES)
Progress line element that shows the simulation speed, e.g. in simulated years per day,
which needs the time step `Δt` [s]. `separator` is printed in front."""
@kwdef struct SimulationSpeed <: AbstractSimulationElement
    separator::String = ", "
    Δt::Float64 = 0.0
end
bind_element(element::SimulationSpeed, vars, model) = SimulationSpeed(element.separator, Float64(model.time_stepping.Δt))

function print_element(element::SimulationSpeed, p)
    sec_per_iter = (p.tcurrent - p.tinit) / max(1, p.counter - p.start)
    return element.separator * speedstring(sec_per_iter, element.Δt)
end

"""$(TYPEDSIGNATURES)
Format a simulation speed from `sec_per_iter` wallclock seconds per time step of `dt_in_sec` seconds
to simulated days/years/... per day."""
function speedstring(sec_per_iter, dt_in_sec)
    (sec_per_iter == Inf || dt_in_sec <= 0) && return "N/A  days/day"

    sim_time_per_time = dt_in_sec / sec_per_iter

    for (divideby, unit) in (
            (365 * 1_000, "millenia"),
            (365, "years"),
            (1, "days"),
            (1 / 24, "hours"),
        )
        if (sim_time_per_time / divideby) > 2
            return @sprintf "%5.2f %2s/day" (sim_time_per_time / divideby) unit
        end
    end
    return "<2 hours/days"
end

"""$(TYPEDSIGNATURES)
Progress line element that shows the maximum absolute zonal wind `|u|` [m/s] of the current grid-space `u`.
Shows nothing if the model has no `u`."""
struct MaximumWindSpeed{V} <: AbstractSimulationElement
    vars::V
end
MaximumWindSpeed() = MaximumWindSpeed(nothing)
bind_element(::MaximumWindSpeed, vars, model) = MaximumWindSpeed(vars)
print_element(::MaximumWindSpeed{Nothing}, p) = ""

function print_element(element::MaximumWindSpeed, p)
    hasproperty(element.vars.grid, :u) || return ""
    umin, umax = extrema(element.vars.grid.u)
    return @sprintf ", %3d m/s" max(abs(umin), abs(umax))
end

"""$(TYPEDSIGNATURES)
Progress line element that shows the range of the grid-space temperature in ˚C.
Shows nothing if the model has no `temperature`."""
struct TemperatureRange{V} <: AbstractSimulationElement
    vars::V
end
TemperatureRange() = TemperatureRange(nothing)
bind_element(::TemperatureRange, vars, model) = TemperatureRange(vars)
print_element(::TemperatureRange{Nothing}, p) = ""

function print_element(element::TemperatureRange, p)
    hasproperty(element.vars.grid, :temperature) || return ""
    tmin, tmax = extrema(element.vars.grid.temperature)
    return @sprintf ", [%4d, %4d] ˚C" tmin - 273.15f0 tmax - 273.15f0
end

"""$(TYPEDSIGNATURES)
Progress line element that shows the maximum vertical Courant number
`max(|σ̇| Δt / Δσ)` over all grid points and layers, with σ̇ the vertical velocity in σ coordinates
at the layer interfaces, Δσ the layer thickness and Δt the time step. Shows nothing for models
without vertical velocity `vars.dynamics.w`. Not part of the default layout, add it with

```julia
Feedback(elements = (ProgressElements.default_elements()..., ProgressElements.VerticalCourantNumber()))
```"""
struct VerticalCourantNumber{V, T} <: AbstractSimulationElement
    vars::V
    Δσ::T
    Δt::Float64
end
VerticalCourantNumber() = VerticalCourantNumber(nothing, nothing, 0.0)

function bind_element(::VerticalCourantNumber, vars, model)
    hasproperty(model, :geometry) || return VerticalCourantNumber()
    return VerticalCourantNumber(vars, model.geometry.σ_levels_thick, Float64(model.time_stepping.Δt))
end
print_element(::VerticalCourantNumber{Nothing}, p) = ""

function print_element(element::VerticalCourantNumber, p)
    hasproperty(element.vars.dynamics, :w) || return ""
    scale = element.vars.prognostic.scale[]     # divergence, hence w, is scaled by the radius in the dynamical core
    return @sprintf ", Cᵥ = %.2f" vertical_courant_number(element.vars.dynamics.w, element.Δσ, element.Δt / scale)
end

"""$(TYPEDSIGNATURES)
Maximum vertical Courant number of layer `k`, `max(|σ̇ₖ₊₁/₂|, |σ̇ₖ₋₁/₂|) Δt / Δσₖ`, over all layers
and grid points. `w` is the vertical velocity at the layer interfaces `k+1/2` (zero at the surface),
`Δσ` the layer thickness. Pass `Δt / scale` if `w` is radius-scaled as in the dynamical core."""
function vertical_courant_number(w, Δσ, Δt)
    # maximum |σ̇| per interface k+1/2 reduced on the device, then loop over the few layers on the CPU
    w_max = vec(Array(maximum(abs, w.data, dims = 1)))
    Δσ_cpu = Array(Δσ)
    courant = zero(eltype(w_max))
    for k in eachindex(Δσ_cpu)
        w_above = k > 1 ? w_max[k - 1] : zero(courant)     # σ̇ = 0 at the top k = 1/2
        courant = max(courant, max(w_max[k], w_above) / Δσ_cpu[k])
    end
    return Δt * courant
end

"""$(TYPEDSIGNATURES)
Default elements of the progress line of a `Feedback`: description, percentage, bar, ETA and in
parentheses the simulation date, simulation speed, maximum wind speed and temperature range."""
default_elements() = (
    Description(), Percentage(), Bar(), ETA(), " (",
    SimulationTime(), SimulationSpeed(), MaximumWindSpeed(), TemperatureRange(), ")",
)

end # module ProgressElements
