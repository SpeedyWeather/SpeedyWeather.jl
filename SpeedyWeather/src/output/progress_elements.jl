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
import ..SpeedyWeather: variables, ScratchVariable, Vertical1D
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
Progress line element that shows an estimate of the maximum vertical Courant number
`max(|σ̇| Δt / Δσ)` over all grid points and layers, with σ̇ the vertical velocity in σ coordinates
at the layer interfaces, Δσ the layer thickness and Δt the time step. It is meant to indicate how
close the simulation is to the vertical stability limit, not to be exact, see
[`vertical_courant_number!`](@ref). Shows nothing for models without vertical velocity
`vars.dynamics.w`. Not part of the default layout, add it with

```julia
Feedback(elements = (ProgressElements.default_elements()..., ProgressElements.VerticalCourantNumber()))
```"""
struct VerticalCourantNumber{V, W, C, S} <: AbstractSimulationElement
    vars::V
    w_max::W            # 1×nlayers view on the vertical scratch vector, maximum |σ̇| per interface
    w_max_cpu::C        # copy of w_max on the CPU
    Δσ::S               # layer thickness on the CPU
    Δt::Float64
end
VerticalCourantNumber() = VerticalCourantNumber(nothing, nothing, nothing, nothing, 0.0)

# vertical scratch vector the horizontal maximum of |σ̇| is reduced into, avoids allocating on redraw
variables(::VerticalCourantNumber) = (
    ScratchVariable(:vertical_velocity_maximum, Vertical1D(), desc = "Maximum |σ̇| per layer interface for the vertical Courant number", units = "1/s"),
)

function bind_element(::VerticalCourantNumber, vars, model)
    hasproperty(vars.dynamics, :w) || return VerticalCourantNumber()
    # reshape once here (it allocates a wrapper) to reduce over the horizontal of w (npoints × nlayers)
    w_max = reshape(vars.scratch.vertical_velocity_maximum, 1, :)
    w_max_cpu = zeros(eltype(w_max), length(w_max))
    Δσ = Array(model.geometry.σ_levels_thick)
    return VerticalCourantNumber(vars, w_max, w_max_cpu, Δσ, Float64(model.time_stepping.Δt))
end
print_element(::VerticalCourantNumber{Nothing}, p) = ""

function print_element(element::VerticalCourantNumber, p)
    (; vars, w_max, w_max_cpu, Δσ, Δt) = element
    scale = vars.prognostic.scale[]     # divergence, hence w, is scaled by the radius in the dynamical core
    return @sprintf ", Cᵥ = %.2f" vertical_courant_number!(w_max, w_max_cpu, vars.dynamics.w, Δσ, Δt / scale)
end

"""$(TYPEDSIGNATURES)
Estimate of the maximum vertical Courant number `|σ̇ₖ₊₁/₂| Δt / Δσₖ` over all layers `k` and grid
points, without allocations. `w` is the vertical velocity σ̇ at the layer interfaces `k+1/2`, `Δσ` the
layer thickness (on the CPU). `w_max` (1×nlayers, on the device of `w`) and `w_max_cpu` (on the CPU)
are preallocated buffers. Pass `Δt / scale` if `w` is radius-scaled as in the dynamical core."""
function vertical_courant_number!(w_max, w_max_cpu, w, Δσ, Δt)
    # maximum |σ̇| per interface k+1/2 over the horizontal, reduced on the device into the scratch
    # vector, starting from zero (init = false) as |σ̇| ≥ 0
    fill!(w_max, 0)
    maximum!(abs, w_max, w.data; init = false)
    copyto!(w_max_cpu, w_max)           # only nlayers values, loop over them on the CPU

    courant = zero(eltype(w_max_cpu))
    for k in eachindex(w_max_cpu, Δσ)
        # Only the interface below layer k (k+1/2) is used, not also the one above (k-1/2).
        # σ levels vary smoothly so neighbouring layers have similar Δσ, and this is only meant to
        # indicate how close the simulation is to the vertical stability limit, not to be exact.
        # (σ̇ = 0 at the surface, so the bottom layer is only seen through the layer above it.)
        courant = max(courant, w_max_cpu[k] / Δσ[k])
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
