"""
    SpeedyWeather.ProgressElements

Elements of the progress lines printed by [`Feedback`](@ref SpeedyWeather.Feedback) while a simulation is running. Contains
the generic elements of `SpeedyWeatherInternals.ProgressMeters` (`Percentage`, `Bar`, `ETA`,
`ElapsedTime`, `Counter`, `NewLine`) and the ones specific to SpeedyWeather (`SimulationTime`,
`SimulationSpeed`, `MaximumWindSpeed`, `TemperatureRange`, `VerticalCourantNumber`, `Precipitation`,
`EnergyImbalance`). Pass a tuple of them as `Feedback(elements = ...)`, [`default_elements`](@ref)
returns the default ones, e.g. for a single line without the bar

```julia
(; Percentage, ETA, SimulationTime, SimulationSpeed, MaximumWindSpeed) = ProgressElements
feedback = Feedback(elements = (Percentage(), ETA(), SimulationTime(), SimulationSpeed(), MaximumWindSpeed()))
```"""
module ProgressElements

using DocStringExtensions: TYPEDSIGNATURES
import Dates
import Printf: @sprintf
import LinearAlgebra: dot
import SpeedyWeatherInternals.ProgressMeters: ProgressMeters, Progress, AbstractProgressElement,
    print_element, Percentage, Bar, ETA, ElapsedTime, Counter, NewLine
import ..SpeedyWeather: variables, ScratchVariable, Vertical1D, ReactantDevice

export AbstractProgressElement, print_element, Percentage, Bar, ETA, ElapsedTime, Counter, NewLine,
    AbstractSimulationElement, AbstractDiagnosticElement, OnReactant, SimulationTime, SimulationSpeed,
    MaximumWindSpeed, TemperatureRange, VerticalCourantNumber, Precipitation, EnergyImbalance, default_elements

# `Feedback` passes the variables and the model on to the elements on every redraw of the progress
# lines, `print_element(element, progress, vars, model)`, so the elements are created once and hold
# no data. Diagnostics are therefore only computed when they are displayed. Elements print only their
# value, the progress meter puts its separator between them, and an element that does not apply to
# the model (e.g. `TemperatureRange` without temperature) prints "", which leaves it out.

"""Supertype of SpeedyWeather's elements of the progress lines, shown with
`print_element(element, progress, vars, model)`."""
abstract type AbstractSimulationElement <: AbstractProgressElement end

# without the simulation (e.g. a progress meter on its own) there's nothing to show
print_element(::AbstractSimulationElement, p::Progress) = ""

"""Supertype of the elements showing a diagnostic of the simulation (wind speed, temperature, ...),
computed with `print_diagnostic(element, progress, vars, model)`. Not shown with Reactant (computing
them falls back to scalar indexing there), unless wrapped in [`OnReactant`](@ref)."""
abstract type AbstractDiagnosticElement <: AbstractSimulationElement end

print_element(element::AbstractDiagnosticElement, p::Progress, vars, model) =
    model.architecture isa ReactantDevice ? "" : print_diagnostic(element, p, vars, model)

"""$(TYPEDSIGNATURES)
Shows the diagnostic `element` also with Reactant, e.g. `OnReactant(MaximumWindSpeed())`."""
struct OnReactant{E <: AbstractDiagnosticElement} <: AbstractSimulationElement
    element::E
end
print_element(element::OnReactant, p::Progress, vars, model) = print_diagnostic(element.element, p, vars, model)
variables(element::OnReactant) = variables(element.element)

"""$(TYPEDSIGNATURES)
Progress line element that shows the current simulation date."""
struct SimulationTime <: AbstractSimulationElement end
print_element(::SimulationTime, p::Progress, vars, model) = string(Dates.Date(vars.prognostic.clock.time))

"""$(TYPEDSIGNATURES)
Progress line element that shows the simulation speed, e.g. in simulated years per day."""
struct SimulationSpeed <: AbstractSimulationElement end

function print_element(::SimulationSpeed, p::Progress, vars, model)
    sec_per_iter = (p.tcurrent - p.tinit) / max(1, p.counter)
    return speedstring(sec_per_iter, Float64(model.time_stepping.Δt))
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
Left out if the model has no `u`."""
struct MaximumWindSpeed <: AbstractDiagnosticElement end

function print_diagnostic(::MaximumWindSpeed, p::Progress, vars, model)
    hasproperty(vars.grid, :u) || return ""
    umin, umax = extrema(vars.grid.u)
    return @sprintf "%d m/s" max(abs(umin), abs(umax))
end

"""$(TYPEDSIGNATURES)
Progress line element that shows the range of the grid-space temperature in ˚C.
Left out if the model has no `temperature`."""
struct TemperatureRange <: AbstractDiagnosticElement end

function print_diagnostic(::TemperatureRange, p::Progress, vars, model)
    hasproperty(vars.grid, :temperature) || return ""
    tmin, tmax = extrema(vars.grid.temperature)
    return @sprintf "[%d, %d] ˚C" tmin - 273.15f0 tmax - 273.15f0
end

"""$(TYPEDSIGNATURES)
Progress line element that shows an estimate of the maximum vertical Courant number
`max(|σ̇| Δt / Δσ)` over all grid points and layers, with σ̇ the vertical velocity in σ coordinates
at the layer interfaces, Δσ the layer thickness and Δt the time step. It is meant to indicate how
close the simulation is to the vertical stability limit, not to be exact, see
[`vertical_courant_number`](@ref). Left out for models without vertical velocity `vars.dynamics.w`."""
struct VerticalCourantNumber <: AbstractDiagnosticElement end

# vertical scratch vector the horizontal maximum of |σ̇| is reduced into, avoids allocating on the device
variables(::VerticalCourantNumber) = (
    ScratchVariable(:vertical_velocity_maximum, Vertical1D(), desc = "Maximum |σ̇| per layer interface for the vertical Courant number", units = "1/s"),
)

function print_diagnostic(::VerticalCourantNumber, p::Progress, vars, model)
    hasproperty(vars.dynamics, :w) || return ""
    w_max = reshape(vars.scratch.vertical_velocity_maximum, 1, :)   # 1×nlayers to reduce over the horizontal
    Δt = model.time_stepping.Δt / vars.prognostic.scale[]           # w is radius-scaled in the dynamical core
    return @sprintf "Cᵥ = %.2f" vertical_courant_number(w_max, vars.dynamics.w, model.geometry.σ_levels_thick, Δt)
end

"""$(TYPEDSIGNATURES)
Estimate of the maximum vertical Courant number `|σ̇ₖ₊₁/₂| Δt / Δσₖ` over all layers `k` and grid
points. `w` is the vertical velocity σ̇ at the layer interfaces `k+1/2`, `Δσ` the layer thickness, `w_max`
a 1×nlayers buffer on the device of `w`. Only `nlayers` values are copied to the CPU. Pass `Δt / scale`
if `w` is radius-scaled as in the dynamical core."""
function vertical_courant_number(w_max, w, Δσ, Δt)
    # maximum |σ̇| per interface k+1/2 over the horizontal, reduced on the device into the scratch
    # vector, starting from zero (init = false) as |σ̇| ≥ 0
    fill!(w_max, 0)
    maximum!(abs, w_max, w.data; init = false)
    w_max_cpu, Δσ_cpu = Array(w_max), Array(Δσ)     # only nlayers values, loop over them on the CPU

    courant = zero(eltype(w_max_cpu))
    for k in eachindex(w_max_cpu, Δσ_cpu)
        # Only the interface below layer k (k+1/2) is used, not also the one above (k-1/2).
        # σ levels vary smoothly so neighbouring layers have similar Δσ, and this is only meant to
        # indicate how close the simulation is to the vertical stability limit, not to be exact.
        # (σ̇ = 0 at the surface, so the bottom layer is only seen through the layer above it.)
        courant = max(courant, w_max_cpu[k] / Δσ_cpu[k])
    end
    return Δt * courant
end

"""$(TYPEDSIGNATURES)
Progress line element that shows the global mean precipitation rate (rain and snow, large-scale and
convective) in mm/day. Left out if the model has no precipitation (`vars.parameterizations.rain_rate`)."""
struct Precipitation <: AbstractDiagnosticElement end

function print_diagnostic(::Precipitation, p::Progress, vars, model)
    (; parameterizations) = vars
    hasproperty(parameterizations, :rain_rate) || return ""
    (; area_weights) = model.geometry
    rate = dot(area_weights, parameterizations.rain_rate.data)       # m/s
    hasproperty(parameterizations, :snow_rate) && (rate += dot(area_weights, parameterizations.snow_rate.data))
    return @sprintf "precip %.2f mm/day" rate * 1000 * 86400
end

"""$(TYPEDSIGNATURES)
Progress line element that shows the global mean energy imbalance at the top of the atmosphere in W/m²:
incoming shortwave radiation minus outgoing shortwave minus outgoing longwave radiation, positive if the
planet gains energy. The global mean incoming shortwave radiation is `solar_constant / 4` (the planet
intercepts `πr² S₀` and has a surface of `4πr²`), the area-weighted mean of `solar_constant * cos_zenith`
on the grid differs from that by less than 0.1 W/m². Left out if the model has no radiation."""
struct EnergyImbalance <: AbstractDiagnosticElement end

function print_diagnostic(::EnergyImbalance, p::Progress, vars, model)
    (; parameterizations) = vars
    has_radiation = hasproperty(parameterizations, :outgoing_shortwave) &&
        hasproperty(parameterizations, :outgoing_longwave) && hasproperty(model, :planet)
    has_radiation || return ""
    (; area_weights) = model.geometry
    incoming = model.planet.solar_constant / 4
    outgoing = dot(area_weights, parameterizations.outgoing_shortwave.data) +
        dot(area_weights, parameterizations.outgoing_longwave.data)
    return @sprintf "TOA %+.1f W/m²" incoming - outgoing
end

"""$(TYPEDSIGNATURES)
Default elements of the progress lines of a `Feedback`, two lines: percentage, bar, ETA, simulation date
and speed, then maximum wind speed, temperature range, vertical Courant number, precipitation and energy
imbalance (those that apply to the model, none with Reactant, see [`OnReactant`](@ref))."""
default_elements() = (
    Percentage(), Bar(), ETA(), SimulationTime(), SimulationSpeed(), NewLine(indent = 5),
    MaximumWindSpeed(), TemperatureRange(), VerticalCourantNumber(), Precipitation(), EnergyImbalance(),
)

end # module ProgressElements
