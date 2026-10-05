abstract type AbstractFeedback <: AbstractModelComponent end

export Feedback, VerticalCourantNumber, SimulationTime, SimulationSpeed, MaximumWindSpeed, TemperatureRange, default_elements

"""
Feedback struct that contains options and object for command-line feedback
like the progress meter.
$(TYPEDFIELDS)"""
@kwdef mutable struct Feedback <: AbstractFeedback
    "[OPTION] print feedback to REPL?, default is isinteractive(), true in interactive REPL mode"
    verbose::Bool = isinteractive()

    "[OPTION] check for NaNs in the prognostic variables"
    debug::Bool = true

    "[OPTION] Progress description"
    description::String = ""

    "[OPTION] Progress bar length in characters, `nothing` fits it to the terminal width, `0` shows no bar"
    progress_bar_length::Union{Int, Nothing} = 0

    "[OPTION] show speed (e.g. in simulated years per day) in progress meter?"
    showspeed::Bool = true

    "[OPTION] Minimum wallclock time between feedback updates"
    feedback_dt::Float32 = 0.1

    "[OPTION] Show simulation time?"
    show_time::Bool = true

    "[OPTION] Show maximum speed of the simulated flow [m/s]"
    show_umax::Bool = true

    "[OPTION] Show temperature range of simulation [˚C]"
    show_temperature_range::Bool = true

    "[OPTION] Interval in time steps between NaN checks, and between the progress meter diagnostics (maximum speed, temperature range)"
    interval::Int = 50

    """[OPTION] Elements of the progress line (subtypes of `ProgressMeter.AbstractProgressElement`),
    `nothing` builds them from `showspeed`, `show_time`, `show_umax`, `show_temperature_range`,
    see `default_elements`. Elements specific to SpeedyWeather (`SimulationTime`, `SimulationSpeed`,
    `MaximumWindSpeed`, `TemperatureRange`, `VerticalCourantNumber`) are bound to the simulation
    in `initialize!`, e.g. `Feedback(elements = (default_elements(Feedback())..., VerticalCourantNumber()))`"""
    elements::Union{Nothing, Tuple} = nothing

    "[DERIVED] struct containing everything progress related"
    progress_meter::ProgressMeter.Progress =
        ProgressMeter.Progress(1, enabled = verbose)

    "[DERIVED] did NaNs occur in the simulation?"
    nans_detected::Bool = false
end

function Base.show(io::IO, P::ProgressMeter.Progress)
    println(io, "$(typeof(P)) <: ProgressMeter.AbstractProgress")
    keys = propertynames(P)
    return print_fields(io, P, keys)
end

"""
$(TYPEDSIGNATURES)
Initializes the a `Feedback` struct."""
function initialize!(feedback::Feedback, variables::Variables, model::AbstractModel)
    (; clock) = variables.prognostic

    # set to false to recheck for NaNs
    feedback.nans_detected = false

    # reinitalize progress meter, minus one to exclude first_timesteps! which contain compilation
    # only do now for benchmark accuracy
    (; description, verbose, feedback_dt) = feedback
    elements = something(feedback.elements, default_elements(feedback))
    elements = map(element -> bind_element(element, variables, model), elements)
    desc = description * (model.output.active ? " $(model.output.run_folder) " : " ")
    feedback.progress_meter = ProgressMeter.Progress(
        clock.n_steps;          # use time stepper steps (regardless Δt) not time steps of size Δt
        enabled = verbose,
        elements,
        desc,
        color = :blue,
        barlen = feedback.progress_bar_length,
        barglyphs = ProgressMeter.BarGlyphs(" ━━  "),
        dt = feedback_dt,
    )

    return nothing
end

progress!(feedback::Feedback) = ProgressMeter.next!(feedback.progress_meter)

function progress!(feedback::Feedback, vars::Variables, model::AbstractModel)
    (; counter, n) = feedback.progress_meter
    interval = max(1, feedback.interval)

    progress!(feedback)

    last_step = counter == n - 1
    feedback.debug && (mod(counter, interval) == 0 || last_step) &&
        nan_detection!(feedback, vars, model)
    return nothing
end

# fallback for feedback = nothing
progress!(::Nothing, vars::Variables, model::AbstractModel) = nothing

"""
$(TYPEDSIGNATURES)
Finalises the progress meter and the progress txt file."""
finalize!(F::Feedback) = ProgressMeter.finish!(F.progress_meter)

# fallback if feedback is set to nothing
finalize!(::Nothing) = nothing

"""$(TYPEDSIGNATURES)
Detect NaN (Not-a-Number, or Inf) in the prognostic variables."""
function nan_detection!(feedback::Feedback, vars::Variables, model::AbstractModel)
    feedback.nans_detected && return nothing        # escape immediately if nans already detected
    i = feedback.progress_meter.counter             # time step
    vor = get_prognostic_step(vars.prognostic.vorticity, model.time_stepping, feedback)
    GPUArrays.@allowscalar vor0 = vor[2, end]       # only check 1-0 mode of surface vorticity

    # just check first harmonic, spectral transform propagates NaNs globally anyway
    (; time) = vars.prognostic.clock                                # current time for feedback
    nans_detected_here = ~isfinite(vor0)
    nans_detected_here && @warn "NaN or Inf detected at time step $i ($time)"
    return feedback.nans_detected = nans_detected_here
end

# Elements of the progress line. They are created unbound, e.g. `VerticalCourantNumber()`, and
# `bind_element` returns a copy that holds what it needs from the simulation (`Variables`, time
# step, ...) in `initialize!(::Feedback, ...)`. `print_element` is only called when the progress
# meter is redrawn, so diagnostics are only computed when they are displayed.

# elements that do not need anything from the simulation are their own bound version
bind_element(element, vars, model) = element

"""$(TYPEDSIGNATURES)
Progress line element that shows the current simulation date."""
struct SimulationTime{C} <: ProgressMeter.AbstractProgressElement
    clock::C
end
SimulationTime() = SimulationTime(nothing)
bind_element(::SimulationTime, vars, model) = SimulationTime(vars.prognostic.clock)
ProgressMeter.print_element(element::SimulationTime, p, status) = string(Dates.Date(element.clock.time))
ProgressMeter.print_element(::SimulationTime{Nothing}, p, status) = ""

"""$(TYPEDSIGNATURES)
Progress line element that shows the simulation speed, e.g. in simulated years per day,
which needs the time step `Δt` [s]. `separator` is printed in front."""
@kwdef struct SimulationSpeed <: ProgressMeter.AbstractProgressElement
    separator::String = ", "
    Δt::Float64 = 0.0
end
bind_element(element::SimulationSpeed, vars, model) = SimulationSpeed(element.separator, Float64(model.time_stepping.Δt))

function ProgressMeter.print_element(element::SimulationSpeed, p, status)
    sec_per_iter = status.elapsed / max(1, p.counter - p.start)
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
Progress line element that shows the maximum wind speed [m/s] of the current grid-space `u`, `v`.
Shows nothing if the model has no `u`."""
struct MaximumWindSpeed{V} <: ProgressMeter.AbstractProgressElement
    vars::V
end
MaximumWindSpeed() = MaximumWindSpeed(nothing)
bind_element(::MaximumWindSpeed, vars, model) = MaximumWindSpeed(vars)
ProgressMeter.print_element(::MaximumWindSpeed{Nothing}, p, status) = ""

function ProgressMeter.print_element(element::MaximumWindSpeed, p, status)
    hasproperty(element.vars.grid, :u) || return ""
    umin, umax = extrema(element.vars.grid.u)
    return @sprintf ", %3d m/s" max(abs(umin), abs(umax))
end

"""$(TYPEDSIGNATURES)
Progress line element that shows the range of the grid-space temperature in ˚C.
Shows nothing if the model has no `temperature`."""
struct TemperatureRange{V} <: ProgressMeter.AbstractProgressElement
    vars::V
end
TemperatureRange() = TemperatureRange(nothing)
bind_element(::TemperatureRange, vars, model) = TemperatureRange(vars)
ProgressMeter.print_element(::TemperatureRange{Nothing}, p, status) = ""

function ProgressMeter.print_element(element::TemperatureRange, p, status)
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
Feedback(elements = (default_elements(Feedback())..., VerticalCourantNumber()))
```"""
struct VerticalCourantNumber{V, T} <: ProgressMeter.AbstractProgressElement
    vars::V
    Δσ::T
    Δt::Float64
end
VerticalCourantNumber() = VerticalCourantNumber(nothing, nothing, 0.0)

function bind_element(::VerticalCourantNumber, vars, model)
    hasproperty(model, :geometry) || return VerticalCourantNumber()
    return VerticalCourantNumber(vars, model.geometry.σ_levels_thick, Float64(model.time_stepping.Δt))
end
ProgressMeter.print_element(::VerticalCourantNumber{Nothing}, p, status) = ""

function ProgressMeter.print_element(element::VerticalCourantNumber, p, status)
    hasproperty(element.vars.dynamics, :w) || return ""
    scale = element.vars.prognostic.scale[]     # divergence, hence w, is scaled by the radius in the dynamical core
    return @sprintf ", Cᵥ = %.2f" vertical_courant_number(element.vars.dynamics.w, element.Δσ, element.Δt / scale)
end

"""$(TYPEDSIGNATURES)
Maximum vertical Courant number of layer `k`, `max(|σ̇ₖ₊₁/₂|, |σ̇ₖ₋₁/₂|) Δt / Δσₖ`, over all layers
and grid points. `w` is the vertical velocity at the layer interfaces `k+1/2` (zero at the surface),
`Δσ` the layer thickness. Pass `Δt / scale` if `w` is radius-scaled as in the dynamical core."""
function vertical_courant_number(w, Δσ, Δt)
    w_data = w.data
    nlayers = size(w_data, 2)
    courant_below = maximum(abs.(w_data) ./ reshape(Δσ, 1, :))
    courant_above = nlayers > 1 ? maximum(abs.(view(w_data, :, 1:(nlayers - 1))) ./ reshape(view(Δσ, 2:nlayers), 1, :)) : zero(courant_below)
    return Δt * max(courant_below, courant_above)
end

"""$(TYPEDSIGNATURES)
Default elements of the progress line of a `Feedback`, in the order description, percentage, bar, ETA,
and in parenthesis (if `showspeed`) simulation date, speed, maximum wind speed and temperature range,
each depending on the options `showspeed`, `show_time`, `show_umax`, `show_temperature_range`."""
function default_elements(feedback::Feedback)
    elements = Any[ProgressMeter.Description(), ProgressMeter.Percentage(), ProgressMeter.Bar(), ProgressMeter.ETA()]
    if feedback.showspeed
        push!(elements, " (")
        feedback.show_time && push!(elements, SimulationTime())
        push!(elements, SimulationSpeed(separator = feedback.show_time ? ", " : ""))
        feedback.show_umax && push!(elements, MaximumWindSpeed())
        feedback.show_temperature_range && push!(elements, TemperatureRange())
        push!(elements, ")")
    end
    return Tuple(elements)
end

export ParametersTxt

"""ParametersTxt callback. Writes a parameters.txt file with all model parameters.
Options are $(TYPEDFIELDS)"""
@kwdef mutable struct ParametersTxt <: AbstractCallback
    "[OPTION] Path for parameters.txt file, uses model.output.run_path if not specified"
    path::String = ""

    "[OPTION] File name for parameters.txt file"
    filename::String = "parameters.txt"

    "[OPTION] Only write with model.output.active = true?"
    write_only_with_output::Bool = true
end

"""$(TYPEDSIGNATURES)
Initialize ParametersTxt by writing the model parameters (via defined show of model components) into a text file."""
function initialize!(parameters_txt::ParametersTxt, vars, model)

    # escape in case of no output
    parameters_txt.write_only_with_output && (model.output.active || return nothing)

    (; filename) = parameters_txt
    path = parameters_txt.path == "" ? model.output.run_path : parameters_txt.path
    mkpath(path)

    # also export parameters into run????/parameters.txt
    file = open(joinpath(path, filename), "w")
    for property in propertynames(model)
        println(file, "model.$property")
        println(file, getfield(model, property), "\n")
    end
    close(file)

    model.output.active || @info "Parameter summary written to $(joinpath(path, filename)) although output=false"
    return nothing
end

callback!(::ParametersTxt, args...) = nothing
finalize!(::ParametersTxt, arg...) = nothing


export ProgressTxt

"""ProgressTxt callback. Writes a progress.txt file with time stepping progress.
Options are $(TYPEDFIELDS)"""
@kwdef mutable struct ProgressTxt <: AbstractCallback
    "[OPTION] Path for progress.txt file, uses model.output.run_path if not specified"
    path::String = ""

    "[OPTION] File name for progress.txt file"
    filename::String = "progress.txt"

    "[OPTION] Only write with model.output.active = true?"
    write_only_with_output::Bool = true

    "[OPTION] Every n% of time steps write to progress.txt, default is 5%"
    every_n_percent::Int = 5

    "[DERIVED] IOStream for progress.txt file"
    file::IOStream = IOStream("")
end

"""$(TYPEDSIGNATURES)
Initializes the ProgressTxt callback by creating a progress.txt file and writing some initial information to it."""
function initialize!(progress_txt::ProgressTxt, vars, model)
    # escape in case of no output
    progress_txt.write_only_with_output && (model.output.active || return nothing)

    (; filename) = progress_txt
    path = progress_txt.path == "" ? model.output.run_path : progress_txt.path
    mkpath(path)

    (; run_folder, run_path) = model.output
    SG = model.spectral_grid
    L = model.time_stepping
    days = Second(vars.prognostic.clock.period).value / (3600 * 24)

    # create progress.txt file in run_????/
    file = open(joinpath(path, filename), "w")
    s = "Starting SpeedyWeather.jl $run_folder on " *
        Dates.format(Dates.now(), Dates.RFC1123Format)
    write(file, s * "\n")
    write(file, "Integrating:\n")
    write(file, "$SG\n")
    write(file, "Time: $days days at Δt = $(L.Δt)s\n")
    model.output.active && write(file, "\nAll data will be stored in $run_path\n")
    model.output.active || write(file, "\nNo output will be written (output=false)\n")
    progress_txt.file = file

    model.output.active || @info "Progress is being written to $(joinpath(path, filename)) although output=false"
    return nothing
end

"""$(TYPEDSIGNATURES)
Writes the time stepping progress to the progress.txt file every `every_n_percent` % of time steps."""
function callback!(progress_txt::ProgressTxt, vars, model)
    # escape in case of no output
    progress_txt.write_only_with_output && (model.output.active || return nothing)
    isnothing(model.feedback) && return nothing

    (; progress_meter, nans_detected) = model.feedback
    (; counter, n) = progress_meter
    (; file, every_n_percent) = progress_txt

    # occasionally write progress to txt file
    if (counter / n * 100 % 1) > ((counter + 1) / n * 100 % 1)
        percent = round(Int, (counter + 1) / n * 100)       # % of time steps completed
        if (percent % every_n_percent == 0)                 # write every p% step in txt
            write(file, @sprintf("\n%3d%%", percent))
            r = remaining_time(progress_meter)
            write(file, ", ETA: $r")

            time_elapsed = progress_meter.tlast - progress_meter.tinit
            s = speedstring(time_elapsed / counter, model.time_stepping.Δt)
            write(file, ", $s")

            nans_detected && write(file, ", NaN/Inf detected.")
            flush(file)
        end
    end
    return nothing
end

"""$(TYPEDSIGNATURES)
Finalizes the ProgressTxt callback by writing the total time taken to the progress.txt file and closing it."""
function finalize!(progress_txt::ProgressTxt, vars, model)
    # escape in case of no output
    progress_txt.write_only_with_output && (model.output.active || return nothing)
    isnothing(model.feedback) && return nothing

    (; file) = progress_txt
    (; progress_meter) = model.feedback

    time_elapsed = progress_meter.tlast - progress_meter.tinit
    s = "\nIntegration done in $(readable_secs(time_elapsed))."
    write(file, "\n$s\n")
    flush(file)
    close(file)
    return nothing
end

"""$(TYPEDSIGNATURES)
Estimates the remaining time from a `ProgresssMeter.Progress`. Adapted from ProgressMeter.jl"""
function remaining_time(p::ProgressMeter.Progress)
    elapsed_time = time() - p.tinit
    est_total_time = elapsed_time * (p.n - p.start) / (p.counter - p.start)
    if 0 <= est_total_time <= typemax(Int)
        eta_sec = round(Int, est_total_time - elapsed_time)
        eta = ProgressMeter.durationstring(eta_sec)
    else
        eta = "N/A"
    end
    return eta
end
