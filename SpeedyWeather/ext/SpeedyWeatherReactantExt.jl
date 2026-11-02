module SpeedyWeatherReactantExt

using SpeedyWeather
using Reactant
using DocStringExtensions
using Dates

using SpeedyWeather: ReactantDevice, scale!, get_step, unpack, time_step!, time_steps!
using SpeedyWeather: reactant_compile_options

const ReactantDatesExt = Base.get_extension(
    Reactant, :ReactantDatesExt
)

const ReactantSimulation = Union{
    Simulation{V, <:BarotropicModel{SG, <:ReactantDevice}},
    Simulation{V, <:ShallowWaterModel{SG, <:ReactantDevice}}, Simulation{V, <:PrimitiveDryModel{SG, <:ReactantDevice}}, Simulation{V, <:PrimitiveWetModel{SG, <:ReactantDevice}},
} where {V, SG}

# TIME STEPPING WITH REACTANT
#
# IO can't happen inside a compiled program, so the time steps between two outputs are one program
# (a traced loop) and output and feedback happen between the programs on the host. Without output
# the whole run is one program (as Oceananigans' compiled `run!`).

"""
$(TYPEDSIGNATURES)
Time loop of a simulation with Reactant. Without output the whole run is one compiled program
([`time_steps!`](@ref), a traced loop). With output every `m` time steps, the `m` time steps between
two outputs are one program, called once per output interval; the spin-up steps of the time
stepper (which don't count for output) and a remainder of fewer than `m` steps at the end use
the compiled single time step. Output and feedback (progress meter, NaN check) happen between the
programs. Everything is compiled with `reactant_compile_options()` (kernels raised to StableHLO, some
patterns excluded)."""
function SpeedyWeather.time_stepping!(simulation::ReactantSimulation)
    (; variables, model) = simulation
    (; feedback, output) = model
    (; clock) = variables.prognostic

    n_steps = Int(clock.n_steps)                        # time stepper steps incl. spin-up steps
    m = output_every_n_steps(output)

    # without output: one program for the whole run
    if m == 0
        r_time_steps! = @compile compile_options = reactant_compile_options() time_steps!(simulation, n_steps)
        r_time_steps!(simulation, n_steps)
        SpeedyWeather.progress!(feedback, variables, model, n_steps)
        return simulation
    end

    # with output every m time steps: one program for the m steps between two outputs
    spin_up_steps = n_steps - Int(clock.n_time_steps)
    n_outputs, remainder = divrem(n_steps - spin_up_steps, m)
    r_time_steps! = @compile compile_options = reactant_compile_options() time_steps!(simulation, m)
    r_time_step! = spin_up_steps + remainder > 0 ? (@compile compile_options = reactant_compile_options() time_step!(simulation)) : nothing

    single_steps!(r_time_step!, simulation, spin_up_steps)
    for _ in 1:n_outputs
        r_time_steps!(simulation, m)
        SpeedyWeather.progress!(feedback, variables, model, m)
        SpeedyWeather.output!(output, simulation)
        nans_detected(feedback) && return simulation    # don't continue with NaNs
    end
    single_steps!(r_time_step!, simulation, remainder)
    return simulation
end

# n steps with the compiled single time step, each followed by feedback and output (only spin-up/remainder)
function single_steps!(r_time_step!, simulation, n)
    (; variables, model) = simulation
    for _ in 1:n
        r_time_step!(simulation)
        SpeedyWeather.progress!(model.feedback, variables, model, 1)
        SpeedyWeather.output!(model.output, simulation)
    end
    return nothing
end

SpeedyWeather.time_stepping!(simulation::ReactantSimulation, ::Nothing) = SpeedyWeather.time_stepping!(simulation)

"""
$(TYPEDSIGNATURES)
Time loop of a simulation with Reactant that calls a compiled single time step `r_time_step!` for every
time step, with output and feedback after every step. Example:

```julia
simulation = initialize!(model)
initialize!(simulation; steps=10) # don't forget this!
compile_options = SpeedyWeather.reactant_compile_options()
r_time_step! = @compile compile_options = compile_options SpeedyWeather.time_step!(simulation)
SpeedyWeather.time_stepping!(simulation, r_time_step!)
SpeedyWeather.finalize!(simulation)
```"""
function SpeedyWeather.time_stepping!(simulation::ReactantSimulation, r_time_step!)
    (; variables, model) = simulation
    (; clock) = variables.prognostic
    (; output, feedback) = model

    # n_steps (not n_time_steps) to match time_stepping!(::AbstractSimulation),
    # it includes time stepper spin-up steps (e.g. initial Euler step for Leapfrog)
    for _ in 1:Int(clock.n_steps)
        r_time_step!(simulation)
        # feedback and output are skipped inside the compiled time step
        SpeedyWeather.progress!(feedback, variables, model)
        SpeedyWeather.output!(output, simulation)
    end
    return simulation
end

"""
$(TYPEDSIGNATURES)
`n` time steps of `simulation` in one traced loop, compiled with `@compile`/`@jit` one program for
all `n` steps. Reverse-mode differentiation of the loop uses `checkpointing` (`false`,
`Reactant.Periodic(k)` or `Reactant.Binomial(budget)`) and `mincut` as Reactant's `@trace`."""
function SpeedyWeather.time_steps!(simulation::ReactantSimulation, n::Integer; checkpointing = false, mincut::Bool = true)
    # track_numbers = false: only trace Reactant numbers as @compile does, tracing the integers of our views fails
    @trace track_numbers = false mincut = mincut checkpointing = checkpointing for _ in 1:n
        time_step!(simulation)
    end
    return nothing
end

# output interval in time steps, 0 without (active) output
output_every_n_steps(::Nothing) = 0
output_every_n_steps(output::SpeedyWeather.AbstractOutput) = output.active ? output.output_every_n_steps : 0

nans_detected(::Nothing) = false
nans_detected(feedback::SpeedyWeather.Feedback) = feedback.nans_detected

SpeedyWeather.Clock(architecture::ReactantDevice) = Reactant.to_rarray(SpeedyWeather.Clock(), track_numbers = true)
SpeedyWeather.scalar_ref(::ReactantDevice, x::Number) = Ref(Reactant.ConcreteRNumber(x))

# TODO: move upstream to Reactant, ReactantDatesExt defines Dates.year(::ReactantDateTime) via
# Dates.year(days) for which Dates' method branches on the traced month, use yearmonthday instead
Dates.year(days::Reactant.TracedRNumber) = first(Dates.yearmonthday(days))

# TODO: move upstream to Reactant, Dates' method indexes a table with the (traced) month and branches
# with `&&`, used by the seasonal climatologies (`monthly_interpolation`)
function Dates.daysinmonth(year::Reactant.TracedRNumber, month::Reactant.TracedRNumber)
    thirty_days = (month == 4) | (month == 6) | (month == 9) | (month == 11)
    return ifelse(month == 2, 28 + Dates.isleapyear(year), ifelse(thirty_days, 30, 31))
end

# Whether the implicit solver has to be reinitialized (time step changed, or the reference temperature
# profile is scheduled to be recomputed) depends on the traced clock, so `update` is only known at run
# time. Instead of branching (`@trace if` would trace all captured model components) the reference
# temperature profile is only updated if `update`, and the operators, deterministic functions of that
# profile and the time step, are recomputed every time step. Equivalent to `reinitialize!` with a
# `update::Bool` for other architectures.
# TODO: This is certainly a bit wasteful, however it matters less than one could think. It adds
# roughly 0.5ms/step. That makes a T32L8 run 10% slower, but already at T64L8 it's basically within
# measuring noise. We still might revise this in the future to either a) explicilty trace whether
# the time step changed or b) just handle the first few steps sepereatly for Reactant (a bit similar
# to what we used to do in general for the time stepping)
function SpeedyWeather.reinitialize!(
        implicit::ImplicitPrimitiveEquation,
        update::Reactant.TracedRNumber{Bool},
        Δt,
        model::SpeedyWeather.PrimitiveEquation,
        vars::SpeedyWeather.Variables,
    )
    (; geometry, geopotential, atmosphere, adiabatic_conversion) = model
    (; temp_profile) = implicit
    Tₖ = vars.dynamics.average_temperature_profile          # current global-mean temperature profile
    temp_profile .= ifelse.(update, Tₖ, temp_profile)

    scale = vars.prognostic.scale[]         # implicit solver needs to be initialized with scaled time step
    initialize!(implicit, Δt / scale, temp_profile, geometry, geopotential, atmosphere, adiabatic_conversion)
    return nothing
end

# The executions of a schedule are counted in a plain integer of the model, a traced decision can't
# be added to it, so executions aren't counted with Reactant
SpeedyWeather.count_execution!(::SpeedyWeather.AbstractSchedule, ::Reactant.TracedRNumber{Bool}) = nothing

SpeedyWeather.dilate(Δt::Millisecond, dilation::Reactant.TracedRNumber) =
    ReactantDatesExt.ReactantMillisecond(round(Int, Δt.value * dilation))

SpeedyWeather.dilate(Δt::ReactantDatesExt.ReactantMillisecond, dilation) =
    ReactantDatesExt.ReactantMillisecond(round(Int, Δt.value * dilation))

# OUTPUT HANDLING FOR REACTANT
#
# Reactant cannot trace through the output! pipeline: `output!` interpolates onto
# CPU scratch fields and writes to disk, both of which are runtime side effects
# the tracer can't (and shouldn't) capture. We split the behaviour by tracing context:
#   - Inside `@compile` (`within_compile()` == true): `do_output!`, which every output
#     writer calls to decide whether to write, returns false. Nothing gets baked into
#     the compiled `time_step!`, so the compiled function only does numerics.
#   - Outside `@compile`: `time_stepping!` (above) explicitly invokes
#     `output!(output, simulation)` after each compiled step. The simulation's
#     data is then ConcretePJRTArrays, which the standard output path can
#     consume via `on_architecture(CPU(), …)`.
const ReactantClock = SpeedyWeather.Clock{<:ReactantDatesExt.ReactantDateTime}

# same for the feedback (progress meter, NaN detection), skipped inside `@compile`
const ReactantModel = Union{
    BarotropicModel{SG, <:ReactantDevice}, ShallowWaterModel{SG, <:ReactantDevice},
    PrimitiveDryModel{SG, <:ReactantDevice}, PrimitiveWetModel{SG, <:ReactantDevice},
} where {SG}

function SpeedyWeather.progress!(feedback::SpeedyWeather.Feedback, vars::SpeedyWeather.Variables, model::ReactantModel)
    Reactant.ReactantCore.within_compile() && return nothing
    return @invoke SpeedyWeather.progress!(feedback::SpeedyWeather.Feedback, vars::SpeedyWeather.Variables, model::SpeedyWeather.AbstractModel)
end

# Don't trace into the feedback (it is not part of the computation and skipped inside `@compile`),
# its progress meter holds IO buffers that Reactant cannot trace, return it unchanged instead,
# identified only by its type for the compile cache
Reactant.traced_type_inner(
    @nospecialize(T::Type{<:SpeedyWeather.AbstractFeedback}), seen, @nospecialize(mode::Reactant.TraceMode),
    @nospecialize(track_numbers::Type), @nospecialize(ndevices), @nospecialize(runtime),
) = T

function Reactant.make_tracer(
        seen, @nospecialize(prev::SpeedyWeather.AbstractFeedback), @nospecialize(path), mode; kwargs...
    )
    if mode == Reactant.TracedToTypes
        push!(path, Core.Typeof(prev))
        return nothing
    end
    return prev
end

function SpeedyWeather.do_output!(
        core::SpeedyWeather.OutputWriterCore,
        clock::ReactantClock,
        output::SpeedyWeather.AbstractOutput,
    )
    Reactant.ReactantCore.within_compile() && return false
    # check on the host first: reading the clock from the device waits for the time step to finish,
    # so the next time step can't be queued meanwhile (4x slower with a sync every time step)
    output.active || return false
    time_step_counter = Int(clock.time_step_counter)
    do_output = time_step_counter > 0 && time_step_counter % core.output_every_n_steps == 0
    core.output_counter += do_output
    return do_output
end

# materialise into a plain DateTime and forward to the output-specific DateTime method (NetCDF, Zarr, …)
SpeedyWeather.output!(output::SpeedyWeather.AbstractOutput, time::ReactantDatesExt.ReactantDateTime) =
    SpeedyWeather.output!(output, convert(DateTime, time))

end
