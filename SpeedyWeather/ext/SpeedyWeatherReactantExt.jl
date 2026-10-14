module SpeedyWeatherReactantExt

using SpeedyWeather
using Reactant
using DocStringExtensions
using Dates

using SpeedyWeather: ReactantDevice, scale!, get_step, unpack, time_step!

const ReactantDatesExt = Base.get_extension(
    Reactant, :ReactantDatesExt
)

const ReactantSimulation = Union{
    Simulation{V, <:BarotropicModel{SG, <:ReactantDevice}},
    Simulation{V, <:ShallowWaterModel{SG, <:ReactantDevice}}, Simulation{V, <:PrimitiveDryModel{SG, <:ReactantDevice}}, Simulation{V, <:PrimitiveWetModel{SG, <:ReactantDevice}},
} where {V, SG}

# time stepping functions with Reactant, take compiled functions as optional arguments
# in case they are not provided, they are compiled on the fly

"""
$(TYPEDSIGNATURES)

Time-stepping function for a simulation with Reactant, take compiled functions as optional arguments.
In case they are not provided, they are compiled on the fly.

Example usage:

```julia
simulation = initialize!(model)
initialize!(simulation; steps=10) # don't forget this!
r_time_step! = @compile raise = true SpeedyWeather.time_step!(simulation)
SpeedyWeather.time_stepping!(simulation, r_time_step!)
SpeedyWeather.finalize!(simulation)
```

The time step is compiled with `raise = true` by default, i.e. KernelAbstractions kernels are raised
to StableHLO: kernels that aren't raised can be handed operands with a non-default layout on CUDA
and read scrambled data (EnzymeAD/Reactant.jl#3377).
"""
function SpeedyWeather.time_stepping!(simulation::ReactantSimulation, r_time_step! = nothing, enable_checkpointing = true)
    if isnothing(r_time_step!)
        @info "Reactant compiling time_step!"
        r_time_step! = @compile raise = true time_step!(simulation)
    end

    (; variables, model) = simulation
    (; clock) = variables.prognostic
    (; output, feedback) = model

    #TODO: reenable @trace once Reactant issues fixed
    #@trace checkpointing = enable_checkpointing for _ in 1:clock.n_steps
    #    r_time_step!(simulation)
    #end

    # n_steps (not n_time_steps) to match time_stepping!(::AbstractSimulation),
    # it includes time stepper spin-up steps (e.g. initial Euler step for Leapfrog)
    for _ in 1:Int(clock.n_steps)
        r_time_step!(simulation)
        # feedback and output are skipped inside the compiled time step
        SpeedyWeather.progress!(feedback, variables, model)
        SpeedyWeather.output!(output, simulation)
    end
    return
end

SpeedyWeather.Clock(architecture::ReactantDevice) = Reactant.to_rarray(SpeedyWeather.Clock(), track_numbers = true)
SpeedyWeather.scalar_ref(::ReactantDevice, x::Number) = Ref(Reactant.ConcreteRNumber(x))

# TODO: move upstream to Reactant, ReactantDatesExt defines Dates.year(::ReactantDateTime) via
# Dates.year(days) for which Dates' method branches on the traced month, use yearmonthday instead
Dates.year(days::Reactant.TracedRNumber) = first(Dates.yearmonthday(days))

const ReactantPrimitiveEquation = Union{
    PrimitiveDryModel{SG, <:ReactantDevice}, PrimitiveWetModel{SG, <:ReactantDevice},
} where {SG}

# The time step depends on the (traced) clock, so whether the implicit solver has to be
# reinitialized (time step changed) is only known at run time. Instead of branching (`@trace if`
# would trace all captured model components) the temperature profile the operators are built
# from is only updated if the time step changed, and the operators, deterministic functions of
# that profile and the time step, are recomputed every time step. Equivalent to reinitializing
# only when the time step changed as `reinitialize!` does for other architectures.
# TODO: This is certainly a bit wasteful, however it matters less than one could think. It adds 
# roughly 0.5ms/step. That makes a T32L8 run 10% slower, but already at T64L8 it's basically within 
# measuring noise. We still might revise this in the future to either a) explicilty trace whether
# the time step changed or b) just handle the first few steps sepereatly for Reactant (a bit similar
# to what we used to do in general for the time stepping)
function SpeedyWeather.reinitialize!(
        implicit::ImplicitPrimitiveEquation,
        model::ReactantPrimitiveEquation,
        vars::SpeedyWeather.Variables,
    )
    (; time_stepping, geometry, geopotential, atmosphere, adiabatic_conversion) = model
    Δt = SpeedyWeather.time_step(time_stepping, vars.prognostic.clock)
    Δt_changed = implicit.Δt[] != Δt
    implicit.Δt[] = Δt

    (; temp_profile) = implicit
    Tₖ = vars.dynamics.average_temperature_profile
    temp_profile .= ifelse.(Δt_changed, Tₖ, temp_profile)

    scale = vars.prognostic.scale[]         # implicit solver needs to be initialized with scaled time step
    initialize!(implicit, Δt / scale, temp_profile, geometry, geopotential, atmosphere, adiabatic_conversion)
    return nothing
end

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
    time_step_counter = Int(clock.time_step_counter)
    do_output = output.active && time_step_counter > 0 && time_step_counter % core.output_every_n_steps == 0
    core.output_counter += do_output
    return do_output
end

# materialise into a plain DateTime and forward to the output-specific DateTime method (NetCDF, Zarr, …)
SpeedyWeather.output!(output::SpeedyWeather.AbstractOutput, time::ReactantDatesExt.ReactantDateTime) =
    SpeedyWeather.output!(output, convert(DateTime, time))

end
