# Reactant-specific methods of the benchmark hooks in `benchmark_suite.jl`. Included by the driver
# only for the `reactant-*` architectures, after `using Reactant`, as the `@compile` macro must be
# resolvable when this file is parsed.

# Without output and feedback (`build_model`), `run!` runs the whole time loop as one compiled
# program (`time_stepping!` in SpeedyWeatherReactantExt) but compiles it on every call, timing it
# would therefore include the compilation. Instead compile the time loop of the timed `period`
# (`time_steps!` for all its steps, one program) once here and use it for both the warm-up and the
# timed run. The compiled program takes the simulation as argument so it is reused for a
# re-initialized simulation of the same model. The compile time is reported separately, as for
# Reactant it is a headline number in itself. Compiled with `SpeedyWeather.reactant_compile_options()` (kernels raised, some patterns excluded) as in `time_stepping!`.
function build_runner(::ReactantDevice, simulation, period)
    (; model) = simulation
    clock = SpeedyWeather.Clock()       # number of steps of the timed run, incl. spin-up steps
    initialize!(clock, model.time_stepping, model.planet, period)
    n_steps = Int(clock.n_steps)
    compile_time = @elapsed compiled_time_steps! = @compile compile_options = SpeedyWeather.reactant_compile_options() SpeedyWeather.time_steps!(simulation, n_steps)
    @info @sprintf("     compiled the time loop (%d steps) in %.1f s", n_steps, compile_time)
    time_stepping(simulation) = compiled_time_steps!(simulation, n_steps)
    return (; time_stepping, compile_time)
end

# `synchronize(::ReactantDevice)` doesn't wait (KernelAbstractions.synchronize is a no-op for
# Reactant's backend). Reading the clock, an output of every compiled time step, from the device
# waits for the last time step to have finished, which depends on all previous ones.
function synchronize_device(::ReactantDevice, simulation)
    Int(simulation.variables.prognostic.clock.time_step_counter)
    return nothing
end

# `Base.summarysize` counts the wrapper of a Reactant array but not its buffer on the device,
# add the buffers' sizes. Mutable objects (incl. the arrays) are counted once.
memory_footprint(::ReactantDevice, simulation) = Base.summarysize(simulation) + device_bytes(simulation)

function device_bytes(x, seen = Base.IdSet{Any}())
    x isa Union{Module, Function, Type, Symbol, String} && return 0
    isbits(x) && return 0
    if ismutable(x)
        x in seen && return 0
        push!(seen, x)
    end
    x isa Reactant.AbstractConcreteArray && return length(x) * sizeof(eltype(x))
    if x isa DenseArray  # host arrays, wrappers (views, LowerTriangularArray, Field) are recursed into
        isbitstype(eltype(x)) && return 0
        return sum((device_bytes(x[i], seen) for i in eachindex(x) if isassigned(x, i)); init = 0)
    end
    return sum((device_bytes(getfield(x, i), seen) for i in 1:nfields(x) if isdefined(x, i)); init = 0)
end

# Compile `f(args...)` once, `sync = true` so that every call waits for the device to finish
function compiled_call(::ReactantDevice, f, args...)
    compiled_f = @compile compile_options = SpeedyWeather.reactant_compile_options(; sync = true) f(args...)
    return () -> compiled_f(args...)
end
