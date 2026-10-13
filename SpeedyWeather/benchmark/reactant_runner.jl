# Reactant-specific methods of the benchmark hooks in `benchmark_suite.jl`. Included by the driver
# only for the `reactant-*` architectures, after `using Reactant`, as the `@compile` macro must be
# resolvable when this file is parsed.

# `run!` compiles the time step on every call (see `time_stepping!` in SpeedyWeatherReactantExt),
# timing it would therefore include the compilation. Instead compile the time step once here and
# pass it to `time_stepping!` for both the warm-up and the timed run. The compiled time step takes
# the simulation as argument so it is reused for a re-initialized simulation of the same model.
# The compile time is reported separately, as for Reactant it is a headline number in itself.
function build_runner(::ReactantDevice, simulation)
    compile_time = @elapsed compiled_time_step! = @compile SpeedyWeather.time_step!(simulation)
    @info @sprintf("     compiled the time step in %.1f s", compile_time)
    time_stepping(simulation) = SpeedyWeather.time_stepping!(simulation, compiled_time_step!)
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
    compiled_f = @compile sync = true f(args...)
    return () -> compiled_f(args...)
end
