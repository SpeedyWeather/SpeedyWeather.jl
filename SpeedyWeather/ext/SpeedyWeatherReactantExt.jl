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
r_time_step! = @compile SpeedyWeather.time_step!(simulation)
SpeedyWeather.time_stepping!(simulation, r_time_step!)
SpeedyWeather.finalize!(simulation)
```
"""
function SpeedyWeather.time_stepping!(simulation::ReactantSimulation, r_time_step! = nothing, enable_checkpointing = true)
    if isnothing(r_time_step!)
        @info "Reactant compiling time_step!"
        r_time_step! = @compile time_step!(simulation)
    end

    clock = simulation.variables.prognostic.clock

    #TODO: reenable @trace once Reactant issues fixed
    #@trace checkpointing = enable_checkpointing for _ in 1:clock.n_steps
    #    r_time_step!(simulation)
    #end

    # n_steps (not n_time_steps) to match time_stepping!(::AbstractSimulation),
    # it includes time stepper spin-up steps (e.g. initial Euler step for Leapfrog)
    for _ in 1:Int(clock.n_steps)
        r_time_step!(simulation)
    end
    return
end

# that's for Reactant TracableDateTime
SpeedyWeather.secondofday(dt::ReactantDatesExt.ReactantDateTime) = Dates.second(ReactantDatesExt.ReactantTime(dt).instant)

SpeedyWeather.Clock(architecture::ReactantDevice) = Reactant.to_rarray(SpeedyWeather.Clock(), track_numbers = true)

end
