module SpeedyWeatherReactantExt

using SpeedyWeather
using Reactant
using DocStringExtensions
using Dates

using SpeedyWeather: ReactantDevice, scale!, get_step, unpack, timestep!, first_timesteps!, later_timestep!

const ReactantDatesExt = Base.get_extension(
    Reactant, :ReactantDatesExt
)

const ReactantSimulation = Union{
    Simulation{<:BarotropicModel{SG, <:ReactantDevice}},
    Simulation{<:ShallowWaterModel{SG, <:ReactantDevice}}, Simulation{<:PrimitiveDryModel{SG, <:ReactantDevice}}, Simulation{<:PrimitiveWetModel{SG, <:ReactantDevice}},
} where {SG}

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
r_first! = @compile SpeedyWeather.first_timesteps!(simulation)
r_later! = @compile SpeedyWeather.later_timestep!(simulation)
SpeedyWeather.time_stepping!(simulation, r_first!, r_later!)
SpeedyWeather.finalize!(simulation)
```
"""
function SpeedyWeather.time_stepping!(simulation::ReactantSimulation, r_first_timesteps! = nothing, r_later_timestep! = nothing, enable_checkpointing = true)
    if isnothing(r_first_timesteps!)
        @info "Reactant compiling first_timesteps!"
        r_first_timesteps! = @compile first_timesteps!(simulation)
    end

    if isnothing(r_later_timestep!)
        @info "Reactant compiling later_timestep!"
        r_later_timestep! = @compile later_timestep!(simulation)
    end

    (; clock) = simulation.prognostic_variables

    #TODO: reenable @trace once Reactant issues fixed
    #@trace checkpointing = enable_checkpointing for _ in clock.timestep_counter:clock.n_timesteps
    #    r_later_timestep!(simulation)
    #end

    r_first_timesteps!(simulation)

    for _ in (Int(clock.timestep_counter) + 1):Int(clock.n_timesteps)
        r_later_timestep!(simulation)
    end
    return
end

# that's for Reactant TracableDateTime
SpeedyWeather.secondofday(dt::ReactantDatesExt.ReactantDateTime) = Dates.second(ReactantDatesExt.ReactantTime(dt).instant)

#TODO: move those the the ReactantDatesExt once I am sure it's all additional functionality I need to add
#Dates.firstdayofmonth(dt::ReactantDatesExt.ReactantDateTime) = ReactantDatesExt.ReactantDate(dt).date
Dates.firstdayofmonth(dt::ReactantDatesExt.ReactantDate) = ReactantDatesExt.ReactantDate(Dates.UTD(Dates.value(dt) - Dates.day(dt) + 1))
Dates.firstdayofmonth(dt::ReactantDatesExt.ReactantDateTime) = ReactantDatesExt.ReactantDateTime(Dates.firstdayofmonth(convert(ReactantDatesExt.ReactantDate, dt)))

# ReactantDatesExt versions of Base.convert methods (mirrors Dates.jl conversions)
Base.convert(::Type{ReactantDatesExt.ReactantDateTime}, dt::ReactantDatesExt.ReactantDate) =
    ReactantDatesExt.ReactantDateTime(Dates.UTInstant(ReactantDatesExt.ReactantMillisecond(Dates.value(dt) * 86400000)))
Base.convert(::Type{ReactantDatesExt.ReactantDate}, dt::ReactantDatesExt.ReactantDateTime) =
    ReactantDatesExt.ReactantDate(Dates.UTInstant(ReactantDatesExt.ReactantDay(div(Dates.value(dt), 86400000))))
Base.convert(::Type{ReactantDatesExt.ReactantTime}, dt::ReactantDatesExt.ReactantDateTime) =
    ReactantDatesExt.ReactantTime(ReactantDatesExt.ReactantNanosecond((Dates.value(dt) % 86400000) * 1000000))

Base.convert(::Type{ReactantDatesExt.ReactantDateTime}, x::Millisecond) =
    ReactantDatesExt.ReactantDateTime(Dates.UTInstant(ReactantDatesExt.ReactantMillisecond(Dates.value(x))))
Base.convert(::Type{Millisecond}, dt::ReactantDatesExt.ReactantDateTime) =
    Millisecond(Dates.value(dt))
Base.convert(::Type{ReactantDatesExt.ReactantDate}, x::Day) =
    ReactantDatesExt.ReactantDate(Dates.UTInstant(ReactantDatesExt.ReactantDay(Dates.value(x))))
Base.convert(::Type{Day}, dt::ReactantDatesExt.ReactantDate) =
    Day(Dates.value(dt))

SpeedyWeather.Clock(architecture::ReactantDevice) = Reactant.to_rarray(SpeedyWeather.Clock(), track_numbers = true)

end
