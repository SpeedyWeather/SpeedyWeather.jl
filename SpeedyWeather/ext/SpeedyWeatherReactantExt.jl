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
r_first! = @compile SpeedyWeather.first_timesteps!(simulation)
r_later! = @compile SpeedyWeather.later_timestep!(simulation)
SpeedyWeather.time_stepping!(simulation, r_first!, r_later!)
SpeedyWeather.finalize!(simulation)
```
"""
function SpeedyWeather.time_stepping!(simulation::ReactantSimulation, r_first_timesteps! = nothing, r_later_timestep! = nothing, enable_checkpointing = true)

    clock = simulation.variables.prognostic.clock
    output = simulation.model.output

    if isnothing(r_first_timesteps!)
        @info "Reactant compiling first_timesteps!"
        r_first_timesteps! = @compile first_timesteps!(simulation)
    end

    #TODO: reenable @trace once Reactant issues fixed
    #@trace checkpointing = enable_checkpointing for _ in clock.timestep_counter:clock.n_timesteps
    #    r_later_timestep!(simulation)
    #end

    r_first_timesteps!(simulation)
    # Output is a no-op inside the compiled function (see overrides below); call it
    # explicitly here on the now-concrete simulation state so file-writing actually
    # happens for Reactant simulations.
    SpeedyWeather.output!(output, simulation)

    if isnothing(r_later_timestep!)
        @info "Reactant compiling later_timestep!"
        r_later_timestep! = @compile later_timestep!(simulation)
    end

    for i in (Int(clock.timestep_counter) + 1):Int(clock.n_timesteps)
        r_later_timestep!(simulation)
        SpeedyWeather.output!(output, simulation)
    end
    return
end

# that's for Reactant TracableDateTime
SpeedyWeather.secondofday(dt::ReactantDatesExt.ReactantDateTime) = Dates.second(convert(ReactantDatesExt.ReactantTime, dt).instant)

# For Reactant tracing, dispatch on the ReactantDateTime argument and forward
# to the internal `_year_angle`/`_solar_hour_angle` implementations. The
# `length_of_day`/`length_of_year` arguments may be regular `Second` (from the
# untraced SolarZenith struct fields) or Reactant period types, so we accept any.
# Both signatures are spelled out to avoid ambiguity with the SpeedyWeather methods
# that dispatch on `Dates.AbstractDateTime` + `Second`.
@inline SpeedyWeather.year_angle(::Type{T}, time::ReactantDatesExt.ReactantDateTime, length_of_day::Second, length_of_year::Second) where {T} = SpeedyWeather._year_angle(T, time, length_of_day, length_of_year)
@inline SpeedyWeather.year_angle(::Type{T}, time::ReactantDatesExt.ReactantDateTime, length_of_day::ReactantDatesExt.ReactantSecond, length_of_year::ReactantDatesExt.ReactantSecond) where {T} = SpeedyWeather._year_angle(T, time, length_of_day, length_of_year)
@inline SpeedyWeather.solar_hour_angle(::Type{T}, time::ReactantDatesExt.ReactantDateTime, λ, length_of_day::Second) where {T} = SpeedyWeather._solar_hour_angle(T, time, λ, length_of_day)
@inline SpeedyWeather.solar_hour_angle(::Type{T}, time::ReactantDatesExt.ReactantDateTime, λ, length_of_day::ReactantDatesExt.ReactantSecond) where {T} = SpeedyWeather._solar_hour_angle(T, time, λ, length_of_day)

Base.:-(x::ReactantDatesExt.ReactantDateTime, y::ReactantDatesExt.ReactantDateTime) =
    ReactantDatesExt.ReactantMillisecond(Dates.value(x) - Dates.value(y))

# These function extend those defined in SpeedyWeather/src/dynamics/clock.jl
# They will not move to ReactantDatesExt as they aren't part of stdlib Dates.jl
Dates.second(x::ReactantDatesExt.ReactantNanosecond) = round(Int, x.value * 1.0e-9)
Dates.second(x::ReactantDatesExt.ReactantMicrosecond) = round(Int, x.value * 1.0e-6)
Dates.second(x::ReactantDatesExt.ReactantMillisecond) = round(Int, x.value * 1.0e-3)

# construct components with ReactantDates types
const DATE_TYPE = Int64

SpeedyWeather.Clock(architecture::ReactantDevice) = Reactant.to_rarray(SpeedyWeather.Clock(), track_numbers = true)

SpeedyWeather.SolarZenith(SG::SpectralGrid{<:ReactantDevice}; kwargs...) = SolarZenith{SG.NF, SinSolarDeclination{typeof(Earth(SG))}, Base.RefValue{ReactantDatesExt.ReactantDateTime{DATE_TYPE}}, Bool, ReactantDatesExt.ReactantSecond{DATE_TYPE}}(; kwargs...)
SpeedyWeather.SolarZenithSeason(SG::SpectralGrid{<:ReactantDevice}; kwargs...) = SolarZenithSeason{SG.NF, SinSolarDeclination{typeof(Earth(SG))}, Base.RefValue{ReactantDatesExt.ReactantDateTime{DATE_TYPE}}, Bool,  ReactantDatesExt.ReactantSecond{DATE_TYPE}}(; kwargs...)
SpeedyWeather.Earth(SG::SpectralGrid{<:ReactantDevice}; kwargs...) = Earth{SG.NF, ReactantDatesExt.ReactantSecond{DATE_TYPE}, ReactantDatesExt.ReactantDateTime{DATE_TYPE}, Bool}(kwargs...)

function SpeedyWeather.WhichZenith(SG::SpectralGrid{<:ReactantDevice}, P::SpeedyWeather.AbstractPlanet; kwargs...)
    (; NF) = SG
    (; daily_cycle, seasonal_cycle, length_of_day, length_of_year) = P
    solar_declination = SpeedyWeather.SinSolarDeclination(P)

    if daily_cycle
        return SolarZenith{NF, typeof(solar_declination), Base.RefValue{ReactantDatesExt.ReactantDateTime{DATE_TYPE}}, Bool, ReactantDatesExt.ReactantSecond{DATE_TYPE}}(;
            length_of_day, length_of_year, solar_declination, seasonal_cycle, kwargs...
        )

    else
        return SolarZenithSeason{NF, typeof(solar_declination), Base.RefValue{ReactantDatesExt.ReactantDateTime{DATE_TYPE}}, Bool,  ReactantDatesExt.ReactantSecond{DATE_TYPE}}(;
            length_of_day, length_of_year, solar_declination, seasonal_cycle, kwargs...
        )
    end
end

Base.convert(::Type{Base.RefValue{ReactantDatesExt.ReactantDateTime{DATE_TYPE}}}, dt::Base.RefValue{DateTime}) = Base.RefValue{ReactantDatesExt.ReactantDateTime{DATE_TYPE}}(ReactantDatesExt.ReactantDateTime(dt[]))

# OUTPUT HANDLING FOR REACTANT
#
# TODO: This is a hacky solution that we might revise if we take output out of the model
# and into the `Simulation`
# 
# Reactant cannot trace through the output! pipeline: `output!` interpolates onto
# CPU scratch fields and writes to disk, both of which are runtime side effects
# the tracer can't (and shouldn't) capture. We split the behaviour by tracing
# context:
#   - Inside `@compile` (`within_compile()` == true): every output! method is a
#     hard no-op. Nothing gets baked into the compiled `first_timesteps!` /
#     `later_timestep!` graphs, so the compiled functions only do numerics.
#   - Outside `@compile`: `time_stepping!` (above) explicitly invokes
#     `output!(output, simulation)` after each compiled step. The simulation's
#     data is now ConcretePJRTArrays, which the standard CPU output path can
#     consume via `on_architecture(CPU(), …)` and the `ConcretePJRT → Array`
#     transfer rules in SpeedyWeatherInternalsReactantExt.
#
# We need overrides for both 2-arg variants of `output!` (the per-simulation
# dispatcher and the time-write helper); the 3-arg per-variable method does not
# need an override because it's only called from the 2-arg simulation method,
# which already gates on the tracing context.

# 2-arg: output!(output, simulation::ReactantSimulation)
# Inside trace → skip; outside trace → run the standard AbstractSimulation path.
function SpeedyWeather.output!(output::SpeedyWeather.AbstractOutput, simulation::ReactantSimulation)
    Reactant.ReactantCore.within_compile() && return nothing
    return @invoke SpeedyWeather.output!(output::SpeedyWeather.AbstractOutput, simulation::SpeedyWeather.AbstractSimulation)
end

# 2-arg: output!(output, time::ReactantDateTime)
# Inside trace → skip; outside trace → materialise into a plain DateTime and
# forward to the backend-specific DateTime method (NetCDF, Zarr, …).
function SpeedyWeather.output!(output::SpeedyWeather.AbstractOutput, time::ReactantDatesExt.ReactantDateTime)
    Reactant.ReactantCore.within_compile() && return nothing
    dt = DateTime(Dates.UTInstant(Millisecond(Dates.value(time))))
    return SpeedyWeather.output!(output, dt)
end

end
