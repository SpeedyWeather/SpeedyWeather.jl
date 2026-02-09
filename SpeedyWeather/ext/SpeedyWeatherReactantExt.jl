module SpeedyWeatherReactantExt

using SpeedyWeather
using Reactant
using Dates, DocStringExtensions

using SpeedyWeather: ReactantDevice, scale!, get_step, unpack, timestep!, first_timesteps!, later_timestep!

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
function SpeedyWeather.time_stepping!(simulation::ReactantSimulation, r_first_timesteps! = nothing, r_later_timestep! = nothing)
    if isnothing(r_first_timesteps!)
        @info "Reactant compiling first_timesteps!"
        r_first_timesteps! = @compile first_timesteps!(simulation)
    end

    if isnothing(r_later_timestep!)
        @info "Reactant compiling later_timestep!"
        r_later_timestep! = @compile later_timestep!(simulation)
    end

    (; clock) = simulation.prognostic_variables
    for _ in 1:clock.n_timesteps        # MAIN LOOP
        timestep!(simulation, r_first_timesteps!, r_later_timestep!)
    end
    return
end

function SpeedyWeather.timestep!(simulation::ReactantSimulation, r_first_timesteps!, r_later_timestep!)
    (; clock) = simulation.prognostic_variables

    return if clock.timestep_counter == 0
        r_first_timesteps!(simulation)
    else
        r_later_timestep!(simulation)
    end
end

end
