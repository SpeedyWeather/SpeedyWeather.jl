module SpeedyWeatherReactantExt

using SpeedyWeather
using Reactant
using Dates, DocStringExtensions

using SpeedyWeather: ReactantDevice, scale!, get_step, unpack, timestep!, first_timesteps!, later_timestep!

const ReactantSimulation = Union{
    Simulation{<:BarotropicModel{SG, <:ReactantDevice}},
    Simulation{<:ShallowWaterModel{SG, <:ReactantDevice}}, Simulation{<:PrimitiveDryModel{SG, <:ReactantDevice}}, Simulation{<:PrimitiveWetModel{SG, <:ReactantDevice}},
} where {SG}

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

function SpeedyWeather.run!(
        simulation::ReactantSimulation, 
        r_first_timesteps! = nothing, 
        r_later_timestep! = nothing;
        period::Period = SpeedyWeather.DEFAULT_PERIOD,
        steps::Int = SpeedyWeather.DEFAULT_TIMESTEPS,
        output::Bool = false,
    )
    SpeedyWeather.initialize!(simulation; period, steps, output)                      # scaling, initialize output, store initial conditions
    SpeedyWeather.time_stepping!(simulation, r_first_timesteps!, r_later_timestep!)   # run it, yeah!
    SpeedyWeather.finalize!(simulation)                                               # unscale, finalize output, write restart file, finalize callbacks
    return SpeedyWeather.unicodeplot(simulation)
end

end
