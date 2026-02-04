module SpeedyWeatherReactantExt

using SpeedyWeather
using Reactant
using Dates, DocStringExtensions

using SpeedyWeather: ReactantDevice, scale!, get_step, unpack, timestep!, first_timesteps!, later_timestep!

const ReactantSimulation = Union{
    Simulation{<:BarotropicModel{SG, <:ReactantDevice}},
    Simulation{<:ShallowWaterModel{SG, <:ReactantDevice}}, Simulation{<:PrimitiveDryModel{SG, <:ReactantDevice}}, Simulation{<:PrimitiveWetModel{SG, <:ReactantDevice}},
} where {SG}

# initialize a BarotropicModel with ReactantDevice, only use @jit on some parts
# TODO: might define a custom @jit_or_not macro for this to use in main code
function SpeedyWeather.initialize!(model::BarotropicModel{SG, <:ReactantDevice}; time::DateTime = SpeedyWeather.DEFAULT_DATE) where {SG <: SpectralGrid}
    (; spectral_grid) = model

    spectral_grid.nlayers > 1 && @error "Only nlayers=1 supported for BarotropicModel, \
        SpectralGrid with nlayers=$(spectral_grid.nlayers) provided."

    # initialize components
    @jit initialize!(model.geometry, model)
    @jit initialize!(model.time_stepping, model)
    @jit initialize!(model.coriolis, model)
    @jit initialize!(model.forcing, model)
    @jit initialize!(model.drag, model)
    @jit initialize!(model.horizontal_diffusion, model)
    @jit initialize!(model.random_process, model)
    @jit initialize!(model.particle_advection, model)

    # allocate prognostic and diagnostic variables
    prognostic_variables = PrognosticVariables(model)
    diagnostic_variables = DiagnosticVariables(model)
    # initialize particles (or other non-atmosphere prognostic variables)
    initialize!(prognostic_variables.particles, prognostic_variables, diagnostic_variables, model)

    # set the initial conditions
    @jit initialize!(prognostic_variables, model.initial_conditions, model)
    (; clock) = prognostic_variables
    clock.time = time       # set the current time
    clock.start = time      # and store the start time

    return Simulation(prognostic_variables, diagnostic_variables, model)
end

function SpeedyWeather.initialize!(
        simulation::ReactantSimulation;
        period::Period = SpeedyWeather.DEFAULT_PERIOD,
        steps::Int = SpeedyWeather.DEFAULT_TIMESTEPS,
        output::Bool = false,
    )
    progn, diagn, model = unpack(simulation)

    # SET THE CLOCK
    (; clock) = progn
    (; time_stepping) = model
    if steps != SpeedyWeather.DEFAULT_TIMESTEPS
        # sets the steps, calculate period from it, store the start date, reset counter
        @assert period == SpeedyWeather.DEFAULT_PERIOD "Period and steps cannot be set simultaneously"
        initialize!(clock, time_stepping, steps)
    else
        # set period = how long to integrate for, tore the start date, reset counter
        initialize!(clock, time_stepping, period)
    end

    # OUTPUT, enable/disable output
    set!(simulation.model.output, active = output, reset_path = true)

    # SCALING: we use vorticity*radius, divergence*radius in the dynamical core
    @jit scale!(progn, diagn, model.planet.radius)

    # OUTPUT INITIALISATION AND STORING INITIAL CONDITIONS + FEEDBACK
    # propagate spectral state to grid variables for initial condition output
    lf = model.time_stepping.first_step_euler ? 1 : 2       # use 2nd leapfrog index when restarting

    # raise a warning if starting with leapfrog but there's zero vorticity
    vor = get_step(progn.vor, lf)
    lf == 2 && all(vor .== 0) && @warn "Vorticity is zero on 2nd leapfrog index though you use it to calculate tendencies." *
        " You may wanted to continue with a leapfrog step without data for it in the 2nd step."

    @jit transform!(diagn, progn, lf, model, initialize = true)
    initialize!(diagn, progn.particles, progn, model)
    initialize!(model.output, model.feedback, progn, diagn, model)
    return @jit initialize!(model.callbacks, progn, diagn, model)
end

function SpeedyWeather.time_stepping!(simulation::ReactantSimulation, r_first_timesteps! = nothing, r_later_timestep! = nothing)
    (; clock) = simulation.prognostic_variables
    for _ in 1:clock.n_timesteps        # MAIN LOOP
        timestep!(simulation, r_first_timesteps!, r_later_timestep!)
    end
    return
end

function SpeedyWeather.timestep!(simulation::ReactantSimulation, r_first_timesteps! = nothing, r_later_timestep! = nothing)
    if isnothing(r_first_timesteps!)
        @info "Reactant compiling first_timesteps!"
        r_first_timesteps! = @compile first_timesteps!(simulation)
    end

    if isnothing(r_later_timestep!)
        @info "Reactant compiling later_timestep!"
        r_later_timestep! = @compile later_timestep!(simulation)
    end

    return if clock.timestep_counter == 0
        r_first_timesteps!(simulation)
    else
        r_later_timestep!(simulation)
    end
end

end
