"""$(TYPEDSIGNATURES) Main time loop that loops over all time steps."""
function time_stepping!(simulation::AbstractSimulation)
    (; clock) = simulation.variables.prognostic
    (; time_stepping) = simulation.model
    for _ in 1:clock.n_timesteps + spin_up_steps(time_stepping)
        time_step!(simulation)
    end
    return simulation
end

# dispatch over time stepper
time_step!(simulation::AbstractSimulation) = time_step!(simulation, simulation.model.time_stepping)

"""$(TYPEDSIGNATURES) Perform one time step of `simulation`."""
function time_step!(simulation::AbstractSimulation, time_stepping::AbstractTimeStepper)
    (; variables, model) = simulation
    (; feedback, output) = model
    (; clock) = variables.prognostic

    # TODO add a possible initialize implicit step here that defaults to nothing

    time_step!(variables, time_stepping, model)     # calculate tendencies and step forward
    time_step!(clock, time_stepping)                # then step the clock forward

    progress!(feedback, variables, model)           # updates the progress meter bar
    output!(output, simulation)                     # do output?
    callback!(model.callbacks, variables, model)    # any callbacks?
    return nothing
end

"""$(TYPEDSIGNATURES)
Calculate a single time step for the barotropic model."""
function time_step!(
        vars::Variables,                        # all variables
        time_stepping::AbstractTimeStepper,     # time stepping parameters
        model::Union{Barotropic, ShallowWater}, # everything that's constant at runtime
    )
    # exit immediately if NaNs/Infs already present
    (!isnothing(model.feedback) && model.feedback.nans_detected) && return nothing
    reset_tendencies!(vars, time_stepping)      # set the tendencies back to zero for accumulation

    dynamics_tendencies!(vars, model)
    diffusion_and_implicit!(vars, model)        # dispatch over time stepper and implicit so that the order can be changed
    update_prognostic!(vars, model)             # step prognostic variables forward
    transform!(vars, model)                     # new spectral state to grid
    particle_advection!(vars, model)            # TODO move up?

    return nothing
end

# dispatch over time stepper here so that other time stepper can change the order
diffusion_and_implicit!(vars, model) = 
    diffusion_and_implicit!(vars, model.time_stepping, model.implicit, model)

# implicit = nothing just call diffusion
diffusion_and_implicit!(vars, ::AbstractTimeStepper, ::Nothing, model) = horizontal_diffusion!(vars, model)

# default order is diffusion then implicit
function diffusion_and_implicit!(vars, ::AbstractTimeStepper, ::AbstractImplicit, model)
    horizontal_diffusion!(vars, model)
    implicit_correction!(vars, model)
    return nothing
end

"""
$(TYPEDSIGNATURES)
Calculate a single time step for the primitive equation model."""
function time_step!(
        vars::Variables,                # all variables
        dt::Real,                       # time step (mostly =2Δt, but for first_timesteps! =Δt, Δt/2)
        model::PrimitiveEquation,       # everything that's constant at runtime
        lf1::Integer = 2,               # leapfrog index 1 (dis/enables Robert+Williams filter)
        lf2::Integer = 2,               # leapfrog index 2 (time step used for tendencies)
    )
    # exit immediately if NaNs/Infs already present
    (!isnothing(model.feedback) && model.feedback.nans_detected) && return nothing
    reset_tendencies!(vars)             # set the tendencies back to zero for accumulation

    if ~model.dynamics_only             # switch on/off all physics parameterizations
        # calculate all parameterizations
        parameterization_tendencies!(vars, model)
        ocean_timestep!(vars, model)    # sea surface temperature and maybe in the future sea ice
        sea_ice_timestep!(vars, model)  # sea ice
        land_timestep!(vars, model)     # soil moisture and temperature, vegetation, maybe rivers
    end

    if model.dynamics                                       # switch on/off all dynamics
        dynamics_tendencies!(vars, lf2, model)              # dynamical core
        implicit_correction!(vars, model.implicit, model)   # semi-implicit time stepping corrections
    else    # just transform physics tendencies to spectral space
        parameterization_tendencies_only!(vars, model)
    end

    # APPLY DIFFUSION, STEP FORWARD IN TIME, AND TRANSFORM NEW TIME STEP TO GRID
    horizontal_diffusion!(vars, model)
    leapfrog!(vars, dt, lf1, model)
    transform!(vars, lf2, model)

    # PARTICLE ADVECTION (always skip 1st step of first_timesteps!)
    not_first_timestep = lf2 == 2
    not_first_timestep && particle_advection!(vars, model)

    return nothing
end

# dispatch via time stepping
update_prognostic!(var::Variables, model::AbstractModel) =
    update_prognostic!(var, model.time_stepping, model)

"""$(TYPEDSIGNATURES)
Leapfrog time stepping for all prognostic variables in `vars` using their tendencies.
Depending on `model` decides which variables to time step."""
function update_prognostic!(
        vars::Variables,
        time_stepping::AbstractTimeStepper,
        model::AbstractModel,
    )
    (; prognostic, tendencies) = vars

    count_step!(time_stepping)

    # atmospheric variables
    for varname in keys(tendencies)
        if !(tendencies[varname] isa NamedTuple)
            var = getfield(prognostic, varname)
            tendency = getfield(tendencies, varname)
            update_prognostic!(var, tendency, vars, time_stepping, model.implicit, model)
        end
    end

    # and time stepping for tracers if active
    for (name, tracer) in model.tracers
        if tracer.active
            var = prognostic.tracers[name]
            tendency = tendencies.tracers[name]
            update_prognostic!(var, tendency, vars, time_stepping, model.implicit, model)
        end
    end

    # ocean variables
    if haskey(tendencies, :ocean) && tendencies.ocean isa NamedTuple
        for varname in keys(tendencies.ocean)
            var = getfield(prognostic.ocean, varname)
            tendency = getfield(tendencies.ocean, varname)
            update_prognostic!(var, tendency, vars, time_stepping, model.implicit, model)
        end
    end

    # land variables
    if haskey(tendencies, :land) && tendencies.land isa NamedTuple
        for varname in keys(tendencies.land)
            var = getfield(prognostic.land, varname)
            tendency = getfield(tendencies.land, varname)
            update_prognostic!(var, tendency, vars, time_stepping, model.implicit, model)
        end
    end

    # evolve the random pattern in time
    random_process!(vars, model.random_process)
    return nothing
end