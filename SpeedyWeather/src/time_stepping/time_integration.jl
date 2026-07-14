"""$(TYPEDSIGNATURES) Main time loop that loops over all time steps."""
function time_stepping!(simulation::AbstractSimulation)
    (; clock) = simulation.variables.prognostic
    for _ in 1:clock.n_steps
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

    # re-initialize model components if needed, e.g. implicit with changing time step
    reinitialize!(model, variables)                      

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

    if model.dynamics
        dynamics_tendencies!(vars, model)
        diffusion_and_implicit!(vars, model)    # dispatch over time stepper and implicit so that the order can be changed
        update_prognostic!(vars, model)         # step prognostic variables forward
        transform!(vars, model)                 # new spectral state to grid
    end
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
        vars::Variables,                    # all variables
        time_stepping::AbstractTimeStepper, # dispatch over time stepper
        model::PrimitiveEquation,           # everything that's constant at runtime
    )
    # exit immediately if NaNs/Infs already present
    (!isnothing(model.feedback) && model.feedback.nans_detected) && return nothing
    reset_tendencies!(vars, time_stepping)  # set the tendencies back to zero for accumulation

    if ~model.dynamics_only                 # switch on/off all physics parameterizations
        greenhouse_gases_time_step!(vars, model)
        parameterization_tendencies!(vars, model)
        ocean_timestep!(vars, model)
        sea_ice_timestep!(vars, model)
        land_timestep!(vars, model)         # soil moisture and temperature,...
    end

    if model.dynamics                           # switch on/off all dynamics
        dynamics_tendencies!(vars, model)       # dynamical core
        diffusion_and_implicit!(vars, model)    # semi-implicit time stepping corrections
    else    # just transform physics tendencies to spectral space
        parameterization_tendencies_only!(vars, model)
        horizontal_diffusion!(vars, model)
    end

    update_prognostic!(vars, model)             # the actual time step
    transform!(vars, model)                     # transform spectral back to grid
    particle_advection!(vars, model)

    return nothing
end

# dispatch via time stepping
update_prognostic!(var::Variables, model::AbstractModel) =
    update_prognostic!(var, model.time_stepping, model)

"""$(TYPEDSIGNATURES)
Time stepping for all prognostic variables in `vars` using their tendencies.
Decides on which variables to time step based on the presence of their tendencies in `vars.tendencies`."""
function update_prognostic!(
        vars::Variables,
        time_stepping::AbstractTimeStepper,
        model::AbstractModel,
    )
    (; prognostic, tendencies) = vars
    (; clock) = prognostic
    scale = prognostic.scale[]

    # atmospheric variables
    for varname in tendency_names(vars)
        var = getfield(prognostic, varname)
        tendency = getfield(tendencies, varname)
        update_prognostic!(var, tendency, clock, time_stepping, model.implicit, model, scale)
    end

    # and time stepping for tracers if active
    for (name, tracer) in model.tracers
        if tracer.active
            var = prognostic.tracers[name]
            tendency = tendencies.tracers[name]
            update_prognostic!(var, tendency, clock, time_stepping, model.implicit, model, scale)
        end
    end

    # ocean variables, use unscaled time step
    for varname in ocean_tendency_names(vars)
        var = getfield(prognostic.ocean, varname)
        tendency = getfield(tendencies.ocean, varname)
        update_prognostic!(var, tendency, clock, time_stepping, model.implicit, model)
    end

    # land variables, use unscaled time step
    for varname in land_tendency_names(vars)
        var = getfield(prognostic.land, varname)
        tendency = getfield(tendencies.land, varname)
        update_prognostic!(var, tendency, clock, time_stepping, model.implicit, model)
    end

    # evolve the random pattern in time
    random_process!(vars, model.random_process)
    return nothing
end