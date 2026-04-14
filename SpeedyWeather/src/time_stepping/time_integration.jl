"""$(TYPEDSIGNATURES)
Computes the time step in [ms]. `Δt_at_T31` is always scaled with the resolution `trunc` 
of the model. In case `adjust_Δt_with_output` is true, the `Δt_at_T31` is additionally 
adjusted to the closest divisor of `output_dt` so that the output time axis is keeping
`output_dt` exactly."""
function get_Δt_millisec(
        Δt_at_T31::TimePeriod,
        trunc,
        radius,
        adjust_with_output::Bool,
        output_dt::TimePeriod = DEFAULT_OUTPUT_DT,
    )
    # linearly scale Δt with trunc+1 (which are often powers of two)
    resolution_factor = (DEFAULT_TRUNC + 1) / (trunc + 1)

    # radius also affects grid spacing, scale proportionally
    radius_factor = radius / DEFAULT_RADIUS

    # maybe rename to _at_trunc_and_radius?
    Δt_at_trunc = Second(Δt_at_T31).value * resolution_factor * radius_factor

    if adjust_with_output && (output_dt > Millisecond(0))
        k = round(Int, Second(output_dt).value / Δt_at_trunc)
        divisors = Primes.divisors(Millisecond(output_dt).value)
        sort!(divisors)
        i = findfirst(x -> x >= k, divisors)
        k_new = isnothing(i) ? k : divisors[i]
        Δt_millisec = Millisecond(round(Int, Millisecond(output_dt).value / k_new))

        # provide info when time step is significantly shortened or lengthened
        Δt_millisec_unadjusted = round(Int, 1000 * Δt_at_trunc)
        Δt_ratio = Δt_millisec.value / Δt_millisec_unadjusted

        if abs(Δt_ratio - 1) > 0.05     # print info only when +-5% changes
            p = round(Int, (Δt_ratio - 1) * 100)
            ps = p > 0 ? "+" : ""
            @info "Time step changed from $Δt_millisec_unadjusted to $Δt_millisec ($ps$p%) to match output frequency."
        end
    else
        Δt_millisec = Millisecond(round(Int, 1000 * Δt_at_trunc))
    end

    return Δt_millisec
end

"""$(TYPEDSIGNATURES)
Change time step of timestepper `L` to `Δt` (unscaled)
and disables adjustment to output frequency."""
function set!(
        L::AbstractTimeStepper,
        Δt::Period,                 # unscaled time step in Second, Minute, ...
    )
    L.Δt_millisec = Millisecond(Δt)         # recalculate all Δt fields
    L.Δt_sec = L.Δt_millisec.value / 1000
    L.Δt = L.Δt_sec / L.radius

    # recalculate the default time step at resolution T31 to be consistent
    resolution_factor = (L.trunc + 1) / (DEFAULT_TRUNC + 1)
    L.Δt_at_T31 = Second(round(Int, L.Δt_sec * resolution_factor))

    # given Δt was manually set disallow adjustment to output frequency
    L.adjust_with_output = false
    return L
end

# also allow for keyword arguments
set!(L::AbstractTimeStepper; Δt::Period) = set!(L, Δt)

"""$(TYPEDSIGNATURES)
Calculate a single time step for the barotropic model."""
function timestep!(
        vars::Variables,                # all variables
        dt::Real,                       # time step (mostly =2Δt, but for first_timesteps! =Δt, Δt/2)
        model::Barotropic,              # everything that's constant at runtime
        lf1::Integer = 2,               # leapfrog index 1 (dis/enables Robert+Williams filter)
        lf2::Integer = 2,               # leapfrog index 2 (time step used for tendencies)
    )
    # exit immediately if NaNs/Infs already present
    (!isnothing(model.feedback) && model.feedback.nans_detected) && return nothing

    reset_tendencies!(vars)             # set the tendencies back to zero for accumulation

    # TENDENCIES, DIFFUSION, LEAPFROGGING AND TRANSFORM SPECTRAL STATE TO GRID
    dynamics_tendencies!(vars, lf2, model)
    horizontal_diffusion!(vars, model.horizontal_diffusion, model)
    leapfrog!(vars, dt, lf1, model)
    transform!(vars, lf2, model)

    # PARTICLE ADVECTION (always skip 1st step of first_timesteps!)
    not_first_timestep = lf2 == 2
    not_first_timestep && particle_advection!(vars, model)

    return nothing
end

"""
$(TYPEDSIGNATURES)
Calculate a single time step for the `model <: ShallowWater`."""
function timestep!(
        vars::Variables,                # all variables
        dt::Real,                       # time step (mostly =2Δt, but for first_timesteps! =Δt, Δt/2)
        model::ShallowWater,            # everything that's constant at runtime
        lf1::Integer = 2,               # leapfrog index 1 (dis/enables Robert+Williams filter)
        lf2::Integer = 2,               # leapfrog index 2 (time step used for tendencies)
    )
    # exit immediately if NaNs/Infs already present
    (!isnothing(model.feedback) && model.feedback.nans_detected) && return nothing
    reset_tendencies!(vars)             # set the tendencies back to zero for accumulation

    # GET TENDENCIES, CORRECT THEM FOR SEMI-IMPLICIT INTEGRATION
    dynamics_tendencies!(vars, lf2, model)
    implicit_correction!(vars, model.implicit, model)

    # APPLY DIFFUSION, STEP FORWARD IN TIME, AND TRANSFORM NEW TIME STEP TO GRID
    horizontal_diffusion!(vars, model.horizontal_diffusion, model)
    leapfrog!(vars, dt, lf1, model)
    transform!(vars, lf2, model)

    # PARTICLE ADVECTION (always skip 1st step of first_timesteps!)
    not_first_timestep = lf2 == 2
    not_first_timestep && particle_advection!(vars, model)

    return nothing
end

"""
$(TYPEDSIGNATURES)
Calculate a single time step for the `model<:PrimitiveEquation`"""
function timestep!(
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
    horizontal_diffusion!(vars, model.horizontal_diffusion, model)
    leapfrog!(vars, dt, lf1, model)
    transform!(vars, lf2, model)

    # PARTICLE ADVECTION (always skip 1st step of first_timesteps!)
    not_first_timestep = lf2 == 2
    not_first_timestep && particle_advection!(vars, model)

    return nothing
end

"""$(TYPEDSIGNATURES)
Perform one single time step of `simulation` including
possibly output and callbacks."""
function timestep!(simulation::AbstractSimulation)
    (; clock) = simulation.variables.prognostic
    @trace if clock.timestep_counter == 0
        first_timesteps!(simulation)
    else
        later_timestep!(simulation)
    end

    return nothing
end

"""$(TYPEDSIGNATURES)
Perform one single "normal" time step of `simulation`, after `first_timesteps!`."""
function later_timestep!(simulation::AbstractSimulation)
    (; variables, model) = simulation
    (; feedback, output) = model
    (; Δt, Δt_millisec) = model.time_stepping
    (; clock) = variables.prognostic

    timestep!(variables, 2Δt, model)                # calculate tendencies and leapfrog forward
    timestep!(clock, Δt_millisec)                   # time of lf=2 and variables after timestep!

    progress!(feedback, variables)                  # updates the progress meter bar
    output!(output, simulation)                     # do output?
    callback!(model.callbacks, variables, model)    # any callbacks?
    return nothing
end

"""$(TYPEDSIGNATURES)
Main time loop that loops over all time steps."""
function time_stepping!(simulation::AbstractSimulation)
    (; clock) = simulation.variables.prognostic
    for _ in 1:clock.n_timesteps        # MAIN LOOP
        timestep!(simulation)
    end
    return simulation
end
