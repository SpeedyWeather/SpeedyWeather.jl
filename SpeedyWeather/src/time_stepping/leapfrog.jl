const DEFAULT_NSTEPS = 2
export Leapfrog

"""Leapfrog time stepping defined by the following fields
$(TYPEDFIELDS)"""
@kwdef mutable struct Leapfrog{NF, IntType, S, MS, B} <: AbstractTimeStepper
    "[DERIVED] Spectral resolution (max degree of spherical harmonics)"
    trunc::IntType

    "[CONST] Number of time steps stored simultaneously in prognostic variables"
    nsteps::IntType = 2

    "[OPTION] Time step in minutes for T31, scale linearly to `trunc`"
    Δt_at_T31::S = Minute(40)

    "[OPTION] Adjust `Δt_at_T31` with the `output_dt` to reach `output_dt` exactly in integer time steps"
    adjust_with_output::B = true

    "[OPTION] Start integration with (1) Euler step with dt/2, (2) Leapfrog step with dt"
    start_with_euler::B = true

    "[OPTION] Sets `first_step_euler=false` after first step to continue with leapfrog after 1st `run!` call"
    continue_with_leapfrog::B = true

    "[DERIVED] Use Euler on first time step? (controlled by `start_with_euler` and `continue_with_leapfrog`)"
    first_step_euler::B = start_with_euler

    "[OPTION] Robert (1966) time filter coefficient to suppress the computational mode"
    robert_filter::NF = 0.1

    "[OPTION] Williams time filter (Amezcua 2011) coefficient for 3rd order acc"
    williams_filter::NF = 0.53

    "[DERIVED] Radius of sphere [m], used for scaling, set in `initialize!` to `planet.radius`"
    radius::NF = DEFAULT_RADIUS

    "[DERIVED] Time step Δt in milliseconds at specified resolution"
    Δt_millisec::MS = get_Δt_millisec(Second(Δt_at_T31), trunc, radius, adjust_with_output)

    "[DERIVED] Time step Δt [s] at specified resolution"
    Δt_sec::NF = Δt_millisec.value / 1000

    "[DERIVED] Time step Δt [s/m] at specified resolution, scaled by 1/radius"
    Δt::NF = Δt_sec / radius
end

function Adapt.adapt_structure(to, L::Leapfrog)
    return (; Δt = L.Δt, Δt_sec = L.Δt_sec, Δt_millisec = L.Δt_millisec)
end

get_prognostic_steps(::Leapfrog) = 2
get_tendency_steps(::Leapfrog) = 1


"""$(TYPEDSIGNATURES)
Generator function for a Leapfrog struct using `spectral_grid`
for the resolution information."""
function Leapfrog(spectral_grid::SpectralGrid; kwargs...)
    (; NF, trunc) = spectral_grid
    return Leapfrog{NF, typeof(trunc), Dates.Second, Dates.Millisecond, Bool}(; trunc, kwargs...)
end

"""$(TYPEDSIGNATURES)
Initialize leapfrogging `L` by recalculating the time step given the output time step
`output_dt` from `model.output`. Recalculating will slightly adjust the time step to
be a divisor such that an integer number of time steps matches exactly with the output
time step."""
function initialize!(L::Leapfrog, model::AbstractModel)
    (; radius) = model.planet
    output_dt = get_output_dt(model.output)

    # take radius from planet and recalculate time step and possibly adjust with output dt
    L.Δt_millisec = get_Δt_millisec(L.Δt_at_T31, L.trunc, radius, L.adjust_with_output, output_dt)
    L.Δt_sec = L.Δt_millisec.value / 1000
    L.Δt = L.Δt_sec / radius

    # check how time steps from time integration and output align
    n = round(Int, Millisecond(output_dt).value / L.Δt_millisec.value)
    nΔt = n * L.Δt_millisec
    if nΔt != output_dt
        @warn "$n steps of Δt = $(L.Δt_millisec.value)ms yield output every " *
            "$(nΔt.value)ms (=$(nΔt.value / 1000)s), but output_dt = $(output_dt.value)s"
    end
    if L.start_with_euler
        L.first_step_euler = true
    end

    return nothing
end

"""$(TYPEDSIGNATURES)
Performs one leapfrog time step with (`lf=2`) or without (`lf=1`) Robert+Williams filter
(see Williams (2009), Montly Weather Review, Eq. 7-9)."""
function leapfrog!(
        A_old::LowerTriangularArray,        # prognostic variable at t
        A_new::LowerTriangularArray,        # prognostic variable at t+dt
        tendency::LowerTriangularArray,     # tendency (dynamics+physics) of A
        dt::Real,                           # time step (=2Δt, but for init steps =Δt, Δt/2)
        lf::Int,                            # leapfrog index to dis/enable Williams filter
        L::Leapfrog{NF},                    # struct with constants
    ) where {NF}                              # number format NF

    @boundscheck lf == 1 || lf == 2 || throw(BoundsError())         # index lf picks leapfrog dim
    @boundscheck size(A_old) == size(A_new) == size(tendency) || throw(BoundsError())

    A_lf = lf == 1 ? A_old : A_new              # view on either t or t+dt to dis/enable Williams filter
    (; robert_filter, williams_filter) = L      # coefficients for the Robert and Williams filter
    dt_NF = convert(NF, dt)                     # time step dt in number format NF

    # LEAP FROG time step with or without Robert+Williams filter
    # Robert time filter to compress computational mode, Williams filter for 3rd order accuracy
    # see Williams (2009), Eq. 7-9
    # for lf == 1 (initial time step) no filter applied (w1=w2=0)
    # for lf == 2 (later steps) Robert+Williams filter is applied
    w1 = lf == 1 ? zero(NF) : robert_filter * williams_filter / 2       # = ν*α/2 in Williams (2009, Eq. 8)
    w2 = lf == 1 ? zero(NF) : robert_filter * (1 - williams_filter) / 2 # = ν(1-α)/2 in Williams (2009, Eq. 9)

    launch!(
        architecture(tendency), SpectralWorkOrder, size(tendency), leapfrog_kernel!,
        A_old, A_new, A_lf, tendency, dt_NF, w1, w2
    )

    return nothing
end

@kernel inbounds = true function leapfrog_kernel!(A_old, A_new, A_lf, tendency, dt, w1, w2)

    lmk = @index(Global, Linear)    # every harmonic lm, every vertical layer k

    a_old = A_old[lmk]
    a_new = a_old + dt * tendency[lmk]
    a_update = a_old - 2A_lf[lmk] + a_new
    A_old[lmk] = A_lf[lmk] + w1 * a_update
    A_new[lmk] = a_new - w2 * a_update
end

# variables that are leapfrogged in the respective models, e.g. :vor_tend, :div_tend, etc...
tendency_names(model::AbstractModel) = tuple((Symbol(var, :_tend) for var in prognostic_variables(model))...)

"""$(TYPEDSIGNATURES)
Leapfrog time stepping for all prognostic variables in `vars` using their tendencies.
Depending on `model` decides which variables to time step."""
function leapfrog!(
        vars::Variables,
        dt::Real,               # time step (mostly =2Δt, but for init steps =Δt, Δt/2)
        lf::Int,                # leapfrog index to dis/enable Williams filter
        model::AbstractModel,
    )
    (; prognostic, tendencies) = vars

    for varname in keys(tendencies)
        if !(tendencies[varname] isa NamedTuple)
            var = getfield(prognostic, varname)
            var_old, var_new = get_steps(var)
            var_tend = getfield(tendencies, varname)
            SpeedyTransforms.spectral_truncation!(var_tend)
            leapfrog!(var_old, var_new, var_tend, dt, lf, model.time_stepping)
        end
    end

    # and time stepping for tracers if active
    for (name, tracer) in model.tracers
        if tracer.active
            var_old, var_new = get_steps(prognostic.tracers[name])
            var_tend = tendencies.tracers[name]
            SpeedyTransforms.spectral_truncation!(var_tend)
            leapfrog!(var_old, var_new, var_tend, dt, lf, model.time_stepping)
        end
    end

    # evolve the random pattern in time
    random_process!(vars, model.random_process)
    return nothing
end

"""$(TYPEDSIGNATURES)
First 1 or 2 time steps of `simulation`. If `model.time_stepping.start_with_euler` is true,
then start with one Euler step with dt/2, followed by one Leapfrog step with dt.
If false, continue with leapfrog steps at 2Δt (e.g. restart)."""
function first_timesteps!(simulation::AbstractSimulation)
    (; variables, model) = simulation
    (; time_stepping) = model
    (; Δt) = time_stepping

    # decide whether to start with 1x Euler then 1x Leapfrog at Δt
    @trace if time_stepping.first_step_euler
        first_timesteps!(variables, model)
        time_stepping.first_step_euler = !time_stepping.continue_with_leapfrog   # after first run! continue with leapfrog
    else # or continue with leaprog steps at 2Δt (e.g. restart)
        # but make sure that implicit solver is initialized in that situation
        initialize!(model.implicit, 2Δt, variables, model)
        set_initialized!(model.implicit)            # mark implicit as initialized
        later_timestep!(simulation)
    end

    # only now initialise feedback for benchmark accuracy
    (; clock) = variables.prognostic
    initialize!(model.feedback, clock, model)
    return nothing
end

"""
$(TYPEDSIGNATURES)
Performs the first two initial time steps (Euler forward, unfiltered leapfrog) to populate the
prognostic variables with two time steps (t=0, Δt) that can then be used in the normal leap frogging."""
function first_timesteps!(
        vars::Variables,        # all variables
        model::AbstractModel,   # everything that is constant at runtime
    )
    (; clock) = vars.prognostic
    # TODO: deactaved that check for Reactant, but it doesn't seem neccary anyway as the toplevle time_stepping! should take care of this
    #clock.n_timesteps == 0 && return nothing    # exit immediately for no time steps

    (; implicit) = model
    (; Δt, Δt_millisec) = model.time_stepping
    Δt_millisec_half = Millisecond(Δt_millisec.value ÷ 2)   # this might be 1ms off

    # FIRST TIME STEP (EULER FORWARD with dt=Δt/2)
    lf1 = 1                             # without Robert+Williams filter
    lf2 = 1                             # evaluates all tendencies at t=0,
    # the first leapfrog index (=>Euler forward)
    initialize!(implicit, Δt / 2, vars, model)      # update precomputed implicit terms with time step Δt/2
    timestep!(vars, Δt / 2, model, lf1, lf2)        # update time by half the leapfrog time step Δt used here
    timestep!(clock, Δt_millisec_half, increase_counter = false)

    # output, callbacks not called after the first Euler step as it's a half-step (i=0 to i=1/2)
    # populating the second leapfrog index to perform the second time step

    # SECOND TIME STEP (UNFILTERED LEAPFROG with dt=Δt, leapfrogging from t=0 over t=Δt/2 to t=Δt)
    initialize!(implicit, Δt, vars, model)  # update precomputed implicit terms with time step Δt
    lf1 = 1                                 # without Robert+Williams filter
    lf2 = 2                                 # evaluate all tendencies at t=dt/2,
    # the 2nd leapfrog index (=>Leapfrog)
    timestep!(vars, Δt, model, lf1, lf2)
    # remove prev Δt/2 in case not even milliseconds, otherwise time is off by 1ms
    timestep!(clock, -Δt_millisec_half, increase_counter = false)
    timestep!(clock, Δt_millisec)

    # do output and callbacks after the first proper (from i=0 to i=1) time step
    callback!(model.callbacks, vars, model)
    output!(model.output, Simulation(vars, model))

    # from now on precomputed implicit terms with 2Δt
    initialize!(implicit, 2Δt, vars, model)
    set_initialized!(implicit)      # mark implicit as initialized

    return nothing
end