export Leapfrog

abstract type AbstractLeapfrog <: AbstractTimeStepper end

"""Leapfrog time stepping defined by the following fields
$(TYPEDFIELDS)"""
mutable struct Leapfrog{NF, S, B, MS} <: AbstractLeapfrog
    "[OPTION] Time step in minutes for T31, scale linearly to `trunc`"
    Δt_at_T31::S

    "[OPTION] Adjust `Δt_at_T31` with the `output_dt` to reach `output_dt` exactly in integer time steps"
    adjust_with_output::B

    "[OPTION] Start integration with (1) Euler step with dt/2, (2) Leapfrog step with dt"
    start_with_euler::B

    "[OPTION] Sets `first_step_euler=false` after first step to continue with leapfrog after 1st `run!` call"
    continue_with_leapfrog::B

    "[DERIVED] Use Euler on first time step? (controlled by `start_with_euler` and `continue_with_leapfrog`)"
    first_step_euler::B

    "[OPTION] Robert (1966) time filter coefficient to suppress the computational mode"
    robert_filter::NF

    "[OPTION] Williams time filter (Amezcua 2011) coefficient for 3rd order acc"
    williams_filter::NF

    "[DERIVED] Time step Δt in milliseconds at specified resolution"
    Δt_millisec::MS

    "[DERIVED] Time step Δt [s] at specified resolution"
    Δt_sec::NF

    "[DERIVED] Time step Δt [s/m] at specified resolution, scaled by 1/radius"
    Δt::NF
end

Adapt.adapt_structure(to, L::Leapfrog) = LeapfrogCore(L.Δt_millisec, L.Δt_sec, L.Δt, L.step_counter)

# HOW MANY STEPS DO VARIABLES NEED?
# leapfrogging always needs 2 steps in spectral
prognostic_spectral_steps(::AbstractLeapfrog) = 2
# but in 2D only 1 step in grid space
prognostic_grid_steps(::AbstractLeapfrog, ::Union{<:Barotropic, <:ShallowWater}) = 1
# but the parameterizations are evaluated at the previous step so 2
prognostic_grid_steps(::AbstractLeapfrog, ::PrimitiveEquation) = 2
# always only one step for tendencies
tendency_steps(::AbstractLeapfrog) = 1

# WHICH STEP TO READ WHEN
@inline which_prognostic_step(var::LowerTriangularArray, ::AbstractLeapfrog, ::SpeedyTransforms.AbstractSpectralTransform) = 2
@inline which_prognostic_step(var::LowerTriangularArray, ::AbstractLeapfrog, ::AbstractForcing) = 2
@inline which_prognostic_step(var::LowerTriangularArray, ::AbstractLeapfrog, ::AbstractDrag) = 2
@inline which_prognostic_step(var::LowerTriangularArray, ::AbstractLeapfrog, ::AbstractDynamicalCoreComponent) = 2
@inline which_prognostic_step(var::LowerTriangularArray, ::AbstractLeapfrog, ::AbstractHorizontalDiffusion) = 1

# dispatch over timestepper to decide between implicit or explicit diffusion
@inline implicit_diffusion(::AbstractHorizontalDiffusion, ::Nothing, ::AbstractLeapfrog) = true
@inline implicit_diffusion(::AbstractHorizontalDiffusion, ::AbstractImplicit, ::AbstractLeapfrog) = true

# copy prognostic variables' 1st step to 2nd, that way which_prognostic_step can always be 2
function initialize!(vars::Variables, ::AbstractLeapfrog, ::AbstractModel)
    (; prognostic) = vars
    for varname in keys(prognostic)
        if prognostic[varname] isa AbstractArray
            var_old, var_new = get_steps(getfield(prognostic, varname))
            var_new .= var_old
        end
    end
end

# for leapfrog do first semi-implicit corrections then horizontal diffusion
function diffusion_and_implicit!(vars, ::AbstractLeapfrog, ::AbstractImplicit, model)
    implicit_correction!(vars, model)
    horizontal_diffusion!(vars, model)
    return nothing
end

"""($TYPEDSIGNATURES)
Immutable core struct used to adapt only the time step fields for use in kernels"""
struct LeapfrogCore{NF, MS} <: AbstractLeapfrog
    "[DERIVED] Time step Δt in milliseconds at specified resolution"
    Δt_millisec::MS

    "[DERIVED] Time step Δt [s] at specified resolution"
    Δt_sec::NF

    "[DERIVED] Time step Δt [s/m] at specified resolution, scaled by 1/radius"
    Δt::NF
end

Adapt.@adapt_structure LeapfrogCore

"""$(TYPEDSIGNATURES)
Generator function for a Leapfrog struct using `spectral_grid`
for the resolution information."""
function Leapfrog(
        spectral_grid::SpectralGrid;
        Δt_at_T31 = Minute(40),
        adjust_with_output = true,
        start_with_euler = true,
        continue_with_leapfrog = true,
        robert_filter = 0.1,
        williams_filter = 0.53,
        radius = DEFAULT_RADIUS,
    )
    (; NF, trunc) = spectral_grid

    # compute time step
    Δt_millisec::Millisecond = get_Δt_millisec(Second(Δt_at_T31), trunc, DEFAULT_RADIUS, adjust_with_output)
    Δt_sec::NF = Δt_millisec.value / 1000
    Δt::NF = Δt_sec / radius

    # derived and mutated, controlled by start_with_euler and continue_with_leapfrog
    first_step_euler = start_with_euler

    return Leapfrog(
        Second(Δt_at_T31), adjust_with_output, start_with_euler, continue_with_leapfrog, first_step_euler,
        NF(robert_filter), NF(williams_filter), Δt_millisec, Δt_sec, Δt,
    )
end

"""$(TYPEDSIGNATURES)
Initialize leapfrogging `L` by recalculating the time step given the output time step
`output_dt` from `model.output`. Recalculating will slightly adjust the time step to
be a divisor such that an integer number of time steps matches exactly with the output
time step."""
function initialize!(L::Leapfrog, model::AbstractModel)
    calculate_Δt!(L, model)  # common among several time steppers
    if L.start_with_euler
        L.first_step_euler = true
    end
    return nothing
end

"""$(TYPEDSIGNATURES) Leapfrog is spun up with 1 Euler forward step that doesn't count for clock + output"""
spin_up_steps(::Leapfrog) = 1

function time_step!(clock::Clock, time_stepping::Leapfrog)
    Δt = time_stepping.Δt_millisec  # ::Millisecond, integer based hence ÷ not / below
    i = clock.step_counter          # 0-based as the clock is only stepped below
    if i == 0                       # first Euler step at Δt/2
        # i counts every time step, for the clock the first Euler step does not count
        # hence after this the time_stepping will be 1 ahead of clock step counter
        time_step!(clock, Δt ÷ 2, increase_counter = false)
    elseif i == 1                   # second step: Leapfrog at Δt
        # subtract the Δt/2 again as otherwise the time can be 1ms off due to rounding
        clock.time -= Δt ÷ 2
        time_step!(clock, Δt)
    else                            # later steps: Leapfrog at 2Δt but increase clock by Δt
        time_step!(clock, Δt)
    end
    return nothing
end

function time_step(L::Leapfrog, clock::Clock)
    (; Δt) = L
    clock.step_counter == 0 && return Δt / 2    # first step Euler with Δt/2
    clock.step_counter == 1 && return Δt        # 2nd step leapfrog with Δt
    return 2Δt                                  # later steps leapfrog with 2Δt
end

function prognostic_step(::Leapfrog, clock::Clock)
    clock.step_counter == 0 && return 1         # first Euler step disable filters
    clock.step_counter == 1 && return 1         # 2nd step: Leapfrog also disable filters
    return 2                                    # later steps: with RAW filters
end

function update_prognostic!(
        var::AbstractArray,
        tendency::AbstractArray,
        vars::Variables,
        time_stepping::Leapfrog,
        implicit::Union{Nothing, AbstractImplicit},
        ::AbstractModel,
    )
    (; clock) = vars.prognostic
    Δt = time_step(time_stepping, clock)
    lf = prognostic_step(time_stepping, clock)  # leapfrog prognostic step index
    var_old, var_new = get_steps(var)
    var_lf = get_step(var, lf)                  # view on either t or t+dt to dis/enable Williams filter
    var_tend = get_tendency_step(tendency, time_stepping, time_stepping)

    @boundscheck lf == 1 || lf == 2 || throw(BoundsError())
    @boundscheck size(var_old) == size(var_new) == size(var_tend) || throw(BoundsError())

    (; robert_filter, williams_filter) = time_stepping          # coefficients for filters

    # LEAP FROG time step with or without Robert+Williams filter
    # Robert time filter to compress computational mode, Williams filter for 3rd order accuracy
    # see Williams (2009), Eq. 7-9
    # for lf == 1 (time steps 1 or 2) no filters applied (w1=w2=0)
    # for lf == 2 (later steps) Robert+Williams filter is applied
    w1 = (lf - 1) * robert_filter * williams_filter / 2         # = ν*α/2 in Williams (2009, Eq. 8)
    w2 = (lf - 1) * robert_filter * (1 - williams_filter) / 2   # = ν(1-α)/2 in Williams (2009, Eq. 9)

    launch!(
        architecture(tendency), SpectralWorkOrder, size(tendency), leapfrog_kernel!,
        var_old, var_new, var_lf, tendency, Δt, w1, w2
    )
    return nothing
end

@kernel inbounds = true function leapfrog_kernel!(var_old, var_new, var_lf, tendency, Δt, w1, w2)
    lmk = @index(Global, Linear)    # every harmonic lm, every vertical layer k
    old = var_old[lmk]
    new = old + Δt * tendency[lmk]
    update = old - 2var_lf[lmk] + new
    var_old[lmk] = var_lf[lmk] + w1 * update
    var_new[lmk] = new - w2 * update
end

# """$(TYPEDSIGNATURES)
# Perform one single time step of `simulation` including
# possibly output and callbacks."""
# function timestep!(simulation::AbstractSimulation, ::Leapfrog)
#     (; clock) = simulation.variables.prognostic
#     @trace if clock.step_counter == 0
#         leapfrog_first_timesteps!(simulation)
#     else
#         later_timestep!(simulation)
#     end

#     return nothing
# end

# """$(TYPEDSIGNATURES)
# First 1 or 2 time steps of `simulation`. If `model.time_stepping.start_with_euler` is true,
# then start with one Euler step with dt/2, followed by one Leapfrog step with dt.
# If false, continue with leapfrog steps at 2Δt (e.g. restart)."""
# function leapfrog_first_timesteps!(simulation::AbstractSimulation)
#     (; variables, model) = simulation
#     (; time_stepping) = model
#     (; Δt) = time_stepping

#     # decide whether to start with 1x Euler then 1x Leapfrog at Δt
#     @trace if time_stepping.first_step_euler
#         leapfrog_first_timesteps!(variables, model)
#         time_stepping.first_step_euler = !time_stepping.continue_with_leapfrog   # after first run! continue with leapfrog
#     else # or continue with leaprog steps at 2Δt (e.g. restart)
#         # but make sure that implicit solver is initialized in that situation
#         initialize!(model.implicit, 2Δt, variables, model)
#         set_initialized!(model.implicit)            # mark implicit as initialized
#         leapfrog_later_timestep!(simulation)
#     end

#     # only now initialise feedback for benchmark accuracy
#     (; clock) = variables.prognostic
#     initialize!(model.feedback, clock, model)
#     return nothing
# end

# """
# $(TYPEDSIGNATURES)
# Performs the first two initial time steps (Euler forward, unfiltered leapfrog) to populate the
# prognostic variables with two time steps (t=0, Δt) that can then be used in the normal leap frogging."""
# function first_timesteps!(
#         vars::Variables,        # all variables
#         model::AbstractModel,   # everything that is constant at runtime
#     )
#     (; clock) = vars.prognostic
#     # TODO: deactaved that check for Reactant, but it doesn't seem neccary anyway as the toplevle time_stepping! should take care of this
#     # clock.n_time_steps == 0 && return nothing    # exit immediately for no time steps

#     (; implicit) = model
#     (; Δt, Δt_millisec) = model.time_stepping
#     Δt_millisec_half = Millisecond(Δt_millisec.value ÷ 2)   # this might be 1ms off

#     # FIRST TIME STEP (EULER FORWARD with dt=Δt/2)
#     lf1 = 1                             # without Robert+Williams filter
#     lf2 = 1                             # evaluates all tendencies at t=0,
#     # the first leapfrog index (=>Euler forward)
#     initialize!(implicit, Δt / 2, vars, model)      # update precomputed implicit terms with time step Δt/2
#     timestep!(vars, Δt / 2, model, lf1, lf2)        # update time by half the leapfrog time step Δt used here
#     timestep!(clock, Δt_millisec_half, increase_counter = false)

#     # output, callbacks not called after the first Euler step as it's a half-step (i=0 to i=1/2)
#     # populating the second leapfrog index to perform the second time step

#     # SECOND TIME STEP (UNFILTERED LEAPFROG with dt=Δt, leapfrogging from t=0 over t=Δt/2 to t=Δt)
#     initialize!(implicit, Δt, vars, model)  # update precomputed implicit terms with time step Δt
#     lf1 = 1                                 # without Robert+Williams filter
#     lf2 = 2                                 # evaluate all tendencies at t=dt/2,
#     # the 2nd leapfrog index (=>Leapfrog)
#     timestep!(vars, Δt, model, lf1, lf2)
#     # remove prev Δt/2 in case not even milliseconds, otherwise time is off by 1ms
#     timestep!(clock, -Δt_millisec_half, increase_counter = false)
#     timestep!(clock, Δt_millisec)

#     # do output and callbacks after the first proper (from i=0 to i=1) time step
#     callback!(model.callbacks, vars, model)
#     output!(model.output, Simulation(vars, model))

#     # from now on precomputed implicit terms with 2Δt
#     initialize!(implicit, 2Δt, vars, model)
#     set_initialized!(implicit)      # mark implicit as initialized

#     return nothing
# end
