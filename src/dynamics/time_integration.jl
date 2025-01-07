const DEFAULT_NSTEPS = 2
export Leapfrog

"""
Leapfrog time stepping defined by the following fields
$(TYPEDFIELDS)
"""
@kwdef mutable struct Leapfrog{NF<:AbstractFloat} <: AbstractTimeStepper

    # DIMENSIONS
    "spectral resolution (max degree of spherical harmonics)"
    trunc::Int                      

    "Number of timesteps stored simultaneously in prognostic variables"
    nsteps::Int = 2

    # OPTIONS
    "Time step in minutes for T31, scale linearly to `trunc`"
    Δt_at_T31::Second = Minute(40)

    "Radius of sphere [m], used for scaling"
    radius::NF = DEFAULT_RADIUS

    "Adjust Δt_at_T31 with the output_dt to reach output_dt exactly in integer time steps"
    adjust_with_output::Bool = true

    # NUMERICS
    "Robert (1966) time filter coefficeint to suppress comput. mode"
    robert_filter::NF = 0.1

    "Williams time filter (Amezcua 2011) coefficient for 3rd order acc"
    williams_filter::NF = 0.53

    # DERIVED FROM OPTIONS    
    "time step Δt [ms] at specified resolution"
    Δt_millisec::Millisecond = get_Δt_millisec(Second(Δt_at_T31), trunc, radius, adjust_with_output)

    "time step Δt [s] at specified resolution"
    Δt_sec::NF = Δt_millisec.value/1000

    "time step Δt [s/m] at specified resolution, scaled by 1/radius"
    Δt::NF = Δt_sec/radius  
end

"""
$(TYPEDSIGNATURES)
Computes the time step in [ms]. `Δt_at_T31` is always scaled with the resolution `trunc` 
of the model. In case `adjust_Δt_with_output` is true, the `Δt_at_T31` is additionally 
adjusted to the closest divisor of `output_dt` so that the output time axis is keeping
`output_dt` exactly.
"""
function get_Δt_millisec(
    Δt_at_T31::Dates.TimePeriod,
    trunc,
    radius,
    adjust_with_output::Bool,
    output_dt::Dates.TimePeriod = DEFAULT_OUTPUT_DT,
)
    # linearly scale Δt with trunc+1 (which are often powers of two)
    resolution_factor = (DEFAULT_TRUNC+1)/(trunc+1)

    # radius also affects grid spacing, scale proportionally
    radius_factor = radius/DEFAULT_RADIUS

    # maybe rename to _at_trunc_and_radius?
    Δt_at_trunc = Second(Δt_at_T31).value * resolution_factor * radius_factor

    if adjust_with_output && (output_dt > Millisecond(0))
        k = round(Int, Second(output_dt).value / Δt_at_trunc)
        divisors = Primes.divisors(Millisecond(output_dt).value)
        sort!(divisors)
        i = findfirst(x -> x>=k, divisors)
        k_new = isnothing(i) ? k : divisors[i]
        Δt_millisec = Millisecond(round(Int, Millisecond(output_dt).value/k_new))

        # provide info when time step is significantly shortened or lengthened
        Δt_millisec_unadjusted = round(Int, 1000*Δt_at_trunc)
        Δt_ratio = Δt_millisec.value/Δt_millisec_unadjusted

        if abs(Δt_ratio - 1) > 0.05     # only when +-5% changes
            p = round(Int, (Δt_ratio - 1)*100)
            ps = p > 0 ? "+" : ""
            @info "Time step changed from $Δt_millisec_unadjusted to $Δt_millisec ($ps$p%) to match output frequency."
        end
    else 
        Δt_millisec = Millisecond(round(Int, 1000*Δt_at_trunc))
    end

    return Δt_millisec
end 

"""
$(TYPEDSIGNATURES)
Generator function for a Leapfrog struct using `spectral_grid`
for the resolution information."""
function Leapfrog(spectral_grid::SpectralGrid; kwargs...)
    (; NF, trunc, radius) = spectral_grid
    return Leapfrog{NF}(; trunc, radius, kwargs...)
end

"""
$(TYPEDSIGNATURES)
Initialize leapfrogging `L` by recalculating the timestep given the output time step
`output_dt` from `model.output`. Recalculating will slightly adjust the time step to
be a divisor such that an integer number of time steps matches exactly with the output
time step."""
function initialize!(L::Leapfrog, model::AbstractModel)
    (; output_dt) = model.output

    if L.adjust_with_output
        # take actual output dt from model.output and recalculate timestep
        L.Δt_millisec = get_Δt_millisec(L.Δt_at_T31, L.trunc, L.radius, L.adjust_with_output, output_dt)
        L.Δt_sec = L.Δt_millisec.value/1000
        L.Δt = L.Δt_sec/L.radius
    end

    # check how time steps from time integration and output align
    n = round(Int, Millisecond(output_dt).value/L.Δt_millisec.value)
    nΔt = n*L.Δt_millisec
    if nΔt != output_dt
        @warn "$n steps of Δt = $(L.Δt_millisec.value)ms yield output every $(nΔt.value)ms (=$(nΔt.value/1000)s), but output_dt = $(output_dt.value)s"
    end
end

"""$(TYPEDSIGNATURES)
Change time step of timestepper `L` to `Δt` (unscaled)
and disables adjustment to output frequency."""
function set!(
    L::AbstractTimeStepper,
    Δt::Period,                 # unscaled time step in Second, Minute, ...
)
    L.Δt_millisec = Millisecond(Δt)         # recalculate all Δt fields
    L.Δt_sec = L.Δt_millisec.value/1000
    L.Δt = L.Δt_sec/L.radius

    # recalculate the default time step at resolution T31 to be consistent
    resolution_factor = (L.trunc+1)/(DEFAULT_TRUNC+1)
    L.Δt_at_T31 = Second(round(Int, L.Δt_sec*resolution_factor))

    # given Δt was manually set disallow adjustment to output frequency
    L.adjust_with_output = false
    return L
end

# also allow for keyword arguments
set!(L::AbstractTimeStepper; Δt::Period) = set!(L, Δt)

"""
$(TYPEDSIGNATURES)
Performs one leapfrog time step with (`lf=2`) or without (`lf=1`) Robert+Williams filter
(see Williams (2009), Montly Weather Review, Eq. 7-9)."""
function leapfrog!(
    A_old::LowerTriangularArray,        # prognostic variable at t
    A_new::LowerTriangularArray,        # prognostic variable at t+dt
    tendency::LowerTriangularArray,     # tendency (dynamics+physics) of A
    dt::Real,                           # time step (=2Δt, but for init steps =Δt, Δt/2)
    lf::Int,                            # leapfrog index to dis/enable Williams filter
    L::Leapfrog{NF},                    # struct with constants
) where NF                              # number format NF

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
    w1 = lf == 1 ? zero(NF) : robert_filter*williams_filter/2       # = ν*α/2 in Williams (2009, Eq. 8)
    w2 = lf == 1 ? zero(NF) : robert_filter*(1-williams_filter)/2   # = ν(1-α)/2 in Williams (2009, Eq. 9)

    @inbounds for lm in eachindex(A_old, A_new, tendency)
        a_old = A_old[lm]                       # double filtered value from previous time step (t-Δt)
        a_new = a_old + dt_NF*tendency[lm]      # Leapfrog/Euler step depending on dt=Δt, 2Δt (unfiltered at t+Δt)
        a_update = a_old - 2A_lf[lm] + a_new    # Eq. 8&9 in Williams (2009), calculate only once
        A_old[lm] = A_lf[lm] + w1*a_update      # Robert's filter: A_old[lm] becomes 2xfiltered value at t
        A_new[lm] = a_new - w2*a_update         # Williams filter: A_new[lm] becomes 1xfiltered value at t+Δt
    end
end

# variables that are leapfrogged in the respective models, e.g. :vor_tend, :div_tend, etc...
tendency_names(model::AbstractModel) = tuple((Symbol(var, :_tend) for var in prognostic_variables(model))...)

function leapfrog!(
    progn::PrognosticVariables,
    tend::Tendencies,
    dt::Real,               # time step (mostly =2Δt, but for init steps =Δt, Δt/2)
    lf::Int,                # leapfrog index to dis/enable Williams filter
    model::AbstractModel,
)
    for (varname, tendname) in zip(prognostic_variables(model), tendency_names(model))
        var_old, var_new = getfield(progn, varname)
        var_tend = getfield(tend, tendname)
        spectral_truncation!(var_tend)
        leapfrog!(var_old, var_new, var_tend, dt, lf, model.time_stepping)
    end

    # and time stepping for tracers if active
    for (name, tracer) in model.tracers
        if tracer.active
            var_old, var_new = progn.tracers[name]
            var_tend = tend.tracers_tend[name]
            spectral_truncation!(var_tend)
            leapfrog!(var_old, var_new, var_tend, dt, lf, model.time_stepping)
        end
    end

    # evolve the random pattern in time
    random_process!(progn, model.random_process)
    return nothing
end

"""
$(TYPEDSIGNATURES)
Performs the first two initial time steps (Euler forward, unfiltered leapfrog) to populate the
prognostic variables with two time steps (t=0, Δt) that can then be used in the normal leap frogging."""
function first_timesteps!(  
    progn::PrognosticVariables,         # all prognostic variables
    diagn::DiagnosticVariables,         # all pre-allocated diagnostic variables
    model::AbstractModel,                  # everything that is constant at runtime
)
    (; clock) = progn
    clock.n_timesteps == 0 && return nothing    # exit immediately for no time steps
    
    (; implicit) = model
    (; Δt, Δt_millisec) = model.time_stepping
    Δt_millisec_half = Dates.Millisecond(Δt_millisec.value÷2)   # this might be 1ms off

    # FIRST TIME STEP (EULER FORWARD with dt=Δt/2)
    lf1 = 1                             # without Robert+Williams filter
    lf2 = 1                             # evaluates all tendencies at t=0,
                                        # the first leapfrog index (=>Euler forward)
    initialize!(implicit, Δt/2, diagn, model)       # update precomputed implicit terms with time step Δt/2
    timestep!(progn, diagn, Δt/2, model, lf1, lf2)  # update time by half the leapfrog time step Δt used here
    timestep!(clock, Δt_millisec_half, increase_counter=false)      

    # output, callbacks not called after the first Euler step as it's a half-step (i=0 to i=1/2)
    # populating the second leapfrog index to perform the second time step

    # SECOND TIME STEP (UNFILTERED LEAPFROG with dt=Δt, leapfrogging from t=0 over t=Δt/2 to t=Δt)
    initialize!(implicit, Δt, diagn, model)    # update precomputed implicit terms with time step Δt
    lf1 = 1                             # without Robert+Williams filter
    lf2 = 2                             # evaluate all tendencies at t=dt/2,
                                        # the 2nd leapfrog index (=>Leapfrog)
    timestep!(progn, diagn, Δt, model, lf1, lf2)
    # remove prev Δt/2 in case not even milliseconds, otherwise time is off by 1ms
    timestep!(clock, -Δt_millisec_half, increase_counter=false) 
    timestep!(clock, Δt_millisec) 
    
    # do output and callbacks after the first proper (from i=0 to i=1) time step
    output!(model.output, progn, diagn, model)
    callback!(model.callbacks, progn, diagn, model)

    # from now on precomputed implicit terms with 2Δt
    initialize!(implicit, 2Δt, diagn, model) 

    return nothing
end

"""$(TYPEDSIGNATURES)
Calculate a single time step for the barotropic model."""
function timestep!( 
    progn::PrognosticVariables,     # all prognostic variables
    diagn::DiagnosticVariables,     # all pre-allocated diagnostic variables
    dt::Real,                       # time step (mostly =2Δt, but for first_timesteps! =Δt, Δt/2)
    model::Barotropic,              # everything that's constant at runtime
    lf1::Integer = 2,               # leapfrog index 1 (dis/enables Robert+Williams filter)
    lf2::Integer = 2,               # leapfrog index 2 (time step used for tendencies)
)
    model.feedback.nars_detected && return nothing  # exit immediately if NaNs/Infs already present
    
    # set the tendencies back to zero for accumulation
    fill!(diagn.tendencies, 0, Barotropic)

    # TENDENCIES, DIFFUSION, LEAPFROGGING AND TRANSFORM SPECTRAL STATE TO GRID
    dynamics_tendencies!(diagn, progn, lf2, model)
    horizontal_diffusion!(diagn, progn, model.horizontal_diffusion, model)
    leapfrog!(progn, diagn.tendencies, dt, lf1, model)
    transform!(diagn, progn, lf2, model)

    # PARTICLE ADVECTION (always skip 1st step of first_timesteps!)
    not_first_timestep = lf2 == 2
    not_first_timestep && particle_advection!(progn, diagn, model.particle_advection)
end

"""
$(TYPEDSIGNATURES)
Calculate a single time step for the `model <: ShallowWater`."""
function timestep!( 
    progn::PrognosticVariables,     # all prognostic variables
    diagn::DiagnosticVariables,     # all pre-allocated diagnostic variables
    dt::Real,                       # time step (mostly =2Δt, but for first_timesteps! =Δt, Δt/2)
    model::ShallowWater,            # everything that's constant at runtime
    lf1::Integer = 2,               # leapfrog index 1 (dis/enables Robert+Williams filter)
    lf2::Integer = 2,               # leapfrog index 2 (time step used for tendencies)
)
    model.feedback.nars_detected && return nothing  # exit immediately if NaRs already present

    # set the tendencies back to zero for accumulation
    fill!(diagn.tendencies, 0, ShallowWater)

    # GET TENDENCIES, CORRECT THEM FOR SEMI-IMPLICIT INTEGRATION
    dynamics_tendencies!(diagn, progn, lf2, model)
    implicit_correction!(diagn, progn, model.implicit)
    
    # APPLY DIFFUSION, STEP FORWARD IN TIME, AND TRANSFORM NEW TIME STEP TO GRID
    horizontal_diffusion!(diagn, progn, model.horizontal_diffusion, model)
    leapfrog!(progn, diagn.tendencies, dt, lf1, model)
    transform!(diagn, progn, lf2, model)
    
    # PARTICLE ADVECTION (always skip 1st step of first_timesteps!)
    not_first_timestep = lf2 == 2
    not_first_timestep && particle_advection!(progn, diagn, model.particle_advection)
end

"""
$(TYPEDSIGNATURES)
Calculate a single time step for the `model<:PrimitiveEquation`"""
function timestep!( 
    progn::PrognosticVariables,     # all prognostic variables
    diagn::DiagnosticVariables,     # all pre-allocated diagnostic variables
    dt::Real,                       # time step (mostly =2Δt, but for first_timesteps! =Δt, Δt/2)
    model::PrimitiveEquation,       # everything that's constant at runtime
    lf1::Integer = 2,               # leapfrog index 1 (dis/enables Robert+Williams filter)
    lf2::Integer = 2,               # leapfrog index 2 (time step used for tendencies)
)

    model.feedback.nars_detected && return nothing  # exit immediately if NaRs already present
    (; time) = progn.clock                           # current time

    # set the tendencies back to zero for accumulation
    fill!(diagn.tendencies, 0, typeof(model))

    if model.physics                                # switch on/off all physics parameterizations
        # calculate all parameterizations
        parameterization_tendencies!(diagn, progn, time, model)
        
        # time step ocean (temperature and TODO sea ice) and land (temperature and soil moisture)
        # with fluxes from parameterizations
        ocean_timestep!(progn, diagn, model)
        land_timestep!(progn, diagn, model)
        soil_moisture_availability!(diagn, progn, model)
    end

    if model.dynamics                                           # switch on/off all dynamics
        forcing!(diagn, progn, lf2, model)
        drag!(diagn, progn, lf2, model)
        dynamics_tendencies!(diagn, progn, lf2, model)          # dynamical core
        implicit_correction!(diagn, model.implicit, progn)      # semi-implicit time stepping corrections
    else    # just transform physics tendencies to spectral space
        physics_tendencies_only!(diagn, model)
    end

    # APPLY DIFFUSION, STEP FORWARD IN TIME, AND TRANSFORM NEW TIME STEP TO GRID
    horizontal_diffusion!(diagn, progn, model.horizontal_diffusion, model)
    leapfrog!(progn, diagn.tendencies, dt, lf1, model)
    transform!(diagn, progn, lf2, model)

    # PARTICLE ADVECTION (always skip 1st step of first_timesteps!)
    not_first_timestep = lf2 == 2
    not_first_timestep && particle_advection!(progn, diagn, model.particle_advection)
end

"""
$(TYPEDSIGNATURES)
Main time loop that that initializes output and feedback, loops over all time steps
and calls the output and feedback functions."""
function time_stepping!(
    progn::PrognosticVariables,         # all prognostic variables
    diagn::DiagnosticVariables,         # all pre-allocated diagnostic variables
    model::AbstractModel,               # all model components
)          
    
    (; clock) = progn
    (; Δt, Δt_millisec) = model.time_stepping

    # SCALING: we use vorticity*radius, divergence*radius in the dynamical core
    scale!(progn, diagn, model.spectral_grid.radius)

    # OUTPUT INITIALISATION AND STORING INITIAL CONDITIONS + FEEDBACK
    # propagate spectral state to grid variables for initial condition output
    (; output, feedback) = model
    lf = 1                                  # use first leapfrog index
    transform!(diagn, progn, lf, model, initialize=true)
    initialize!(progn.particles, progn, diagn, model.particle_advection)
    initialize!(output, feedback, progn, diagn, model)
    initialize!(model.callbacks, progn, diagn, model)
    
    # FIRST TIMESTEPS: EULER FORWARD THEN 1x LEAPFROG
    # considered part of the model initialisation
    first_timesteps!(progn, diagn, model)
    
    # only now initialise feedback for benchmark accuracy
    initialize!(feedback, clock, model)

    # MAIN LOOP
    for _ in 2:clock.n_timesteps            # start at 2 as first Δt in first_timesteps!
        timestep!(progn, diagn, 2Δt, model) # calculate tendencies and leapfrog forward
        timestep!(clock, Δt_millisec)       # time of lf=2 and diagn after timestep!

        progress!(feedback, progn)          # updates the progress meter bar
        output!(output, progn, diagn, model)
        callback!(model.callbacks, progn, diagn, model)
    end
    
    # UNSCALE, CLOSE, FINALIZE
    finalize!(feedback)                     # finish the progress meter, do first for benchmark accuracy
    unscale!(progn)                         # undo radius-scaling for vor, div from the dynamical core
    unscale!(diagn)                         # undo radius-scaling for vor, div from the dynamical core
    finalize!(output, progn, diagn, model)  # possibly post-process output, then close netCDF file
    write_restart_file(output, progn)       # as JLD2 
    finalize!(model.callbacks, progn, diagn, model)

    # return a UnicodePlot of surface vorticity
    surface_vorticity = diagn.grid.vor_grid[:, end]
    return plot(surface_vorticity, title="Surface relative vorticity [1/s]")
end 