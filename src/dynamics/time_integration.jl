const DEFAULT_NSTEPS = 2
export Leapfrog

# ============================================================================
# LEAPFROG TIME STEPPING
# ============================================================================

"""Leapfrog time stepping defined by the following fields
$(TYPEDFIELDS)"""
@kwdef mutable struct Leapfrog{NF<:AbstractFloat} <: AbstractTimeStepper
    "[DERIVED] Spectral resolution (max degree of spherical harmonics)"
    trunc::Int                      

    "[CONST] Number of time steps stored simultaneously in prognostic variables"
    nsteps::Int = 2

    "[OPTION] Time step in minutes for T31, scale linearly to `trunc`"
    Δt_at_T31::Second = Minute(40)

    "[OPTION] Adjust `Δt_at_T31` with the `output_dt` to reach `output_dt` exactly in integer time steps"
    adjust_with_output::Bool = true

    "[OPTION] Start integration with (1) Euler step with dt/2, (2) Leapfrog step with dt"
    start_with_euler::Bool = true

    "[OPTION] Sets `first_step_euler=false` after first step to continue with leapfrog after 1st `run!` call"
    continue_with_leapfrog::Bool = true

    "[DERIVED] Use Euler on first time step? (controlled by `start_with_euler` and `continue_with_leapfrog`)"
    first_step_euler::Bool = start_with_euler

    "[OPTION] Robert (1966) time filter coefficient to suppress the computational mode"
    robert_filter::NF = 0.1

    "[OPTION] Williams time filter (Amezcua 2011) coefficient for 3rd order acc"
    williams_filter::NF = 0.53

    "[DERIVED] Radius of sphere [m], used for scaling, set in `initialize!` to `planet.radius`"
    radius::NF = DEFAULT_RADIUS

    "[DERIVED] Time step Δt in milliseconds at specified resolution"
    Δt_millisec::Millisecond = get_Δt_millisec(Second(Δt_at_T31), trunc, radius, adjust_with_output)

    "[DERIVED] Time step Δt [s] at specified resolution"
    Δt_sec::NF = Δt_millisec.value/1000

    "[DERIVED] Time step Δt [s/m] at specified resolution, scaled by 1/radius"
    Δt::NF = Δt_sec/radius  
end

"""$(TYPEDSIGNATURES)
Computes the time step in [ms]. `Δt_at_T31` is always scaled with the resolution `trunc` 
of the model. In case `adjust_Δt_with_output` is true, the `Δt_at_T31` is additionally 
adjusted to the closest divisor of `output_dt` so that the output time axis is keeping
`output_dt` exactly."""
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

        if abs(Δt_ratio - 1) > 0.05     # print info only when +-5% changes
            p = round(Int, (Δt_ratio - 1)*100)
            ps = p > 0 ? "+" : ""
            @info "Time step changed from $Δt_millisec_unadjusted to $Δt_millisec ($ps$p%) to match output frequency."
        end
    else 
        Δt_millisec = Millisecond(round(Int, 1000*Δt_at_trunc))
    end

    return Δt_millisec
end 

"""$(TYPEDSIGNATURES)
Generator function for a Leapfrog struct using `spectral_grid`
for the resolution information."""
function Leapfrog(spectral_grid::SpectralGrid; kwargs...)
    (; NF, trunc) = spectral_grid
    return Leapfrog{NF}(; trunc, kwargs...)
end

"""$(TYPEDSIGNATURES)
Initialize leapfrogging `L` by recalculating the time step given the output time step
`output_dt` from `model.output`. Recalculating will slightly adjust the time step to
be a divisor such that an integer number of time steps matches exactly with the output
time step."""
function initialize!(L::Leapfrog, model::AbstractModel)
    (; output_dt) = model.output
    (; radius) = model.planet

    # take radius from planet and recalculate time step and possibly adjust with output dt
    L.Δt_millisec = get_Δt_millisec(L.Δt_at_T31, L.trunc, radius, L.adjust_with_output, output_dt)
    L.Δt_sec = L.Δt_millisec.value/1000
    L.Δt = L.Δt_sec/radius

    # check how time steps from time integration and output align
    n = round(Int, Millisecond(output_dt).value/L.Δt_millisec.value)
    nΔt = n*L.Δt_millisec
    if nΔt != output_dt
        @warn "$n steps of Δt = $(L.Δt_millisec.value)ms yield output every $(nΔt.value)ms (=$(nΔt.value/1000)s), but output_dt = $(output_dt.value)s"
    end
    if L.start_with_euler
        L.first_step_euler = true
    end

    return nothing
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

    launch!(architecture(tendency), SpectralWorkOrder, size(tendency), leapfrog_kernel!, A_old, A_new, A_lf, tendency, dt_NF, w1, w2)

    return nothing
end

@kernel inbounds=true function leapfrog_kernel!(A_old, A_new, A_lf, tendency, @Const(dt), @Const(w1), @Const(w2))

    lmk = @index(Global, Linear)    # every harmonic lm, every vertical layer k

    a_old = A_old[lmk]
    a_new = a_old + dt*tendency[lmk]
    a_update = a_old - 2A_lf[lmk] + a_new
    A_old[lmk] = A_lf[lmk] + w1*a_update
    A_new[lmk] = a_new - w2*a_update
end

# variables that are leapfrogged in the respective models, e.g. :vor_tend, :div_tend, etc...
tendency_names(model::AbstractModel) = tuple((Symbol(var, :_tend) for var in prognostic_variables(model))...)

"""$(TYPEDSIGNATURES)
Leapfrog time stepping for all prognostic variables in `progn` using their tendencies in `tend`.
Depending on `model` decides which variables to time step."""
function leapfrog!(
    progn::PrognosticVariables,
    tend::Tendencies,
    dt::Real,               # time step (mostly =2Δt, but for init steps =Δt, Δt/2)
    lf::Int,                # leapfrog index to dis/enable Williams filter
    model::AbstractModel,
)
    for (varname, tendname) in zip(prognostic_variables(model), tendency_names(model))
        var = getfield(progn, varname)
        var_old, var_new = get_steps(var)
        var_tend = getfield(tend, tendname)
        spectral_truncation!(var_tend)
        leapfrog!(var_old, var_new, var_tend, dt, lf, model.time_stepping)
    end

    # and time stepping for tracers if active
    for (name, tracer) in model.tracers
        if tracer.active
            var_old, var_new = get_steps(progn.tracers[name])
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
    model::AbstractModel,               # everything that is constant at runtime
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
    callback!(model.callbacks, progn, diagn, model)
    output!(model.output, Simulation(progn, diagn, model))

    # from now on precomputed implicit terms with 2Δt
    initialize!(implicit, 2Δt, diagn, model)
    set_initialized!(implicit)      # mark implicit as initialized

    return nothing
end

# ============================================================================
# LORENZ N-CYCLE TIME STEPPING
# ============================================================================

export NCycleLorenz, NCycleLorenzA, NCycleLorenzB, NCycleLorenzAB, NCycleLorenzABBA

"""Abstract type for Lorenz N-cycle variants following Hotta et al. (2016)."""
abstract type NCycleLorenzVariant end

"""Version A: weights w_k = N/(N-k) for k=1,...,N-1; w_0 = 1"""
struct NCycleLorenzA <: NCycleLorenzVariant end

"""Version B: weights w_k = N/k for k=1,...,N-1; w_0 = 1"""
struct NCycleLorenzB <: NCycleLorenzVariant end

"""Version AB: alternates A and B every N steps"""
struct NCycleLorenzAB <: NCycleLorenzVariant end

"""Version ABBA: uses A-B-B-A sequence (only for N=4, provides 4th order accuracy)"""
struct NCycleLorenzABBA <: NCycleLorenzVariant end

"""
    NCycleLorenz{NF, V} <: AbstractTimeStepper

A semi-implicit Lorenz N-cycle time integration scheme following Hotta et al. (2016).

# Algorithm (per substep)
1. G = w*F_E(x) + (1-w)*G     (weighted tendency accumulation)
2. dx = (I - α*Δt*L_I)^(-1) * (G + L_I*x)  (implicit solve)
3. x = x + Δt*dx              (state update)
$(TYPEDFIELDS)
"""
@kwdef mutable struct NCycleLorenz{NF<:AbstractFloat, V<:NCycleLorenzVariant} <: AbstractTimeStepper
    "[DERIVED] Spectral resolution (max degree of spherical harmonics)"
    trunc::Int

    "[CONST] Number of time steps stored (2: current state + tendency accumulator)"
    nsteps::Int = 2
    
    "[OPTION] Number of cycles N (3 or 4 recommended, 4 is more stable)"
    cycles::Int = 4
    
    "[OPTION] Variant: NCycleLorenzA() (default), B, AB, or ABBA"
    variant::V = NCycleLorenzA()
    
    "[OPTION] Centering parameter: 0.5=Crank-Nicolson (2nd order), 1.0=Backward Euler (1st order)"
    α::NF = 0.5
    
    "[OPTION] Time step in minutes for T31, scale linearly to `trunc`"
    Δt_at_T31::Second = Minute(30)

    "[OPTION] Adjust `Δt_at_T31` with the `output_dt` to reach `output_dt` exactly"
    adjust_with_output::Bool = true

    "[OPTION] No Euler first step needed"
    first_step_euler::Bool = false           # Not used for Lorenz, but needed for compatibility

    "[DERIVED] Radius of sphere [m], set in `initialize!` to `planet.radius`"
    radius::NF = DEFAULT_RADIUS

    "[DERIVED] Time step Δt in milliseconds at specified resolution"
    Δt_millisec::Millisecond = get_Δt_millisec(Second(Δt_at_T31), trunc, radius, adjust_with_output)

    "[DERIVED] Time step Δt [s] at specified resolution"
    Δt_sec::NF = Δt_millisec.value/1000

    "[DERIVED] Time step Δt [s/m] at specified resolution, scaled by 1/radius"
    Δt::NF = Δt_sec/radius
end

"""$(TYPEDSIGNATURES)
Generator function for NCycleLorenz struct using `spectral_grid` for resolution."""
function NCycleLorenz(spectral_grid::SpectralGrid; kwargs...)
    (; NF, trunc) = spectral_grid
    return NCycleLorenz{NF, NCycleLorenzA}(; trunc, kwargs...)
end

"""$(TYPEDSIGNATURES)
Get current substep within N-cycle (0 to N-1) from the clock."""
current_substep(L::NCycleLorenz, clock) = mod(clock.timestep_counter, L.cycles)

"""$(TYPEDSIGNATURES)
Initialize NCycleLorenz time stepper."""
function initialize!(L::NCycleLorenz, model::AbstractModel)
    (; output_dt) = model.output
    (; radius) = model.planet

    # Validate compatibility - runs ONCE, not every timestep
    if L.variant isa NCycleLorenzABBA && L.cycles != 4
        @warn "ABBA variant designed for N=4 (4th order accurate), but N=$(L.cycles). Consider cycles=4 or variant A/B/AB."
    end

    L.Δt_millisec = get_Δt_millisec(L.Δt_at_T31, L.trunc, radius, L.adjust_with_output, output_dt)
    L.Δt_sec = L.Δt_millisec.value/1000
    L.Δt = L.Δt_sec/radius

    n = round(Int, Millisecond(output_dt).value/L.Δt_millisec.value)
    nΔt = n*L.Δt_millisec
    if nΔt != output_dt
        @warn "$n steps of Δt = $(L.Δt_millisec.value)ms yield output every $(nΔt.value)ms, but output_dt = $(output_dt.value)ms"
    end

    return nothing
end

# Weight coefficient functions (Hotta et al. 2016, Eqs 2-3, 7-8)
"""$(TYPEDSIGNATURES)
Compute weight coefficient w for current substep."""
function weight_coefficient(L::NCycleLorenz{NF}, clock) where NF
    k = current_substep(L, clock)
    return weight_coefficient(NF, L.variant, k, L.cycles)
end

# Type-stable versions with explicit NF
@inline function weight_coefficient(::Type{NF}, ::NCycleLorenzA, k::Int, N::Int) where NF
    return k == 0 ? one(NF) : convert(NF, N) / convert(NF, N - k)
end

@inline function weight_coefficient(::Type{NF}, ::NCycleLorenzB, k::Int, N::Int) where NF
    return k == 0 ? one(NF) : convert(NF, N) / convert(NF, k)
end

# Fallback for when NF is not specified (uses DEFAULT_NF)
@inline weight_coefficient(V::NCycleLorenzVariant, k::Int, N::Int) = 
    weight_coefficient(DEFAULT_NF, V, k, N)


# ============================================================================
# LORENZ N-CYCLE HELPER KERNELS
# ============================================================================

"""$(TYPEDSIGNATURES)
Accumulate weighted tendency: G = w*F_E + (1-w)*G (Hotta et al. 2016, Eq. 19)"""
function accumulate_tendency!(
    G::LowerTriangularArray,
    F_explicit::LowerTriangularArray,
    w::Real,
    L::NCycleLorenz{NF},
) where NF
    @boundscheck size(G) == size(F_explicit) || throw(BoundsError())
    
    w_NF = convert(NF, w)
    
    launch!(
        architecture(G), 
        SpectralWorkOrder, 
        size(G), 
        accumulate_tendency_kernel!, 
        G, F_explicit, w_NF
    )
    
    return nothing
end

@kernel inbounds=true function accumulate_tendency_kernel!(
    G, F_explicit, @Const(w)
)
    lmk = @index(Global, Linear)
    G[lmk] = w * F_explicit[lmk] + (1 - w) * G[lmk]
end

"""$(TYPEDSIGNATURES)
Update state: x = x + Δt*G (Hotta et al. 2016, Eq. 21)"""
function update_state!(
    x::LowerTriangularArray,
    G::LowerTriangularArray,
    dt::Real,
    L::NCycleLorenz{NF},
) where NF
    @boundscheck size(x) == size(G) || throw(BoundsError())
    
    dt_NF = convert(NF, dt)
    
    launch!(
        architecture(x), 
        SpectralWorkOrder, 
        size(x), 
        update_state_kernel!, 
        x, G, dt_NF
    )
    
    return nothing
end

@kernel inbounds=true function update_state_kernel!(
    x, G, @Const(dt)
)
    lmk = @index(Global, Linear)
    x[lmk] = x[lmk] + dt * G[lmk]
end

"""$(TYPEDSIGNATURES)
Lorenz N-cycle time stepping for all prognostic variables following Hotta et al. (2016).

This function performs ONE substep of the N-cycle. The clock tracks which substep we're on.

Algorithm per substep:
1. Accumulate weighted tendency: G = w*F_E + (1-w)*G (Eq. 19)
2. Apply implicit correction: G → (I - α*Δt*L_I)^(-1) * (G + L_I*x) (Eq. 20)
3. Update state: x = x + Δt*G (Eq. 21)
"""
function lorenz_ncycle_timestep!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    tend::Tendencies,
    dt::Real,
    model::AbstractModel,
)
    L = model.time_stepping::NCycleLorenz
    clock = progn.clock
    
    # Get current weight from the clock
    w_current = weight_coefficient(L, clock)
    
    # Step 1: Accumulate weighted explicit tendencies (Eq. 19)
    for (varname, tendname) in zip(prognostic_variables(model), tendency_names(model))
        var = getfield(progn, varname)
        G = get_step(var, 2)
        
        F_explicit = getfield(tend, tendname)
        spectral_truncation!(F_explicit)
        
        accumulate_tendency!(G, F_explicit, w_current, L)
    end
    
    # Accumulate for tracers
    for (name, tracer) in model.tracers
        if tracer.active
            var = progn.tracers[name]
            G = get_step(var, 2)
            
            F_explicit = tend.tracers_tend[name]
            spectral_truncation!(F_explicit)
            
            accumulate_tendency!(G, F_explicit, w_current, L)
        end
    end
    
    # Step 2: Apply implicit correction (Eq. 20)
    lorenz_implicit_correction!(diagn, progn, model.implicit, model)
    
    # Step 3: Update states with corrected tendencies (Eq. 21)
    for varname in prognostic_variables(model)
        var = getfield(progn, varname)
        x = get_step(var, 1)
        G = get_step(var, 2)
        
        update_state!(x, G, dt, L)
    end
    
    # Update tracers
    for (name, tracer) in model.tracers
        if tracer.active
            var = progn.tracers[name]
            x = get_step(var, 1)
            G = get_step(var, 2)
            
            update_state!(x, G, dt, L)
        end
    end
    
    
    # Evolve random pattern
    random_process!(progn, model.random_process)
    
    return nothing
end



# ============================================================================
# MODEL-SPECIFIC TIMESTEP FUNCTIONS
# These compute tendencies and call the appropriate time stepping scheme
# The original timestep! functions remain unchanged for Leapfrog
# We add new methods that check model.time_stepping type at runtime
# ============================================================================

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
    model.feedback.nans_detected && return nothing  # exit immediately if NaNs/Infs already present

    # set the tendencies back to zero for accumulation
    fill!(diagn.tendencies, 0, Barotropic)

    # Check which time stepper we're using and dispatch accordingly
    if model.time_stepping isa NCycleLorenz
        # Lorenz N-cycle: always use index 1 for current state
        dynamics_tendencies!(diagn, progn, 1, model)
        horizontal_diffusion!(diagn, progn, model.horizontal_diffusion, model)
        lorenz_ncycle_timestep!(progn, diagn, diagn.tendencies, dt, model)
        transform!(diagn, progn, 1, model)
        
        # Particle advection (skip first step)
        progn.clock.timestep_counter > 0 && particle_advection!(progn, diagn, model)
    else
        # Leapfrog (original code path)
        dynamics_tendencies!(diagn, progn, lf2, model)
        horizontal_diffusion!(diagn, progn, model.horizontal_diffusion, model)
        leapfrog!(progn, diagn.tendencies, dt, lf1, model)
        transform!(diagn, progn, lf2, model)
        
        # PARTICLE ADVECTION (always skip 1st step of first_timesteps!)
        not_first_timestep = lf2 == 2
        not_first_timestep && particle_advection!(progn, diagn, model)
    end

    return nothing 
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
    model.feedback.nans_detected && return nothing  # exit immediately if NaNs already present

    # set the tendencies back to zero for accumulation
    fill!(diagn.tendencies, 0, ShallowWater)

    if model.time_stepping isa NCycleLorenz
        # Lorenz N-cycle path
        dynamics_tendencies!(diagn, progn, 1, model)
        horizontal_diffusion!(diagn, progn, model.horizontal_diffusion, model)
        # Note: implicit correction is called inside lorenz_ncycle_timestep!
        lorenz_ncycle_timestep!(progn, diagn, diagn.tendencies, dt, model)
        transform!(diagn, progn, 1, model)
        
        progn.clock.timestep_counter > 0 && particle_advection!(progn, diagn, model)
    else
        # Leapfrog path (original)
        dynamics_tendencies!(diagn, progn, lf2, model)
        implicit_correction!(diagn, progn, model.implicit, model)
        horizontal_diffusion!(diagn, progn, model.horizontal_diffusion, model)
        leapfrog!(progn, diagn.tendencies, dt, lf1, model)
        transform!(diagn, progn, lf2, model)
        
        not_first_timestep = lf2 == 2
        not_first_timestep && particle_advection!(progn, diagn, model)
    end

    return nothing
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

    model.feedback.nans_detected && return nothing  # exit immediately if NaNs already present
    (; time) = progn.clock                           # current time

    # set the tendencies back to zero for accumulation
    fill!(diagn.tendencies, 0, typeof(model))

    if model.physics                                # switch on/off all physics parameterizations
        # calculate all parameterizations
        parameterization_tendencies!(diagn, progn, time, model)
        ocean_timestep!(progn, diagn, model)    # sea surface temperature and maybe in the future sea ice
        sea_ice_timestep!(progn, diagn, model)  # sea ice
        land_timestep!(progn, diagn, model)     # soil moisture and temperature, vegetation, maybe rivers
    end

    if model.time_stepping isa NCycleLorenz
        # Lorenz N-cycle path
        if model.dynamics
            forcing!(diagn, progn, 1, model)
            drag!(diagn, progn, 1, model)
            dynamics_tendencies!(diagn, progn, 1, model)
            # Note: implicit correction is called inside lorenz_ncycle_timestep!
        else
            physics_tendencies_only!(diagn, model)
        end
        
        horizontal_diffusion!(diagn, progn, model.horizontal_diffusion, model)
        lorenz_ncycle_timestep!(progn, diagn, diagn.tendencies, dt, model)
        transform!(diagn, progn, 1, model)
        
        progn.clock.timestep_counter > 0 && particle_advection!(progn, diagn, model)
    else
        # Leapfrog path 
        if model.dynamics
            forcing!(diagn, progn, lf2, model)
            drag!(diagn, progn, lf2, model)
            dynamics_tendencies!(diagn, progn, lf2, model)
            implicit_correction!(diagn, progn, model.implicit, model)
        else
            physics_tendencies_only!(diagn, model)
        end
        
        horizontal_diffusion!(diagn, progn, model.horizontal_diffusion, model)
        leapfrog!(progn, diagn.tendencies, dt, lf1, model)
        transform!(diagn, progn, lf2, model)
        
        not_first_timestep = lf2 == 2
        not_first_timestep && particle_advection!(progn, diagn, model)
    end

    return nothing 
end

# ============================================================================
# MAIN TIME STEPPING LOOP
# ============================================================================

"""$(TYPEDSIGNATURES)
First 1 or 2 time steps of `simulation`. If `model.time_stepping.start_with_euler` is true,
then start with one Euler step with dt/2, followed by one Leapfrog step with dt.
If false, continue with leapfrog steps at 2Δt (e.g. restart)."""
function first_timesteps!(simulation::AbstractSimulation)
    progn, diagn, model = unpack(simulation)
    (; time_stepping) = model
    
    # Check which time stepper we're using
    if time_stepping isa NCycleLorenz
        # Lorenz N-cycle: self-starting, no Euler initialization needed
        (; Δt, α) = time_stepping
        
        # Ensure implicit solver uses the same α as the time stepper
        if !isnothing(model.implicit) && hasfield(typeof(model.implicit), :α)
            model.implicit.α = α
        end
        
        # Initialize implicit solver with Δt (not 2Δt like leapfrog)
        initialize!(model.implicit, Δt, diagn, model)
        set_initialized!(model.implicit)
        
        # Initialize feedback
        initialize!(model.feedback, progn.clock, model)
        
        # Just do a normal timestep since Lorenz N-cycle is self-starting!
        later_timestep!(simulation)
    else
        # Leapfrog (original code path)
        (; Δt) = time_stepping
        
        # decide whether to start with 1x Euler then 1x Leapfrog at Δt
        if time_stepping.first_step_euler
            first_timesteps!(progn, diagn, model)
            time_stepping.first_step_euler = !time_stepping.continue_with_leapfrog   # after first run! continue with leapfrog
            
        else    # or continue with leaprog steps at 2Δt (e.g. restart)
                # but make sure that implicit solver is initialized in that situation
            initialize!(model.implicit, 2Δt, diagn, model)
            set_initialized!(model.implicit)            # mark implicit as initialized
            later_timestep!(simulation)
        end

        # only now initialise feedback for benchmark accuracy
        (; clock) = progn
        initialize!(model.feedback, clock, model)
    end
    
    return nothing
end

"""$(TYPEDSIGNATURES)
Perform one single time step of `simulation` including
possibly output and callbacks."""
function timestep!(simulation::AbstractSimulation)
    (; clock) = simulation.prognostic_variables

    if clock.timestep_counter == 0
        first_timesteps!(simulation)
    else
        later_timestep!(simulation)
    end
end

"""$(TYPEDSIGNATURES)
Perform one single "normal" time step of `simulation`, after `first_timesteps!`."""
function later_timestep!(simulation::AbstractSimulation)
    progn, diagn, model = unpack(simulation)
    (; feedback, output) = model
    (; time_stepping) = model
    (; clock) = progn

    if time_stepping isa NCycleLorenz
        # Lorenz N-cycle uses Δt per substep
        (; Δt, Δt_millisec) = time_stepping
        timestep!(progn, diagn, Δt, model)
        timestep!(clock, Δt_millisec)
    else
        # Leapfrog uses 2Δt
        (; Δt, Δt_millisec) = time_stepping
        timestep!(progn, diagn, 2Δt, model)
        timestep!(clock, Δt_millisec)
    end

    progress!(feedback, progn)                      # updates the progress meter bar
    output!(output, simulation)                     # do output?
    callback!(model.callbacks, progn, diagn, model) # any callbacks?
    return nothing
end

"""$(TYPEDSIGNATURES)
Main time loop that loops over all time steps."""
function time_stepping!(simulation::AbstractSimulation)          
    (; clock) = simulation.prognostic_variables
    for _ in 1:clock.n_timesteps        # MAIN LOOP
        timestep!(simulation)
    end
end