# ============================================================================
# LORENZ N-CYCLE TIME STEPPING
# Based on Hotta et al. (2016), Monthly Weather Review
# ============================================================================

export NCycleLorenz, NCycleLorenzA, NCycleLorenzB, NCycleLorenzAB, NCycleLorenzABBA

# ============================================================================
# TYPE DEFINITIONS
# ============================================================================

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

# Key differences from Leapfrog
- Self-starting (no special first steps needed)
- Two time levels: current state (index 1) + tendency accumulator (index 2)
- No computational mode (no filtering needed)
- Stable for dissipative systems

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

    "[DERIVED] Radius of sphere [m], set in `initialize!` to `planet.radius`"
    radius::NF = DEFAULT_RADIUS

    "[DERIVED] Time step Δt in milliseconds at specified resolution"
    Δt_millisec::Millisecond = get_Δt_millisec(Second(Δt_at_T31), trunc, radius, adjust_with_output)

    "[DERIVED] Time step Δt [s] at specified resolution"
    Δt_sec::NF = Δt_millisec.value/1000

    "[DERIVED] Time step Δt [s/m] at specified resolution, scaled by 1/radius"
    Δt::NF = Δt_sec/radius
end

# ============================================================================
# CONSTRUCTOR AND INITIALIZATION
# ============================================================================

"""$(TYPEDSIGNATURES)
Generator function for NCycleLorenz struct using `spectral_grid` for resolution."""
function NCycleLorenz(spectral_grid::SpectralGrid; kwargs...)
    (; NF, trunc) = spectral_grid
    return NCycleLorenz{NF, NCycleLorenzA}(; trunc, kwargs...)
end

"""$(TYPEDSIGNATURES)
Initialize NCycleLorenz time stepper.

Unlike Leapfrog, this is straightforward:
- Calculate final time step based on resolution and output frequency
- Initialize implicit solver with Δt (not 2Δt)
- No special first-step handling needed
"""
function initialize!(L::NCycleLorenz, model::AbstractModel)
    (; output_dt) = model.output
    (; radius) = model.planet

    # Validate compatibility - runs ONCE, not every timestep
    if L.variant isa NCycleLorenzABBA && L.cycles != 4
        @warn "ABBA variant designed for N=4 (4th order accurate), but N=$(L.cycles). Consider cycles=4 or variant A/B/AB."
    end

    # Calculate time step (may be adjusted to match output_dt)
    L.Δt_millisec = get_Δt_millisec(L.Δt_at_T31, L.trunc, radius, L.adjust_with_output, output_dt)
    L.Δt_sec = L.Δt_millisec.value/1000
    L.Δt = L.Δt_sec/radius

    # Verify alignment with output frequency
    n = round(Int, Millisecond(output_dt).value/L.Δt_millisec.value)
    nΔt = n*L.Δt_millisec
    if nΔt != output_dt
        @warn "$n steps of Δt = $(L.Δt_millisec.value)ms yield output every $(nΔt.value)ms, but output_dt = $(output_dt.value)ms"
    end

    return nothing
end

# ============================================================================
# WEIGHT COEFFICIENTS (Hotta et al. 2016, Eqs 2-3, 7-8)
# ============================================================================

"""$(TYPEDSIGNATURES)
Get current substep within N-cycle (0 to N-1) from the clock."""
current_substep(L::NCycleLorenz, clock) = mod(clock.timestep_counter, L.cycles)

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

# Fallback for when NF is not specified
@inline weight_coefficient(V::NCycleLorenzVariant, k::Int, N::Int) = 
    weight_coefficient(DEFAULT_NF, V, k, N)

# ============================================================================
# HELPER KERNELS FOR STATE UPDATES
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

# ============================================================================
# CORE LORENZ N-CYCLE TIMESTEP ALGORITHM
# ============================================================================

"""$(TYPEDSIGNATURES)
Lorenz N-cycle time stepping for all prognostic variables following Hotta et al. (2016).

This function performs ONE substep of the N-cycle. The clock tracks which substep we're on.

Algorithm per substep:
1. Accumulate weighted tendency: G = w*F_E + (1-w)*G (Eq. 19)
2. Apply implicit correction: G → (I - α*Δt*L_I)^(-1) * (G + L_I*x) (Eq. 20)
3. Update state: x = x + Δt*G (Eq. 21)

Notes:
- G (tendency accumulator) is stored in progn[:,:,2]
- x (current state) is stored in progn[:,:,1]
- No separate leapfrog indices or special first-step handling
"""
function lorenz_ncycle_step!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    tend::Tendencies,
    dt::Real,
    model::AbstractModel,
)
    L = model.time_stepping::NCycleLorenz
    clock = progn.clock
    
    # RESET tendency accumulator at start of each N-cycle
    if current_substep(L, clock) == 0
        for varname in prognostic_variables(model)
            var = getfield(progn, varname)
            G = get_step(var, 2)
            fill!(G, 0)
        end

        for (name, tracer) in model.tracers
            tracer.active || continue
            G = get_step(progn.tracers[name], 2)
            fill!(G, 0)
        end
    end

    # ----------------------------------------------------------------------
    # STEP 1: Accumulate weighted explicit tendencies (Hotta et al. Eq. 19)
    # G ← w F_E + (1 − w) G
    # ----------------------------------------------------------------------
    w_current = weight_coefficient(L, clock)

    for (varname, tendname) in zip(prognostic_variables(model), tendency_names(model))
        var = getfield(progn, varname)
        G = get_step(var, 2)

        F_explicit = getfield(tend, tendname)
        spectral_truncation!(F_explicit)

        accumulate_tendency!(G, F_explicit, w_current, L)
    end

    for (name, tracer) in model.tracers
        tracer.active || continue
        var = progn.tracers[name]
        G = get_step(var, 2)

        F_explicit = tend.tracers_tend[name]
        spectral_truncation!(F_explicit)

        accumulate_tendency!(G, F_explicit, w_current, L)
    end

    # ----------------------------------------------------------------------
    # STEP 2: Semi-implicit correction (Hotta et al. Eq. 20)
    # G ← (I − ξ L_I)⁻¹ (G + L_I x)
    # ----------------------------------------------------------------------
    lorenz_implicit_correction!(diagn, progn, model.implicit, model)

    # ----------------------------------------------------------------------
    # STEP 3: State update (Hotta et al. Eq. 21)
    # x ← x + Δt G
    # ----------------------------------------------------------------------
    for varname in prognostic_variables(model)
        var = getfield(progn, varname)
        x = get_step(var, 1)
        G = get_step(var, 2)

        update_state!(x, G, dt, L)
    end

    for (name, tracer) in model.tracers
        tracer.active || continue
        var = progn.tracers[name]
        x = get_step(var, 1)
        G = get_step(var, 2)

        update_state!(x, G, dt, L)
    end

    random_process!(progn, model.random_process)
    return nothing
end


# ============================================================================
# MODEL-SPECIFIC TIMESTEP FUNCTIONS
# These compute tendencies and call the Lorenz N-cycle timestepper
# ============================================================================

# Helper to get tendency names
tendency_names(model::AbstractModel) = tuple((Symbol(var, :_tend) for var in prognostic_variables(model))...)

"""$(TYPEDSIGNATURES)
Calculate a single time step for the barotropic model."""
function timestep!( 
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    dt::Real,
    model::Barotropic,
)
    model.feedback.nans_detected && return nothing

    # Reset tendencies to zero
    fill!(diagn.tendencies, 0, Barotropic)

    # Compute dynamics tendencies (always use index 1 for current state)
    dynamics_tendencies!(diagn, progn, 1, model)
    horizontal_diffusion!(diagn, progn, model.horizontal_diffusion, model)
    
    # Perform one substep of the N-cycle
    lorenz_ncycle_step!(progn, diagn, diagn.tendencies, dt, model)
    
    # Transform to grid space
    transform!(diagn, progn, 1, model)
    
    # Particle advection (skip first step)
    progn.clock.timestep_counter > 0 && particle_advection!(progn, diagn, model)

    return nothing 
end

"""$(TYPEDSIGNATURES)
Calculate a single time step for the shallow water model."""
function timestep!( 
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    dt::Real,
    model::ShallowWater,
)
    model.feedback.nans_detected && return nothing

    # Reset tendencies to zero
    fill!(diagn.tendencies, 0, ShallowWater)

    # Compute dynamics tendencies
    dynamics_tendencies!(diagn, progn, 1, model)
    horizontal_diffusion!(diagn, progn, model.horizontal_diffusion, model)
    
    # Perform one substep of the N-cycle (implicit correction is inside)
    lorenz_ncycle_step!(progn, diagn, diagn.tendencies, dt, model)
    
    # Transform to grid space
    transform!(diagn, progn, 1, model)
    
    # Particle advection (skip first step)
    progn.clock.timestep_counter > 0 && particle_advection!(progn, diagn, model)

    return nothing
end

"""$(TYPEDSIGNATURES)
Calculate a single time step for the primitive equation model."""
function timestep!( 
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    dt::Real,
    model::PrimitiveEquation,
)
    model.feedback.nans_detected && return nothing

    # Reset tendencies to zero
    fill!(diagn.tendencies, 0, typeof(model))

    # Physics parameterizations
    if model.physics
        parameterization_tendencies!(diagn, progn, model)  # ← REMOVE 'time' parameter
        ocean_timestep!(progn, diagn, model)
        sea_ice_timestep!(progn, diagn, model)
        land_timestep!(progn, diagn, model)
    end

    # Dynamics
    if model.dynamics
        forcing!(diagn, progn, 1, model)
        drag!(diagn, progn, 1, model)
        dynamics_tendencies!(diagn, progn, 1, model)
    else
        physics_tendencies_only!(diagn, model)
    end
    
    # Horizontal diffusion
    horizontal_diffusion!(diagn, progn, model.horizontal_diffusion, model)
    
    # Perform one substep of the N-cycle
    lorenz_ncycle_step!(progn, diagn, diagn.tendencies, dt, model)
    
    # Transform to grid space
    transform!(diagn, progn, 1, model)
    
    # Particle advection (skip first step)
    progn.clock.timestep_counter > 0 && particle_advection!(progn, diagn, model)

    return nothing 
end
# ============================================================================
# MAIN TIME STEPPING LOOP
# ============================================================================

"""$(TYPEDSIGNATURES)
Perform one single time step of `simulation`.
"""
function timestep!(simulation::AbstractSimulation)
    progn, diagn, model = unpack(simulation)
    (; feedback, output) = model
    (; time_stepping) = model
    (; clock) = progn

    # Initialize on first call - now we have diagn!
    if clock.timestep_counter == 0
        initialize!(model.feedback, clock, model)
        
        # Initialize implicit solver with diagn now available
        if !isnothing(model.implicit)
            initialize!(model.implicit, time_stepping.Δt, time_stepping.α, diagn, model)
            set_initialized!(model.implicit)
        end
    end

    # Perform one substep (uses Δt, not 2Δt like Leapfrog)
    (; Δt, Δt_millisec) = time_stepping
    timestep!(progn, diagn, Δt, model)
    timestep!(clock, Δt_millisec)

    # Progress updates and output
    progress!(feedback, progn)
    output!(output, simulation)
    callback!(model.callbacks, progn, diagn, model)
    
    return nothing
end

"""$(TYPEDSIGNATURES)
Main time loop that loops over all time steps."""
function time_stepping!(simulation::AbstractSimulation)          
    (; clock) = simulation.prognostic_variables
    
    for _ in 1:clock.n_timesteps
        timestep!(simulation)
    end
end

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

"""$(TYPEDSIGNATURES)
Computes the time step in [ms]. `Δt_at_T31` is always scaled with the resolution `trunc` 
of the model. In case `adjust_with_output` is true, `Δt_at_T31` is additionally 
adjusted to the closest divisor of `output_dt` so that the output time axis keeps
`output_dt` exactly."""
function get_Δt_millisec(
    Δt_at_T31::Dates.TimePeriod,
    trunc,
    radius,
    adjust_with_output::Bool,
    output_dt::Dates.TimePeriod = DEFAULT_OUTPUT_DT,
)
    # Linearly scale Δt with trunc+1 (which are often powers of two)
    resolution_factor = (DEFAULT_TRUNC+1)/(trunc+1)

    # Radius also affects grid spacing, scale proportionally
    radius_factor = radius/DEFAULT_RADIUS

    # Scale time step
    Δt_at_trunc = Second(Δt_at_T31).value * resolution_factor * radius_factor

    if adjust_with_output && (output_dt > Millisecond(0))
        k = round(Int, Second(output_dt).value / Δt_at_trunc)
        divisors = Primes.divisors(Millisecond(output_dt).value)
        sort!(divisors)
        i = findfirst(x -> x>=k, divisors)
        k_new = isnothing(i) ? k : divisors[i]
        Δt_millisec = Millisecond(round(Int, Millisecond(output_dt).value/k_new))

        # Provide info when time step is significantly changed
        Δt_millisec_unadjusted = round(Int, 1000*Δt_at_trunc)
        Δt_ratio = Δt_millisec.value/Δt_millisec_unadjusted

        if abs(Δt_ratio - 1) > 0.05
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
Change time step of timestepper `L` to `Δt` (unscaled) and disable adjustment to output frequency."""
function set!(
    L::NCycleLorenz,
    Δt::Period,
)
    L.Δt_millisec = Millisecond(Δt)
    L.Δt_sec = L.Δt_millisec.value/1000
    L.Δt = L.Δt_sec/L.radius

    # Recalculate the default time step at resolution T31 to be consistent
    resolution_factor = (L.trunc+1)/(DEFAULT_TRUNC+1)
    L.Δt_at_T31 = Second(round(Int, L.Δt_sec*resolution_factor))

    # Given Δt was manually set, disallow adjustment to output frequency
    L.adjust_with_output = false
    return L
end

# Also allow for keyword arguments
set!(L::NCycleLorenz; Δt::Period) = set!(L, Δt)