# ============================================================================
# LORENZ N-CYCLE TIME STEPPING
# Based on Hotta et al. (2016), Monthly Weather Review
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

"""Version ABBA: uses A-B-B-A sequence every 2N steps (N=4 only, 4th-order accurate for
nonlinear systems per Hotta et al. 2016 section 5c)"""
struct NCycleLorenzABBA <: NCycleLorenzVariant end

"""
    NCycleLorenz{NF, V} <: AbstractTimeStepper

A semi-implicit Lorenz N-cycle time integration scheme following Hotta et al. (2016).

# Algorithm (per substep)
1. G   ← w*F_E(x) + (1-w)*G              (weighted F_E accumulation, Eq. 19; G preserved)
2. δx  = (I - α*Δt*L_I)^(-1) * (G + L_I*x)  (implicit correction, Eq. 20; G untouched)
3. x   ← x + Δt*δx                       (state update, Eq. 21)

# Key differences from Leapfrog
- Self-starting: no special first steps needed
- Two time levels only: current state (index 1) + tendency accumulator (index 2)
- No computational mode: no Robert-Asselin filtering needed
- Stable for dissipative systems: physics and diffusion can be treated uniformly
- ξ = α*Δt (not α*2Δt as in leapfrog)

$(TYPEDFIELDS)
"""
@kwdef mutable struct NCycleLorenz{NF<:AbstractFloat, V<:NCycleLorenzVariant} <: AbstractTimeStepper
    "[DERIVED] Spectral resolution (max degree of spherical harmonics)"
    trunc::Int

    "[CONST] Number of time steps stored (2: current state + tendency accumulator)"
    nsteps::Int = 2

    "[OPTION] Number of substeps per cycle N (3 or 4 recommended; 4 is more stable)"
    cycles::Int = 4

    "[OPTION] Variant: NCycleLorenzA() (default), NCycleLorenzB(), NCycleLorenzAB(), NCycleLorenzABBA()"
    variant::V = NCycleLorenzA()

    "[OPTION] Centering parameter: 0.5=Crank-Nicolson (2nd order), 1.0=Backward Euler (1st order)"
    α::NF = 0.5

    "[OPTION] Time step in minutes for T31, scaled linearly to `trunc`"
    Δt_at_T31::Second = Minute(30)

    "[OPTION] Adjust `Δt_at_T31` with `output_dt` to reach `output_dt` exactly"
    adjust_with_output::Bool = true

    "[DERIVED] Radius of sphere [m], set in `initialize!` from `planet.radius`"
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
Initialize NCycleLorenz time stepper: calculate final Δt based on resolution and
output frequency. Unlike Leapfrog, no special first-step handling is needed."""
function initialize!(L::NCycleLorenz, model::AbstractModel)
    (; output_dt) = model.output
    (; radius) = model.planet

    if L.variant isa NCycleLorenzABBA && L.cycles != 4
        @warn "ABBA variant is designed for N=4 (Hotta et al. 2016, section 5c), but N=$(L.cycles). Consider cycles=4."
    end

    L.Δt_millisec = get_Δt_millisec(L.Δt_at_T31, L.trunc, radius, L.adjust_with_output, output_dt)
    L.Δt_sec = L.Δt_millisec.value/1000
    L.Δt = L.Δt_sec/radius

    n = round(Int, Millisecond(output_dt).value/L.Δt_millisec.value)
    nΔt = n*L.Δt_millisec
    if nΔt != output_dt
        @warn "$n steps of Δt=$(L.Δt_millisec.value)ms yield output every $(nΔt.value)ms, but output_dt=$(output_dt.value)ms"
    end

    return nothing
end

# ============================================================================
# WEIGHT COEFFICIENTS (Hotta et al. 2016, Eqs 2-3, 7-8)
# ============================================================================

"""$(TYPEDSIGNATURES)
Get the current substep index k within a single N-cycle (0 to N-1)."""
current_substep(L::NCycleLorenz, clock) = mod(clock.timestep_counter, L.cycles)

"""$(TYPEDSIGNATURES)
For alternating variants (AB, ABBA), determine which base variant (A or B) to use
at the current timestep_counter.

- AB:   repeats [A, B, A, B, ...] where each letter is a full N-cycle
- ABBA: repeats [A, B, B, A, A, B, B, A, ...] where each letter is a full N-cycle
        (only meaningful for N=4, Hotta et al. 2016 section 5c)
"""
function current_base_variant(L::NCycleLorenz{NF, NCycleLorenzAB}, clock) where NF
    # Which N-cycle are we in? (integer division)
    cycle_index = div(clock.timestep_counter, L.cycles)
    return iseven(cycle_index) ? NCycleLorenzA() : NCycleLorenzB()
end

function current_base_variant(L::NCycleLorenz{NF, NCycleLorenzABBA}, clock) where NF
    # ABBA repeats with period 4N: cycle indices 0,3 → A; 1,2 → B
    cycle_index = mod(div(clock.timestep_counter, L.cycles), 4)
    return (cycle_index == 0 || cycle_index == 3) ? NCycleLorenzA() : NCycleLorenzB()
end

"""$(TYPEDSIGNATURES)
Compute weight coefficient w for the current substep."""
function weight_coefficient(L::NCycleLorenz{NF}, clock) where NF
    k = current_substep(L, clock)
    return weight_coefficient(NF, L.variant, k, L.cycles, clock)
end

# Version A: w_0 = 1, w_k = N/(N-k)  (Hotta et al. 2016, Eqs 2-3)
@inline function weight_coefficient(::Type{NF}, ::NCycleLorenzA, k::Int, N::Int, clock) where NF
    return k == 0 ? one(NF) : convert(NF, N) / convert(NF, N - k)
end

# Version B: w_0 = 1, w_k = N/k  (Hotta et al. 2016, Eqs 7-8)
@inline function weight_coefficient(::Type{NF}, ::NCycleLorenzB, k::Int, N::Int, clock) where NF
    return k == 0 ? one(NF) : convert(NF, N) / convert(NF, k)
end

# AB and ABBA: delegate to the appropriate base variant for this cycle
@inline function weight_coefficient(::Type{NF}, ::NCycleLorenzAB, k::Int, N::Int, clock) where NF
    # Determine which base variant we're using for this N-cycle
    cycle_index = div(clock.timestep_counter, N)
    base = iseven(cycle_index) ? NCycleLorenzA() : NCycleLorenzB()
    return weight_coefficient(NF, base, k, N, clock)
end

@inline function weight_coefficient(::Type{NF}, ::NCycleLorenzABBA, k::Int, N::Int, clock) where NF
    cycle_index = mod(div(clock.timestep_counter, N), 4)
    base = (cycle_index == 0 || cycle_index == 3) ? NCycleLorenzA() : NCycleLorenzB()
    return weight_coefficient(NF, base, k, N, clock)
end

# ============================================================================
# HELPER KERNELS FOR STATE UPDATES
# ============================================================================

"""$(TYPEDSIGNATURES)
Accumulate weighted explicit tendency in-place: G ← w*F_E + (1-w)*G  (Hotta et al. 2016, Eq. 19).
`_strip_implicit!` must have been called before this to remove the implicit operator L_I(x) from F,
so that only the explicit part F_E is accumulated here."""
function accumulate_tendency!(
    G::LowerTriangularArray,
    F_explicit::LowerTriangularArray,
    w::Real,
    L::NCycleLorenz{NF},
) where NF
    @boundscheck size(G) == size(F_explicit) || throw(BoundsError())
    w_NF = convert(NF, w)
    launch!(architecture(G), SpectralWorkOrder, size(G),
            accumulate_tendency_kernel!, G, F_explicit, w_NF)
    return nothing
end

@kernel inbounds=true function accumulate_tendency_kernel!(G, F_explicit, @Const(w))
    lmk = @index(Global, Linear)
    G[lmk] = w * F_explicit[lmk] + (1 - w) * G[lmk]
end

"""$(TYPEDSIGNATURES)
Update state in-place: x ← x + Δt*G  (Hotta et al. 2016, Eq. 21)."""
function update_state!(
    x::LowerTriangularArray,
    G::LowerTriangularArray,
    dt::Real,
    L::NCycleLorenz{NF},
) where NF
    @boundscheck size(x) == size(G) || throw(BoundsError())
    dt_NF = convert(NF, dt)
    launch!(architecture(x), SpectralWorkOrder, size(x),
            update_state_kernel!, x, G, dt_NF)
    return nothing
end

@kernel inbounds=true function update_state_kernel!(x, G, @Const(dt))
    lmk = @index(Global, Linear)
    x[lmk] = x[lmk] + dt * G[lmk]
end

# ============================================================================
# CORE LORENZ N-CYCLE SUBSTEP
# ============================================================================

"""$(TYPEDSIGNATURES)
Perform ONE substep of the Lorenz N-cycle for all prognostic variables,
following Hotta et al. (2016). The clock tracks which substep we are on within
the current N-cycle.

Steps per substep:
1. Reset G to zero at the start of each new N-cycle (k=0)
2. Accumulate: G ← w*F_E + (1-w)*G                     (Eq. 19)
3. Implicit correction: G ← (I - α*Δt*L_I)⁻¹*(G + L_I*x) (Eq. 20)
4. State update: x ← x + Δt*G                           (Eq. 21)

Storage convention:
  progn variable index 1 → current state x
  progn variable index 2 → tendency accumulator G
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

    # STEP 0: Reset tendency accumulator G at the start of each new N-cycle
    if current_substep(L, clock) == 0
        for varname in prognostic_variables(model)
            fill!(get_step(getfield(progn, varname), 2), 0)
        end
        for (name, tracer) in model.tracers
            tracer.active || continue
            fill!(get_step(progn.tracers[name], 2), 0)
        end
    end

    # STEP 1: Weighted tendency accumulation  G ← w*F_E + (1-w)*G  (Eq. 19)
    # First strip the implicit operator L_I(x) from the full tendency F = F_E + L_I(x)
    # so that only the explicit part F_E is accumulated. The implicit correction (STEP 2)
    # adds L_I(x_current) back once to form the exact Eq. 20 RHS.
    _strip_implicit!(diagn, progn, model.implicit, model)
    w = weight_coefficient(L, clock)

    for (varname, tendname) in zip(prognostic_variables(model), tendency_names(model))
        var = getfield(progn, varname)
        G   = get_step(var, 2)
        F_explicit = getfield(tend, tendname)
        spectral_truncation!(F_explicit)
        accumulate_tendency!(G, F_explicit, w, L)
    end

    for (name, tracer) in model.tracers
        tracer.active || continue
        G          = get_step(progn.tracers[name], 2)
        F_explicit = tend.tracers_tend[name]
        spectral_truncation!(F_explicit)
        accumulate_tendency!(G, F_explicit, w, L)
    end

    # STEP 2: Semi-implicit correction  δx = (I - α*Δt*L_I)⁻¹*(G_E + L_I*x)  (Eq. 20)
    # First copy G_E → tend so that:
    #   - explicit variables (e.g. vor) have δx = G_E already in tend
    #   - G_E in progn index 2 is left unchanged (not overwritten with δx)
    # Then lorenz_implicit_correction! overwrites only the implicit variables
    # (div_tend, pres_tend for SWE) with the corrected δx.
    for (varname, tendname) in zip(prognostic_variables(model), tendency_names(model))
        G  = get_step(getfield(progn, varname), 2)
        δx = getfield(tend, tendname)
        copyto!(δx, G)
    end
    for (name, tracer) in model.tracers
        tracer.active || continue
        copyto!(tend.tracers_tend[name], get_step(progn.tracers[name], 2))
    end
    lorenz_implicit_correction!(diagn, progn, model.implicit, model)

    # STEP 3: State update  x ← x + Δt*δx  (Eq. 21), reading δx from tend
    for (varname, tendname) in zip(prognostic_variables(model), tendency_names(model))
        var = getfield(progn, varname)
        δx  = getfield(tend, tendname)
        update_state!(get_step(var, 1), δx, dt, L)
    end

    for (name, tracer) in model.tracers
        tracer.active || continue
        var = progn.tracers[name]
        update_state!(get_step(var, 1), tend.tracers_tend[name], dt, L)
    end

    random_process!(progn, model.random_process)
    return nothing
end

# ============================================================================
# MODEL-SPECIFIC TIMESTEP! METHODS
# ============================================================================

"""$(TYPEDSIGNATURES)
Returns the tendency field names corresponding to each prognostic variable."""
tendency_names(model::AbstractModel) =
    tuple((Symbol(var, :_tend) for var in prognostic_variables(model))...)

"""$(TYPEDSIGNATURES)
One substep for the barotropic model."""
function timestep!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    dt::Real,
    model::Barotropic,
)
    model.feedback.nans_detected && return nothing

    fill!(diagn.tendencies, 0, Barotropic)
    dynamics_tendencies!(diagn, progn, 1, model)
    horizontal_diffusion!(diagn, progn, model.horizontal_diffusion, model)
    lorenz_ncycle_step!(progn, diagn, diagn.tendencies, dt, model)
    transform!(diagn, progn, 1, model)
    progn.clock.timestep_counter > 0 && particle_advection!(progn, diagn, model)

    return nothing
end

"""$(TYPEDSIGNATURES)
One substep for the shallow water model."""
function timestep!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    dt::Real,
    model::ShallowWater,
)
    model.feedback.nans_detected && return nothing

    fill!(diagn.tendencies, 0, ShallowWater)
    dynamics_tendencies!(diagn, progn, 1, model)
    horizontal_diffusion!(diagn, progn, model.horizontal_diffusion, model)
    lorenz_ncycle_step!(progn, diagn, diagn.tendencies, dt, model)
    transform!(diagn, progn, 1, model)
    progn.clock.timestep_counter > 0 && particle_advection!(progn, diagn, model)

    return nothing
end

"""$(TYPEDSIGNATURES)
One substep for the primitive equation model."""
function timestep!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    dt::Real,
    model::PrimitiveEquation,
)
    model.feedback.nans_detected && return nothing

    fill!(diagn.tendencies, 0, typeof(model))

    if model.physics
        parameterization_tendencies!(diagn, progn, model)
        ocean_timestep!(progn, diagn, model)
        sea_ice_timestep!(progn, diagn, model)
        land_timestep!(progn, diagn, model)
    end

    if model.dynamics
        forcing!(diagn, progn, 1, model)
        drag!(diagn, progn, 1, model)
        dynamics_tendencies!(diagn, progn, 1, model)
    else
        physics_tendencies_only!(diagn, model)
    end

    horizontal_diffusion!(diagn, progn, model.horizontal_diffusion, model)
    lorenz_ncycle_step!(progn, diagn, diagn.tendencies, dt, model)
    transform!(diagn, progn, 1, model)
    progn.clock.timestep_counter > 0 && particle_advection!(progn, diagn, model)

    return nothing
end

# ============================================================================
# SIMULATION-LEVEL TIMESTEP AND MAIN LOOP
# ============================================================================

"""$(TYPEDSIGNATURES)
Perform one substep of `simulation`, handling first-step initialization,
model integration, output, and callbacks."""
function timestep!(simulation::AbstractSimulation)
    progn, diagn, model = unpack(simulation)
    (; feedback, output) = model
    (; time_stepping) = model
    (; clock) = progn

    # On the very first call: initialize feedback and the implicit solver.
    # The implicit solver is initialized here (not in model initialize!) because
    # it needs diagn.temp_average which is only available after the first transform.
    if clock.timestep_counter == 0
        initialize!(feedback, clock, model)

        if !isnothing(model.implicit)
            # Pass Δt (not 2Δt): the Lorenz N-cycle uses a single time level.
            # α is read from model.implicit itself inside initialize!.
            initialize!(model.implicit, time_stepping.Δt, diagn, model)
            set_initialized!(model.implicit)
        end
    end

    (; Δt, Δt_millisec) = time_stepping
    timestep!(progn, diagn, Δt, model)
    timestep!(clock, Δt_millisec)

    progress!(feedback, progn)
    output!(output, simulation)
    callback!(model.callbacks, progn, diagn, model)

    return nothing
end

"""$(TYPEDSIGNATURES)
Main time loop: advance `simulation` for all scheduled timesteps."""
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
Compute the time step in milliseconds. `Δt_at_T31` is scaled linearly with
spectral resolution `trunc`. If `adjust_with_output` is true, it is additionally
snapped to the nearest divisor of `output_dt` so the output axis is exact."""
function get_Δt_millisec(
    Δt_at_T31::Dates.TimePeriod,
    trunc::Int,
    radius::Real,
    adjust_with_output::Bool,
    output_dt::Dates.TimePeriod = DEFAULT_OUTPUT_DT,
)
    resolution_factor = (DEFAULT_TRUNC + 1) / (trunc + 1)
    radius_factor     = radius / DEFAULT_RADIUS
    Δt_at_trunc       = Second(Δt_at_T31).value * resolution_factor * radius_factor

    if adjust_with_output && (output_dt > Millisecond(0))
        k        = round(Int, Second(output_dt).value / Δt_at_trunc)
        divisors = Primes.divisors(Millisecond(output_dt).value)
        sort!(divisors)
        i        = findfirst(x -> x >= k, divisors)
        k_new    = isnothing(i) ? k : divisors[i]
        Δt_millisec = Millisecond(round(Int, Millisecond(output_dt).value / k_new))

        Δt_millisec_unadjusted = round(Int, 1000*Δt_at_trunc)
        Δt_ratio = Δt_millisec.value / Δt_millisec_unadjusted
        if abs(Δt_ratio - 1) > 0.05
            p  = round(Int, (Δt_ratio - 1)*100)
            ps = p > 0 ? "+" : ""
            @info "Time step changed from $(Δt_millisec_unadjusted)ms to $(Δt_millisec.value)ms ($ps$p%) to match output frequency."
        end
    else
        Δt_millisec = Millisecond(round(Int, 1000*Δt_at_trunc))
    end

    return Δt_millisec
end

"""$(TYPEDSIGNATURES)
Manually set the time step of `L` to `Δt` and disable output-frequency adjustment."""
function set!(L::NCycleLorenz, Δt::Period)
    L.Δt_millisec = Millisecond(Δt)
    L.Δt_sec      = L.Δt_millisec.value / 1000
    L.Δt          = L.Δt_sec / L.radius

    resolution_factor = (L.trunc + 1) / (DEFAULT_TRUNC + 1)
    L.Δt_at_T31       = Second(round(Int, L.Δt_sec * resolution_factor))
    L.adjust_with_output = false
    return L
end

set!(L::NCycleLorenz; Δt::Period) = set!(L, Δt)