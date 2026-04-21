const DEFAULT_NSTEPS = 2
export LorenzNCycle

"""Lorenz N-cycle semi-implicit time stepping (Lorenz 1971, Hotta et al. 2016 MWR).
Two arrays are stored per variable: step 1 = current state x, step 2 = accumulated
explicit tendency G. One tendency evaluation per step; self-starting.

Versions A and B differ in their blending weights w^k:
  Version A: w^0 = 1, w^k = N/(N-k) for k = 1..N-1
  Version B: w^0 = 1, w^k = N/k     for k = 1..N-1

At each step: G ← w·F^E(x) + (1-w)·G, then solve implicitly for δx, then x ← x + Δt·δx.
N = 4 is recommended (Hotta et al. 2016).

$(TYPEDFIELDS)"""
@kwdef mutable struct LorenzNCycle{NF<:AbstractFloat} <: AbstractTimeStepper
    "[DERIVED] Spectral resolution (max degree of spherical harmonics)"
    trunc::Int

    "[CONST] Number of time levels stored: 1 = current state x, 2 = accumulated G"
    nsteps::Int = 2

    "[OPTION] Lorenz N-cycle order (N=4 recommended)"
    N::Int = 4

    "[OPTION] Version A (1) or Version B (2) weights"
    version::Int = 2

    "[DERIVED] Current position in the N-cycle, 0-indexed, increments mod N each step"
    k::Int = 0

    "[OPTION] Time step in minutes for T31, scaled linearly to `trunc`"
    Δt_at_T31::Second = Minute(40)

    "[OPTION] Adjust `Δt_at_T31` with `output_dt` to reach it exactly in integer steps"
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

"""$(TYPEDSIGNATURES)
The Lorenz N-cycle blending weight w for the current cycle position k."""
function lorenz_weight(L::LorenzNCycle{NF}) where NF
    (; k, N, version) = L
    k == 0 && return one(NF)
    return version == 1 ? NF(N) / NF(N - k) : NF(N) / NF(k)
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
Generator using `spectral_grid` for resolution information."""
function LorenzNCycle(spectral_grid::SpectralGrid; kwargs...)
    (; NF, trunc) = spectral_grid
    return LorenzNCycle{NF}(; trunc, kwargs...)
end

"""$(TYPEDSIGNATURES)
Initialize `LorenzNCycle`: recalculate the time step from the model's output frequency
and planet radius, and reset the N-cycle counter k to 0."""
function initialize!(L::LorenzNCycle, model::AbstractModel)
    (; output_dt) = model.output
    (; radius) = model.planet

    L.radius = radius
    L.Δt_millisec = get_Δt_millisec(L.Δt_at_T31, L.trunc, radius, L.adjust_with_output, output_dt)
    L.Δt_sec = L.Δt_millisec.value/1000
    L.Δt = L.Δt_sec/radius
    L.k = 0

    n = round(Int, Millisecond(output_dt).value/L.Δt_millisec.value)
    nΔt = n*L.Δt_millisec
    if nΔt != output_dt
        @warn "$n steps of Δt = $(L.Δt_millisec.value)ms yield output every $(nΔt.value)ms, but output_dt = $(Millisecond(output_dt).value/1000)s"
    end
    return nothing
end

"""$(TYPEDSIGNATURES)
Manually set the time step of `L` to `Δt` and disable adjustment to output frequency."""
function set!(L::AbstractTimeStepper, Δt::Period)
    L.Δt_millisec = Millisecond(Δt)
    L.Δt_sec = L.Δt_millisec.value/1000
    L.Δt = L.Δt_sec/L.radius
    resolution_factor = (L.trunc+1)/(DEFAULT_TRUNC+1)
    L.Δt_at_T31 = Second(round(Int, L.Δt_sec*resolution_factor))
    L.adjust_with_output = false
    return L
end

set!(L::AbstractTimeStepper; Δt::Period) = set!(L, Δt)

function Adapt.adapt_structure(to, L::LorenzNCycle)
    return (; Δt=L.Δt, Δt_sec=L.Δt_sec, Δt_millisec=L.Δt_millisec)
end

# ---------------------------------------------------------------------------
# Kernels
# ---------------------------------------------------------------------------

"""Kernel for purely explicit variables (vor, humid, tracers).
Blends the current explicit tendency with stored G, updates state x, stores new G:
  G_new = w·F^E + (1-w)·G_old
  x    += Δt·G_new"""
@kernel inbounds=true function lorenz_step_kernel!(x, G, tendency, @Const(dt), @Const(w))
    lmk = @index(Global, Linear)
    g_new = w * tendency[lmk] + (1 - w) * G[lmk]
    x[lmk] = x[lmk] + dt * g_new
    G[lmk] = g_new
end

"""Kernel for implicitly corrected variables (div, pres, temp for PrimitiveEquation).
`implicit_correction!` performs the G accumulation and implicit solve, writing δx into
`tendency`. This kernel only does: x ← x + Δt·δx."""
@kernel inbounds=true function lorenz_update_kernel!(x, tendency, @Const(dt))
    lmk = @index(Global, Linear)
    x[lmk] = x[lmk] + dt * tendency[lmk]
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

tendency_names(model::AbstractModel) =
    tuple((Symbol(var, :_tend) for var in prognostic_variables(model))...)

"""Variables handled by `implicit_correction!`.
For these, the G blend is done inside that function, not in `lorenz_step!`."""
implicit_variables(::AbstractModel)     = ()
implicit_variables(::ShallowWater)      = (:div, :pres)
implicit_variables(::PrimitiveEquation) = (:div, :temp, :pres)

# ---------------------------------------------------------------------------
# lorenz_step!  —  advance all prognostic variables one Lorenz step
# ---------------------------------------------------------------------------

"""$(TYPEDSIGNATURES)
Advance all prognostic variables of `progn` by one Lorenz N-cycle step.

Explicit variables (vor, humid, tracers): G blend + state update in `lorenz_step_kernel!`.
Implicit variables (div, pres, temp): `implicit_correction!` already did the G blend
and wrote δx into `tend`; here just x ← x + Δt·δx via `lorenz_update_kernel!`."""
function lorenz_step!(
    progn::PrognosticVariables,
    tend::Tendencies,
    dt::Real,
    w::Real,
    model::AbstractModel,
)
    impl_vars = implicit_variables(model)

    for (varname, tendname) in zip(prognostic_variables(model), tendency_names(model))
        var      = getfield(progn, varname)
        x        = get_step(var, 1)
        var_tend = getfield(tend, tendname)
        SpeedyTransforms.spectral_truncation!(var_tend)

        if varname ∈ impl_vars
            launch!(architecture(x), SpectralWorkOrder, size(x),
                    lorenz_update_kernel!, x, var_tend, dt)
        else
            G = get_step(var, 2)
            launch!(architecture(x), SpectralWorkOrder, size(x),
                    lorenz_step_kernel!, x, G, var_tend, dt, w)
        end
    end

    # Tracers are always explicit
    for (name, tracer) in model.tracers
        tracer.active || continue
        var      = progn.tracers[name]
        x        = get_step(var, 1)
        G        = get_step(var, 2)
        var_tend = tend.tracers_tend[name]
        SpeedyTransforms.spectral_truncation!(var_tend)
        launch!(architecture(x), SpectralWorkOrder, size(x),
                lorenz_step_kernel!, x, G, var_tend, dt, w)
    end

    random_process!(progn, model.random_process)
    return nothing
end

# ---------------------------------------------------------------------------
# timestep!  —  one complete model time step
# ---------------------------------------------------------------------------

"""$(TYPEDSIGNATURES)
Single Lorenz N-cycle step for the `Barotropic` model."""
function timestep!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    dt::Real,
    w::Real,
    model::Barotropic,
)
    model.feedback.nans_detected && return nothing
    fill!(diagn.tendencies, 0, Barotropic)

    dynamics_tendencies!(diagn, progn, 1, model)
    horizontal_diffusion!(diagn, progn, model.horizontal_diffusion, model)
    lorenz_step!(progn, diagn.tendencies, dt, w, model)
    transform!(diagn, progn, 1, model)
    particle_advection!(progn, diagn, model)
    return nothing
end

"""$(TYPEDSIGNATURES)
Single Lorenz N-cycle step for the `ShallowWater` model."""
function timestep!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    dt::Real,
    w::Real,
    model::ShallowWater,
)
    model.feedback.nans_detected && return nothing
    fill!(diagn.tendencies, 0, ShallowWater)

    dynamics_tendencies!(diagn, progn, 1, model)
    horizontal_diffusion!(diagn, progn, model.horizontal_diffusion, model)
    implicit_correction!(diagn, progn, model.implicit, model, w)
    lorenz_step!(progn, diagn.tendencies, dt, w, model)
    transform!(diagn, progn, 1, model)
    particle_advection!(progn, diagn, model)
    return nothing
end

"""$(TYPEDSIGNATURES)
Single Lorenz N-cycle step for `PrimitiveEquation` models."""
function timestep!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    dt::Real,
    w::Real,
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
        dynamics_tendencies!(diagn, progn, 1, model)
        horizontal_diffusion!(diagn, progn, model.horizontal_diffusion, model)
        implicit_correction!(diagn, progn, model.implicit, model, w)
    else
        physics_tendencies_only!(diagn, model)
        horizontal_diffusion!(diagn, progn, model.horizontal_diffusion, model)
    end
    lorenz_step!(progn, diagn.tendencies, dt, w, model)
    transform!(diagn, progn, 1, model)
    particle_advection!(progn, diagn, model)
    return nothing
end

# ---------------------------------------------------------------------------
# Simulation-level stepping
# ---------------------------------------------------------------------------

"""$(TYPEDSIGNATURES)
Dispatch to `first_timesteps!` on the first call, `later_timestep!` thereafter."""
function timestep!(simulation::AbstractSimulation)
    (; clock) = simulation.prognostic_variables
    if clock.timestep_counter == 0
        first_timesteps!(simulation)
    else
        later_timestep!(simulation)
    end
end

"""$(TYPEDSIGNATURES)
Initialise the implicit solver and run the first Lorenz N-cycle step.

The Lorenz N-cycle is self-starting: at k=0 the weight w=1, so the first step
is a forward-Euler step with semi-implicit correction — no special initialisation
(Euler half-step, etc.) is required."""
function first_timesteps!(simulation::AbstractSimulation)
    progn, diagn, model = unpack(simulation)
    (; clock) = progn
    clock.n_timesteps == 0 && return nothing

    (; implicit, time_stepping) = model
    (; Δt) = time_stepping

    # Lorenz uses ξ = α·Δt (not 2Δt as in leapfrog)
    initialize!(implicit, Δt, diagn, model)
    set_initialized!(implicit)

    later_timestep!(simulation)
    initialize!(model.feedback, clock, model)
    return nothing
end

"""$(TYPEDSIGNATURES)
One normal Lorenz N-cycle step: compute w from the current cycle position k,
call `timestep!`, increment k mod N, advance the clock, and handle output/callbacks."""
function later_timestep!(simulation::AbstractSimulation)
    progn, diagn, model = unpack(simulation)
    (; feedback, output) = model
    (; Δt, Δt_millisec) = model.time_stepping
    (; time_stepping) = model
    (; clock) = progn

    w = lorenz_weight(time_stepping)
    timestep!(progn, diagn, Δt, w, model)
    time_stepping.k = mod(time_stepping.k + 1, time_stepping.N)
    timestep!(clock, Δt_millisec)

    progress!(feedback, progn)
    output!(output, simulation)
    callback!(model.callbacks, progn, diagn, model)
    return nothing
end

"""$(TYPEDSIGNATURES)
Main time loop over all time steps."""
function time_stepping!(simulation::AbstractSimulation)
    (; clock) = simulation.prognostic_variables
    for _ in 1:clock.n_timesteps
        timestep!(simulation)
    end
end
