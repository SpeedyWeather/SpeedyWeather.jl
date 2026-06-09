export NCycleLorenz, NCycleLorenzA, NCycleLorenzB, NCycleLorenzAB, NCycleLorenzABBA

"""Abstract type for Lorenz N-cycle variants following Hotta et al. (2016)."""
abstract type NCycleLorenzVariant end

"""Version A: weights w_k = N/(N-k) for k=1,...,N-1; w_0 = 1"""
struct NCycleLorenzA <: NCycleLorenzVariant end
@inline subcycles(::NCycleLorenzA) = 1

"""Version B: weights w_k = N/k for k=1,...,N-1; w_0 = 1"""
struct NCycleLorenzB <: NCycleLorenzVariant end
@inline subcycles(::NCycleLorenzB) = 1

"""Version AB: alternates A and B every N steps"""
struct NCycleLorenzAB <: NCycleLorenzVariant end
@inline subcycles(::NCycleLorenzAB) = 2         # one cycle of A one of B

"""Version ABBA: uses A-B-B-A sequence (only for N=4, provides 4th order accuracy)"""
struct NCycleLorenzABBA <: NCycleLorenzVariant end
@inline subcycles(::NCycleLorenzABBA) = 4       # 4 cycles, A, B, B, then A again

abstract type AbstractNCycleLorenz <: AbstractTimeStepper end

"""
    NCycleLorenz{NF, V, ...} <: AbstractTimeStepper

A semi-implicit Lorenz N-cycle time integration scheme following Hotta et al. (2016).

# Algorithm (per substep of a cycle)
1. G = w*F_E(x) + (1-w)*G                   (weighted tendency accumulation)
2. dx = (I - α*Δt*L_I)^(-1) * (G + L_I*x)   (implicit solve)
3. x = x + Δt*dx                            (state update)
$(TYPEDFIELDS)
"""
mutable struct NCycleLorenz{NF, V, IntType, S, MS, B} <: AbstractNCycleLorenz
    "[OPTION] Number of steps N in a cycle (3 or 4 recommended, 4 is more stable)"
    steps::IntType

    "[OPTION] Variant: NCycleLorenzA() (default), B, AB, or ABBA"
    variant::V

    "[OPTION] Time step for T31, scale linearly with resolution"
    Δt_at_T31::S

    "[OPTION] Adjust `Δt_at_T31` with the `interval` to reach `interval` exactly"
    adjust_with_output::B

    "[DERIVED] Time step Δt in milliseconds at specified resolution"
    Δt_millisec::MS

    "[DERIVED] Time step Δt [s] at specified resolution"
    Δt_sec::NF

    "[DERIVED] Time step Δt [s/m] at specified resolution, scaled by 1/radius"
    Δt::NF
end

Adapt.adapt_structure(to, L::NCycleLorenz) = Adapt.adapt_structure(to, NCycleLorenzCore(L.Δt_millisec, L.Δt_sec, L.Δt))

prognostic_steps(::NCycleLorenz) = 1
tendency_grid_steps(::NCycleLorenz) = 1     # the grid tendencies are only for F though, the G term only needs storing in spectral space
tendency_spectral_steps(::NCycleLorenz) = 2 # to store F, G in Hotta et al. 2016, Eqs 5 & 6

# Most components just use the 1st tendency step and only the timestepping itself needs the 2nd so the default step 1 is used throughout
# Exceptions here to be explicit: While two tendencies F, G are used, only G retains memory to the next time step, therefore reset F to zero for accumulation
@inline which_tendency_step(var, ::AbstractNCycleLorenz, ::ResetTendencies) = 1

# dispatch over timestepper to decide between implicit or explicit diffusion
@inline implicit_diffusion(::AbstractHorizontalDiffusion, ::Nothing, ::AbstractNCycleLorenz) = true
@inline implicit_diffusion(::AbstractHorizontalDiffusion, ::AbstractImplicit, ::AbstractNCycleLorenz) = true

# dispatch over time stepper here so that other time stepper can change the order
@noinline function diffusion_and_implicit!(vars, ::AbstractNCycleLorenz, ::AbstractImplicit, model)
    horizontal_diffusion!(vars, model)
    implicit_correction!(vars, model)
    return nothing
end

struct NCycleLorenzCore{NF, MS} <: AbstractNCycleLorenz
    Δt_millisec::MS
    Δt_sec::NF
    Δt::NF
end

Adapt.@adapt_structure NCycleLorenzCore

"""$(TYPEDSIGNATURES)
Generator function for NCycleLorenz struct using `spectral_grid` for resolution."""
function NCycleLorenz(
        spectral_grid::SpectralGrid;
        steps = 3,
        variant = NCycleLorenzA(),
        Δt_at_T31 = Minute(30),
        adjust_with_output = true,
        radius = DEFAULT_RADIUS,
    )
    (; NF, trunc) = spectral_grid

    # compute time step
    Δt_millisec::Millisecond = get_Δt_millisec(Second(Δt_at_T31), trunc, DEFAULT_RADIUS, adjust_with_output)
    Δt_sec::NF = Δt_millisec.value / 1000
    Δt::NF = Δt_sec / radius

    return NCycleLorenz(steps, variant, Second(Δt_at_T31), adjust_with_output, Δt_millisec, Δt_sec, Δt)
end

"""$(TYPEDSIGNATURES)
Get current substep within N-cycle (0 to N-1) from the clock."""
@inline current_substep(L::NCycleLorenz, clock) = mod(clock.step_counter, L.steps)

"""$(TYPEDSIGNATURES)
Initialize NCycleLorenz time stepper."""
function initialize!(L::NCycleLorenz, model::AbstractModel)
    if L.variant isa NCycleLorenzABBA && L.steps != 4          # Validate compatibility
        @warn "N-Cycle Lorenz with ABBA variant is for N=4 (4th order accurate), but N=$(L.steps). Consider steps=4 or variant A/B/AB."
    end

    calculate_Δt!(L, model)
    return nothing
end

# Weight coefficient functions (Hotta et al. 2016, Eqs 2-3, 7-8)
"""$(TYPEDSIGNATURES)
Compute weight coefficient w for current substep."""
function weight_coefficient(L::NCycleLorenz{NF}, clock::Clock) where {NF}
    return weight_coefficient(NF, L.variant, clock.step_counter, L.steps)
end

"""$(TYPEDSIGNATURES) Weight coefficient of the A-variant of the N-Cycle Lorenz time stepping scheme."""
@inline function weight_coefficient(::Type{NF}, V::NCycleLorenzA, i, N::Integer) where {NF}
    k = mod(i, subcycles(V) * N)   # current substep
    return ifelse(k == 0, one(NF), convert(NF, N) / convert(NF, N - k))
end

"""$(TYPEDSIGNATURES) Weight coefficient of the B-variant of the N-Cycle Lorenz time stepping scheme."""
@inline function weight_coefficient(::Type{NF}, V::NCycleLorenzB, i, N::Integer) where {NF}
    k = mod(i, subcycles(V) * N)   # current substep
    return ifelse(k == 0, one(NF), convert(NF, N) / convert(NF, k))
end

"""$(TYPEDSIGNATURES) Weight coefficient of the A-variant of the N-Cycle Lorenz time stepping scheme."""
@inline function weight_coefficient(::Type{NF}, V::NCycleLorenzAB, i::Integer, N::Integer) where {NF}
    k = mod(i, subcycles(V) * N)   # current substep
    variant = k < N ? NCycleLorenzA() : NCycleLorenzB()
    return weight_coefficient(NF, variant, i, N)
end

@inline function weight_coefficient(::Type{NF}, V::NCycleLorenzABBA, i::Integer, N::Integer) where {NF}
    k = mod(i, subcycles(V) * N)   # current substep across all subcycles
    variant = k < N || k >= (subcycles(V) - 1) * N ? NCycleLorenzA() : NCycleLorenzB()
    return weight_coefficient(NF, variant, i, N)
end

function update_prognostic!(
        var::AbstractArray,
        tendency::AbstractArray,
        vars::Variables,
        time_stepping::NCycleLorenz,
        implicit::Union{Nothing, AbstractImplicit},
        ::AbstractModel,
    )
    (; Δt) = time_stepping
    w = weight_coefficient(time_stepping, vars.prognostic.clock)

    # with an implicit solver the tendency_average_kernel! has to be computed
    # before the implicit solver, so the responsibility is left therein
    # and execute the prognostic update here only, dispatched over the type of implicit
    # without an implicit solver we compute both tendency average and update
    # here in one kernel, notation following largely Hotta et al. 2016
    F = get_step(tendency, 1)   # tendency of current time step (or weighted + implicitly corrected tendency)
    G = get_step(tendency, 2)   # accumulated weighted tendencies

    launch!(
        architecture(var), LinearWorkOrder, size(var), ncycle_lorenz_kernel!,
        var, F, G, w, Δt
    )
    return nothing
end

@kernel inbounds = true function ncycle_lorenz_kernel!(
        var, F, G, w, Δt,
    )
    lmk = @index(Global, Linear)
    G[lmk] = w * F[lmk] + (1 - w) * G[lmk]  # Hotta et al. 2016 eq (5)
    var[lmk] = var[lmk] + Δt * G[lmk]       # and equation (6)
end
