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

# Algorithm (per substep k of a cycle)
1. G = (I - w*θₖ*Δt*L)^(-1) * (w*F(x) + (1-w)*G)  (weighted tendency accumulation with implicit solve)
2. x = x + Δt*G                                  (state update)

with F the full tendency, L the linear gravity-wave terms, w the weight and θₖ the implicit
centering of substep k (see `implicit_centering`).
$(TYPEDFIELDS)
"""
mutable struct NCycleLorenz{NF, V, IntType, C, S, MS, B, I, O, L} <: AbstractNCycleLorenz
    "[OPTION] Number of steps N in a cycle (3 or 4 recommended, 4 (default) is more stable and less diffusive)"
    steps::IntType

    "[OPTION] Variant: NCycleLorenzA() (default), B, AB, or ABBA"
    variant::V

    "[OPTION] Implicit centering per substep of the cycle (length `steps`), `nothing` uses `implicit.centering` on every substep"
    implicit_centering::C

    "[OPTION] Time step for T32 resolution, scale linearly with resolution"
    Δt_at_T32::S

    "[OPTION] Adjust `Δt_at_T32` with the `interval` to reach `interval` exactly"
    adjust_with_output::B

    "[DERIVED] Time step Δt in milliseconds at specified resolution"
    Δt_millisec::MS

    "[DERIVED] Time step Δt [s] at specified resolution"
    Δt::NF

    "[OPTION] Implicit solver for the gravity waves, `DefaultImplicit()` is replaced by the model's default"
    implicit::I

    "[OPTION] Time stepper for the ocean (incl sea ice) variables, `nothing` for the same N-cycle as the atmosphere"
    ocean::O

    "[OPTION] Time stepper for the land variables, `nothing` for the same N-cycle as the atmosphere"
    land::L
end

Adapt.adapt_structure(to, L::NCycleLorenz) = Adapt.adapt_structure(to, NCycleLorenzCore(L.Δt_millisec, L.Δt))

prognostic_steps(::NCycleLorenz) = 1
tendency_steps(::NCycleLorenz) = 2          # F and G, e.g. for ocean/land variables stepped directly on the grid
tendency_grid_steps(::NCycleLorenz) = 1     # the grid tendencies are only for F though, the G term only needs storing in spectral space
tendency_spectral_steps(::NCycleLorenz) = 2 # to store F, G in Hotta et al. 2016, Eqs 5 & 6

# ocean and land use the time steppers in the respective fields, `nothing` means the same N-cycle as the atmosphere
time_stepper(L::NCycleLorenz, ::Val{:ocean}) = something(L.ocean, L)
time_stepper(L::NCycleLorenz, ::Val{:land}) = something(L.land, L)

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
    Δt::NF
end

Adapt.@adapt_structure NCycleLorenzCore

"""$(TYPEDSIGNATURES)
Generator function for NCycleLorenz struct using `spectral_grid` for resolution."""
function NCycleLorenz(
        spectral_grid::SpectralGrid;
        steps = 4,
        variant = NCycleLorenzA(),
        implicit_centering = default_implicit_centering(steps, variant),
        Δt_at_T32 = Minute(30),
        adjust_with_output = true,
        implicit = DefaultImplicit(),
        ocean = nothing,    # default can't be NCycleLorenz, creates an infinite loop
        land = nothing,
    )
    (; NF, truncation) = spectral_grid

    # compute time step
    Δt_millisec::Millisecond = get_Δt_millisec(Second(Δt_at_T32), truncation, DEFAULT_RADIUS, adjust_with_output)
    Δt::NF = Δt_millisec.value / 1000

    implicit_centering = isnothing(implicit_centering) ? nothing : NTuple{steps, NF}(implicit_centering)
    return NCycleLorenz(steps, variant, implicit_centering, Second(Δt_at_T32), adjust_with_output, Δt_millisec, Δt, implicit, ocean, land)
end

"""$(TYPEDSIGNATURES)
Default implicit centering θₖ per substep k of the cycle, `nothing` (= `implicit.centering` on
every substep) if none is defined for the cycle length `steps` and `variant`. The centerings
satisfy the second-order condition (effective centering over the cycle α_eff = 1/2, see
`effective_centering`) and were chosen in a linear stability analysis to be stable with
Doppler-shifted gravity waves, Coriolis coupling and a reference layer thickness larger than the
actual one (e.g. over orography) up to the explicit stability limit of the cycle, except for the
3-cycle variant A which is slightly off-centred (α_eff ≈ 0.56) for robustness."""
default_implicit_centering(steps, variant) = nothing
default_implicit_centering(steps, ::NCycleLorenzA) =
    steps == 4 ? (2.06, 0.93, 1.02, 0.09) : steps == 3 ? (-1.03, 0.7, 1.05) : nothing
default_implicit_centering(steps, ::NCycleLorenzB) = steps == 3 ? (1.63, 0.55, -0.55) : nothing

"""$(TYPEDSIGNATURES)
Effective implicit centering α_eff of the Lorenz N-cycle with weights of `variant` and per-substep
centerings `θ` (a tuple of length N, N = 3 or 4) for the linear terms: the cycle's amplification
factor of the linear terms agrees to second order with a θ-method over the whole cycle (time step
NΔt) with centering α_eff. α_eff = 1/2 is second-order accurate (as Crank-Nicolson), α_eff > 1/2
damps gravity waves at second order like an off-centred Crank-Nicolson. From the expansion of the
cycle in Δt: α_eff = 1/2 + κ c⋅θ / N² with c, κ given per N and variant."""
function effective_centering(variant::NCycleLorenzVariant, θ::NTuple{N}) where {N}
    c, κ = effective_centering_coefficients(variant, Val(N))
    return 1 / 2 + κ * sum(c .* θ) / N^2
end

effective_centering_coefficients(::NCycleLorenzA, ::Val{3}) = (1, -1, 2), 3 / 2
effective_centering_coefficients(::NCycleLorenzB, ::Val{3}) = (0, 1, 1), 3 / 2
effective_centering_coefficients(::NCycleLorenzA, ::Val{4}) = (0, 1, -1, 1), 4
effective_centering_coefficients(::NCycleLorenzB, ::Val{4}) = (0, 1, 1, 1), 4 / 3

"""$(TYPEDSIGNATURES)
Implicit centering for the current substep of the cycle, `implicit.centering` if
`time_stepping.implicit_centering` is `nothing`."""
implicit_centering(::AbstractNCycleLorenz, implicit, clock) = implicit.centering
function implicit_centering(L::NCycleLorenz, implicit, clock)
    isnothing(L.implicit_centering) && return implicit.centering
    return L.implicit_centering[current_substep(L, clock) + 1]
end

"""$(TYPEDSIGNATURES)
Get current substep within N-cycle (0 to N-1) from the clock."""
@inline current_substep(L::NCycleLorenz, clock) = mod(clock.step_counter, L.steps)

"""$(TYPEDSIGNATURES)
Initialize NCycleLorenz time stepper."""
function initialize!(L::NCycleLorenz, model::AbstractModel)
    calculate_Δt!(L, model)
    return nothing
end

# Weight coefficient functions (Hotta et al. 2016, Eqs 2-3, 7-8)
"""$(TYPEDSIGNATURES)
Compute weight coefficient w for current substep in the NCycleLorenz time stepping scheme."""
function weight_coefficient(L::NCycleLorenz{NF}, clock::Clock) where {NF}
    return weight_coefficient(NF, L.variant, clock.step_counter, L.steps)
end

"""$(TYPEDSIGNATURES) Weight coefficient of the A-variant of the N-Cycle Lorenz time stepping scheme."""
@inline function weight_coefficient(::Type{NF}, V::NCycleLorenzA, i, N::Integer) where {NF}
    k = one(NF) * mod(i, subcycles(V) * N)   # current substep
    return ifelse(k == 0, one(NF), NF(N) / (NF(N) - k))
end

"""$(TYPEDSIGNATURES) Weight coefficient of the B-variant of the N-Cycle Lorenz time stepping scheme."""
@inline function weight_coefficient(::Type{NF}, V::NCycleLorenzB, i, N::Integer) where {NF}
    k = one(NF) * mod(i, subcycles(V) * N)   # current substep
    return ifelse(k == 0, one(NF), NF(N) / k)
end

"""$(TYPEDSIGNATURES) Weight coefficient of the AB-variant of the N-Cycle Lorenz time stepping scheme."""
@inline function weight_coefficient(::Type{NF}, V::NCycleLorenzAB, i::Integer, N::Integer) where {NF}
    k = mod(i, subcycles(V) * N)   # current substep
    variant = ifelse(k < N, NCycleLorenzA(), NCycleLorenzB())
    return weight_coefficient(NF, variant, i, N)
end

"""$(TYPEDSIGNATURES) Weight coefficient of the ABBA-variant of the N-Cycle Lorenz time stepping scheme."""
@inline function weight_coefficient(::Type{NF}, V::NCycleLorenzABBA, i::Integer, N::Integer) where {NF}
    k = mod(i, subcycles(V) * N)   # current substep across all subcycles
    variant = ifelse((k < N) | (k >= (subcycles(V) - 1) * N), NCycleLorenzA(), NCycleLorenzB())
    return weight_coefficient(NF, variant, i, N)
end

function update_prognostic!(
        var::AbstractArray,
        tendency::AbstractArray,
        clock::Clock,
        time_stepping::NCycleLorenz,
        implicit::Union{Nothing, AbstractImplicit},
        ::AbstractModel,
        scale::Real = 1,
    )
    (; Δt) = time_stepping
    Δt /= oftype(Δt, scale)     # scale on the fly
    w = weight_coefficient(time_stepping, clock)

    # notation following largely Hotta et al. 2016
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
