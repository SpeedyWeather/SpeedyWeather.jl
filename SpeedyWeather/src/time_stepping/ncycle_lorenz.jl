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


abstract type AbstractNCycleLorenz <: AbstractTimeStepper end

"""
    NCycleLorenz{NF, V} <: AbstractTimeStepper

A semi-implicit Lorenz N-cycle time integration scheme following Hotta et al. (2016).

# Algorithm (per substep of a cycle)
1. G = w*F_E(x) + (1-w)*G     (weighted tendency accumulation)
2. dx = (I - α*Δt*L_I)^(-1) * (G + L_I*x)  (implicit solve)
3. x = x + Δt*dx              (state update)
$(TYPEDFIELDS)
"""
mutable struct NCycleLorenz{NF, V, IntType, S, MS, B} <: AbstractNCycleLorenz
    "[OPTION] Number of steps N in a cycles (3 or 4 recommended, 4 is more stable)"
    steps::IntType

    "[OPTION] Variant: NCycleLorenzA() (default), B, AB, or ABBA"
    variant::V

    "[OPTION] Time step for T31, scale linearly with resolution"
    Δt_at_T31::S

    "[OPTION] Adjust `Δt_at_T31` with the `output_dt` to reach `output_dt` exactly"
    adjust_with_output::B

    "[DERIVED] Time step Δt in milliseconds at specified resolution"
    Δt_millisec::MS

    "[DERIVED] Time step Δt [s] at specified resolution"
    Δt_sec::NF

    "[DERIVED] Time step Δt [s/m] at specified resolution, scaled by 1/radius"
    Δt::NF
end

Adapt.adapt_structure(to, L::NCycleLorenz) = NCycleLorenzCore(L.Δt_millisec, L.Δt_sec, L.Δt)

prognostic_steps(::NCycleLorenz) = 1
tendency_grid_steps(::NCycleLorenz) = 1     # the grid tendencies are only for F though, the G term only needs storing in spectral space
tendency_spectral_steps(::NCycleLorenz) = 2 # to store F, G in Hotta et al. 2016, Eqs 5 & 6

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
        steps = 4,
        variant = NCycleLorenzA(),
        Δt_at_T31 = Minute(40),
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
@inline current_substep(L::NCycleLorenz, clock) = mod(clock.timestep_counter, L.steps)

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
    return weight_coefficient(NF, L.variant, clock.timestep_counter, L.steps)
end

# Type-stable versions with explicit NF
@inline function weight_coefficient(::Type{NF}, ::NCycleLorenzA, i::Integer, N::Integer) where {NF}
    k = mod(i, N)   # current substep
    return k == 0 ? one(NF) : convert(NF, N) / convert(NF, N - k)
end

@inline function weight_coefficient(::Type{NF}, ::NCycleLorenzB, i::Integer, N::Integer) where {NF}
    k = mod(i, N)   # current substep
    return k == 0 ? one(NF) : convert(NF, N) / convert(NF, k)
end

function update_prognostic!(
        var::AbstractArray,
        tendency::AbstractArray,
        vars::Variables,
        time_stepping::NCycleLorenz,
        ::AbstractModel,
    )
    (; Δt) = time_stepping
    w = weight_coefficient(time_stepping, vars.prognostic.clock)
    F = get_step(tendency, 1)   # 1st step is for F at current time step
    G = get_step(tendency, 2)   # 2nd step is for G which accumulates the weighted tendencies

    launch!(
        architecture(var), LinearWorkOrder, size(var), ncycle_lorenz_kernel!,
        var, G, F, w, Δt,
    )
    return nothing
end

@kernel inbounds = true function ncycle_lorenz_kernel!(
        var, G, F, w, Δt,
    )
    lmk = @index(Global, Linear)
    G[lmk] = w * F[lmk] + (1 - w) * G[lmk]
    var[lmk] = var[lmk] + Δt * G[lmk]
end