"""Abstract supertype for the time integration of relaxation terms `dx/dt = -(x - x₀)/τ` that relax
a variable `x` towards a target value `x₀` on a time scale `τ`, e.g. the surface fluxes
`flux = -conductance * (x - x_surface)` or the `SpeedLimitDrag` relaxing the wind speed towards
the speed limit. Subtypes define `relaxation(scheme, z)` with `z = Δt/τ` the ratio of
time step to relaxation time scale, a factor that multiplies the explicit flux (or tendency)."""
abstract type AbstractRelaxationIntegration end
Base.broadcastable(scheme::AbstractRelaxationIntegration) = Ref(scheme)

export ExplicitRelaxation, ImplicitRelaxation, ExponentialRelaxation, RationalRelaxation

"""Explicit (forward) time integration of a relaxation, unstable for `Δt/τ > 2`."""
struct ExplicitRelaxation <: AbstractRelaxationIntegration end

"""Implicit (backward Euler) time integration of a relaxation, damping the explicit
flux or tendency by `1/(1 + Δt/τ)`. Unconditionally stable."""
struct ImplicitRelaxation <: AbstractRelaxationIntegration end

"""Exponential time differencing of a relaxation, damping the explicit flux or tendency by
`(1 - exp(-Δt/τ))/(Δt/τ)`. Exact for a fixed target value, unconditionally stable."""
struct ExponentialRelaxation <: AbstractRelaxationIntegration end

"""Rational approximation of `ExponentialRelaxation`, damping the explicit flux or tendency by
`(6 + z)/(6 + 4z + z²)` with `z = Δt/τ`. Matches `(1 - exp(-z))/z` to second order
in `z` and as `1/z` for large `z`, at the cost of `ImplicitRelaxation`. Unconditionally
stable but overshoots the target value slightly (≤ 10%) for `z > 3`."""
struct RationalRelaxation <: AbstractRelaxationIntegration end

# default time integration of the surface fluxes
const DEFAULT_FLUX_INTEGRATION = RationalRelaxation

"""$(TYPEDSIGNATURES)
Factor φ(z) multiplying an explicit relaxation flux or tendency for `z = Δt/τ` given time integration `scheme`."""
@inline relaxation(::ExplicitRelaxation, z) = one(z)
@inline relaxation(::ImplicitRelaxation, z) = 1 / (1 + z)
@inline relaxation(::RationalRelaxation, z) = (6 + z) / (6 + 4z + z^2)

# numerator and denominator guarded separately (1/1 for small z), not with ifelse(small, 1, -expm1(-z)/z),
# whose untaken branch is 0/0 for z = 0 and makes the derivative NaN (z = 0 wherever a drag is inactive)
@inline function relaxation(::ExponentialRelaxation, z)
    small = z < sqrt(eps(z))
    return ifelse(small, one(z), -expm1(-z)) / ifelse(small, one(z), z)
end
