"""Abstract supertype for the time integration of fluxes or drags that relax a variable towards
a target value, e.g. the surface fluxes `flux = -conductance * (x - x_surface)` or the
`SpeedLimitDrag` relaxing the wind speed towards the speed limit. Subtypes define
`flux_damping(scheme, z)` with `z = Δt/τ` the ratio of time step to relaxation time scale,
a factor that multiplies the explicit flux (or tendency)."""
abstract type AbstractFluxIntegration end
Base.broadcastable(scheme::AbstractFluxIntegration) = Ref(scheme)

export ExplicitFluxes, ImplicitFluxes, ExponentialFluxes, RationalFluxes

"""Explicit (forward) time integration of surface fluxes, unstable for `Δt/τ > 2`."""
struct ExplicitFluxes <: AbstractFluxIntegration end

"""Implicit (backward Euler) time integration of surface fluxes, damping the explicit
flux by `1/(1 + Δt/τ)`. Unconditionally stable."""
struct ImplicitFluxes <: AbstractFluxIntegration end

"""Exponential time differencing of surface fluxes, damping the explicit flux by
`(1 - exp(-Δt/τ))/(Δt/τ)`. Exact for a fixed surface value, unconditionally stable."""
struct ExponentialFluxes <: AbstractFluxIntegration end

"""Rational approximation of `ExponentialFluxes`, damping the explicit flux by
`(6 + z)/(6 + 4z + z²)` with `z = Δt/τ`. Matches `(1 - exp(-z))/z` to second order
in `z` and as `1/z` for large `z`, at the cost of `ImplicitFluxes`. Unconditionally
stable but overshoots the surface value slightly (≤ 10%) for `z > 3`."""
struct RationalFluxes <: AbstractFluxIntegration end

const DEFAULT_FLUX_INTEGRATION = RationalFluxes

"""$(TYPEDSIGNATURES)
Factor φ(z) multiplying an explicit relaxation flux for `z = Δt/τ` given time integration `scheme`."""
@inline flux_damping(::ExplicitFluxes, z) = one(z)
@inline flux_damping(::ImplicitFluxes, z) = 1 / (1 + z)
@inline flux_damping(::RationalFluxes, z) = (6 + z) / (6 + 4z + z^2)

# numerator and denominator guarded separately (1/1 for small z), not with ifelse(small, 1, -expm1(-z)/z),
# whose untaken branch is 0/0 for z = 0 and makes the derivative NaN (z = 0 wherever a drag is inactive)
@inline function flux_damping(::ExponentialFluxes, z)
    small = z < sqrt(eps(z))
    return ifelse(small, one(z), -expm1(-z)) / ifelse(small, one(z), z)
end
