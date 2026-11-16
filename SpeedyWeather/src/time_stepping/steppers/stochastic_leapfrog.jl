export StochasticLeapfrog

"""Leapfrog time stepping with additive Wiener noise, i.e. the stochastic differential equation

```math
dx = F(x) \\, dt + a \\, dW
```

with the deterministic tendencies ``F`` of the model, white noise patterns of fixed amplitude
``a`` (see `AdditiveWienerNoise`) and the semi-implicit leapfrog scheme of `Leapfrog` (Euler
start, Robert-Asselin-Williams filter, implicit diffusion and gravity waves).

The scheme is drift-implicit: the noise is explicit and added to the tendencies as
``\\frac{a}{\\sqrt{\\Delta t}} \\xi`` with a new pattern ``\\xi`` (standard deviation
`standard_deviation`) every time step, right before the implicit terms, i.e. after the
parameterizations and the dynamical core. Schematically, with the linear terms ``L`` that are
treated implicitly,

```math
x^{n+1} = x^{n-1} + 2\\Delta t \\, (1 - 2\\Delta t \\, \\alpha L)^{-1}
\\left[ F(x^n, x^{n-1}) + \\frac{a}{\\sqrt{\\Delta t}} \\xi^n \\right]
```

(``\\alpha = 1`` for the horizontal diffusion, ``\\alpha`` = `centering` for the gravity waves of
the semi-implicit correction; for the first Euler step ``\\Delta t/2`` instead of ``2\\Delta t``).
The noise is on the right-hand side of the implicit solves like every explicit tendency, so the
implicit terms damp it within the time step, as the exact solution
``\\int e^{L(t-s)} a \\, dW(s)`` of the linear problem does. This is how stiff linear terms are
treated in the standard schemes for stochastic differential equations:

- the implicit Euler schemes of Kloeden and Platen (1992, chapter 12) and the stochastic
  theta method (Higham, 2000) evaluate (part of) the drift at the new time step and add the
  Wiener increment explicitly,
  ``x^{n+1} = x^n + \\Delta t [\\theta F(x^{n+1}) + (1-\\theta) F(x^n)] + a \\Delta W^n``,
- the stochastic semi-implicit leapfrog scheme of Ewald and Témam (2005) and its weak version
  (Ewald),
  ``x^{n+1} = x^{n-1} + 2\\Delta t F_1(x^n) + 2\\Delta t F_2(x^{n+1}) + a (\\Delta W^{n-1} + \\Delta W^n)``
  with the implicit part ``F_2``.

Here one increment ``2a\\sqrt{\\Delta t} \\, \\xi^n`` is added per leapfrog step instead of Ewald
and Témam's sum of two: the physical mode (the mean of the two interleaved leapfrog sequences)
receives the same ``a\\sqrt{\\Delta t} \\, (\\xi^{n-1} + \\xi^n)``, the computational mode that is
excited in addition is damped by the time filter. The alternative of adding the noise after the
implicit terms is a valid scheme too (split-step, Higham, Mao and Stuart, 2002) and converges to
the same equation for ``\\Delta t \\to 0``, but it leaves the noise in the stiff modes (the
smallest scales, fast gravity waves) undamped by the terms that are implicit because they are
stiff. Ocean and land variables have no implicit terms and are stepped with `EulerForward`, for
them the scheme is the Euler-Maruyama method.

The diffusion this scheme injects per variable is the nominal ``a^2 \\sigma^2 / 2`` (``\\sigma``
the `standard_deviation`) only where neither the implicit terms nor the time filter matter, and
should be measured where it matters (e.g. 1.33× nominal for a T8 barotropic model).

The pattern ``\\xi`` applied in a time step is drawn at the end of the previous time step and stored
in the variables (`vars.prognostic.wiener_pattern`, `vars.prognostic.wiener_pattern_surface`),
so the first time step of a new simulation is deterministic and a simulation continued with
another `run!` continues with the stored pattern. The noise realisation is reproducible
through `seed`, independent of other random processes.

`StochasticLeapfrog` is a `Leapfrog` whose `noise` field holds an `AbstractWienerNoise` (a type
alias, not a new time stepper): every method of `Leapfrog` applies, only the Wiener increment is
added. Use as `time_stepping = StochasticLeapfrog(spectral_grid; seed = 123)`, keyword arguments
other than those of `Leapfrog` are passed on to `AdditiveWienerNoise`, or equivalently
`Leapfrog(spectral_grid; noise = AdditiveWienerNoise(spectral_grid; seed = 123))`.
All amplitudes zero gives the same trajectory as `Leapfrog` without noise.

References

- Ewald, B. D. and R. Témam, 2005. Numerical analysis of stochastic schemes in geophysics.
  SIAM J. Numer. Anal. 42(6), 2257-2276, doi:10.1137/S0036142902418333
- Ewald, B. D. Weak versions of stochastic Adams-Bashforth and semi-implicit leapfrog schemes
  for SDEs (scheme 4.3), https://www.math.fsu.edu/~ewald/
- Higham, D. J., 2000. Mean-square and asymptotic stability of the stochastic theta method.
  SIAM J. Numer. Anal. 38(3), 753-769, doi:10.1137/S003614299834736X
- Higham, D. J., X. Mao and A. M. Stuart, 2002. Strong convergence of Euler-type methods for
  nonlinear stochastic differential equations. SIAM J. Numer. Anal. 40(3), 1041-1063,
  doi:10.1137/S0036142901389530
- Kloeden, P. E. and E. Platen, 1992. Numerical solution of stochastic differential equations.
  Springer, doi:10.1007/978-3-662-12616-5"""
const StochasticLeapfrog{NF, S, B, MS, O, L, N <: AbstractWienerNoise} = Leapfrog{NF, S, B, MS, O, L, N}

"""$(TYPEDSIGNATURES)
Leapfrog with additive Wiener noise. Keyword arguments of `Leapfrog` are passed on to it, all
others to `AdditiveWienerNoise` unless `noise` is provided."""
function StochasticLeapfrog(
        spectral_grid::SpectralGrid;
        Δt_at_T32 = Minute(40),
        adjust_with_output = true,
        robert_filter = 0.1,
        williams_filter = 0.53,
        ocean = EulerForward(spectral_grid; Δt_at_T32, adjust_with_output),
        land = EulerForward(spectral_grid; Δt_at_T32, adjust_with_output),
        noise = AdditiveWienerNoise(spectral_grid),
        kwargs...,
    )
    return Leapfrog(spectral_grid; Δt_at_T32, adjust_with_output, robert_filter, williams_filter, ocean, land, noise)
end

"""$(TYPEDSIGNATURES)
Add the Wiener increment of this time step (pattern drawn at the end of the previous one) to the
tendencies. Called in `time_step!` after the dynamical core and before the implicit terms of the
time stepping (`diffusion_and_implicit!`), which therefore act on the noise like on every explicit
tendency."""
function random_process!(vars::Variables, L::StochasticLeapfrog, model::AbstractModel)
    wiener_tendencies!(vars, L.noise, L, model)
    return nothing
end

"""$(TYPEDSIGNATURES)
Draw the noise patterns for the next time step, called in `time_step!` via `update_random!(vars, model)`."""
update_random!(vars::Variables, L::StochasticLeapfrog) = wiener_draw!(vars, L.noise)

"""$(TYPEDSIGNATURES)
Add the stored noise patterns times `amplitude/√Δt` to the tendencies: of the atmospheric
variables in spectral space, radius-scaled as the dynamical core leaves them, of the ocean and
land variables in grid space, see `AdditiveWienerNoise`."""
function wiener_tendencies!(
        vars::Variables,
        noise::AdditiveWienerNoise,
        L::AbstractTimeStepper,
        model::AbstractModel,
    )
    isdt = 1 / sqrt(L.Δt)               # Δt in [s], the Wiener increment scaling
    R = vars.prognostic.scale[]         # radius scaling
    wiener_atmosphere!(vars, noise, L, isdt, R)
    wiener_surface!(vars, noise, L, isdt)
    return nothing
end

# The drivers below are unrolled over the variable names, which are compile-time constants (they
# come from the type of the `Variables`), so that every `getfield` infers concretely; a loop over
# `tendency_names(vars)` would dispatch at runtime on every variable. Same pattern as
# `update_prognostic_atmosphere!` and `update_prognostic_namespaces!`.

#TODO: might update the function below to use the fused vairbales systems
"""$(TYPEDSIGNATURES)
Add the noise to the atmospheric tendencies. Slots of the spectral pattern: variables with layers
first, then single-layer variables, each in the order of `tendency_names`. After the dynamical core
the tendencies are scaled with the radius R (`scale_tendencies!`, the time step is Δt/R), vorticity
and divergence are stored scaled by R in addition (`scale_prognostic!`), so need R once more."""
@generated function wiener_atmosphere!(vars::Variables{Po, G, T}, noise, L, isdt, R) where {Po, G, T}
    calls = Expr[]
    for layers in (true, false), name in _tendency_names(T)
        # R: tendencies are radius-scaled (time step Δt/R); R²: vorticity, divergence are also stored × R
        a = name in (:vorticity, :divergence) ? :(amplitude(noise, $(QuoteNode(name))) * isdt * R * R) :
            :(amplitude(noise, $(QuoteNode(name))) * isdt * R)
        push!(
            calls, :(
                offset = wiener_slots!(
                    getfield(vars.tendencies, $(QuoteNode(name))), vars.prognostic.wiener_pattern,
                    $a, noise.taper, offset, Val($layers), L, noise,
                )
            )
        )
    end
    return Expr(:block, :(offset = 0), calls..., :(return nothing))
end

# add the noise to the tendency `var` if its layers match the pass `layers` (with layers first, then
# without), return the next slot of the pattern; `has_layers` of a concrete type is a constant
@inline function wiener_slots!(var, pattern, a, taper, offset, ::Val{layers}, L, noise) where {layers}
    tendency = get_tendency_step(var, L, noise)
    has_layers(tendency) == layers || return offset
    wiener!(tendency, pattern, a, layers ? taper : nothing, offset)
    return offset + nslots(tendency)
end

"""$(TYPEDSIGNATURES)
Add the noise to the ocean and land tendencies. Slots of the grid-point pattern only for variables
with a nonzero amplitude, not radius-scaled: ocean and land variables are stepped with the unscaled
time step. Without any such variable there is no `wiener_pattern_surface` and nothing to do."""
@generated function wiener_surface!(vars::Variables{Po, G, T}, noise, L, isdt) where {Po, G, T}
    :wiener_pattern_surface in fieldnames(Po) || return :(return nothing)
    calls = Expr[]
    for namespace in (:ocean, :land), name in _namespace_names(T, namespace)
        push!(
            calls, :(
                offset = wiener_surface_slots!(
                    getfield(getfield(vars.tendencies, $(QuoteNode(namespace))), $(QuoteNode(name))),
                    vars.prognostic.wiener_pattern_surface, amplitude(noise, $(QuoteNode(name))),
                    isdt, offset, time_stepper(L, Val($(QuoteNode(namespace)))), noise,
                )
            )
        )
    end
    return Expr(:block, :(offset = 0), calls..., :(return nothing))
end

# `time_stepping` is the ocean or land time stepper, e.g. EulerForward, which decides the tendency step
@inline function wiener_surface_slots!(var, pattern, a, isdt, offset, time_stepping, noise)
    iszero(a) && return offset
    tendency = get_tendency_step(var, time_stepping, noise)
    wiener!(tendency, pattern, a * isdt, nothing, offset)
    return offset + nslots(tendency)
end

# tendency[i, k] += a * taper[k] * pattern[i, offset + k] for every spectral or grid index i, layer k
function wiener!(tendency, pattern, a, taper, offset)
    iszero(a) && return nothing
    workorder = tendency isa LowerTriangularArray ? SpectralWorkOrder : RingGridWorkOrder
    launch!(architecture(tendency), workorder, size(tendency), wiener_kernel!, tendency, pattern, a, taper, offset)
    return nothing
end

@kernel inbounds = true function wiener_kernel!(tendency, pattern, a, taper, offset)
    I = @index(Global, NTuple)
    i = I[1]
    k = length(I) > 1 ? I[2] : 1
    t = isnothing(taper) ? one(a) : taper[k]
    tendency[I...] += a * t * pattern[i, offset + k]
end

"""$(TYPEDSIGNATURES)
Draw the noise patterns for the next time step."""
function wiener_draw!(vars::Variables, noise::AdditiveWienerNoise)
    update_random!(vars.prognostic.wiener_pattern, noise.spectral)
    haskey(vars.prognostic, :wiener_pattern_surface) &&
        update_random!(vars.prognostic.wiener_pattern_surface, noise.grid)
    return nothing
end
