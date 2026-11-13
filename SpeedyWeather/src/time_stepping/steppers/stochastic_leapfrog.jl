export StochasticLeapfrog, AdditiveWienerNoise

# IMPLEMENTATION NOTES (status 2026-11-13, work in progress, to be revised before the PR)
#
# Origin: ported from the FokkerPlanckSensitivities.jl project (SpectralWhiteNoiseProcess +
# GridWhiteNoiseProcess under a CombinedNoiseProcess driving an additive forcing). Up to commit
# 52c666f4 this file reproduced that code BITWISE (project: speedy/checks/upstream_wiener_equivalence.jl,
# `julia +1.11 --project=speedy-upstream`), with the noise added before the dynamical core. Since
# 2026-11-13 the noise is added right before the implicit terms instead: the same scheme, but
# results agree with the project's only to rounding (the gate compares one time step with a
# tolerance now). Tests in test/dynamics/stochastic_leapfrog.jl and
# test/parameterizations/{random_process,white_noise_process}.jl.
#
# Design
# - `StochasticLeapfrog` is a type alias of `Leapfrog` with an `AbstractWienerNoise` in its `noise`
#   field (no new stepper type). `random_process!(vars, time_stepping, model)` is called in both
#   `time_step!` methods between `dynamics_tendencies!` and `diffusion_and_implicit!` (no-op for
#   other steppers; explicit in `time_step!` on the user's request, 2026-11-13, rather than hidden
#   in a `diffusion_and_implicit!` method); `update_random!` advances random processes.
# - Amplitudes are a NamedTuple keyed by tendency name (any name the model time steps works,
#   e.g. η, soil_temperature); unknown non-default names warn in `initialize!`.
# - Variables are inferred from the model's tendencies. `tendency_names(vars)` & co. need the
#   `Variables`, which do not exist when `variables(noise, model)` is called, hence
#   `tendency_variables(model)` collects the same definitions from the model and its components
#   (excluding the time stepper to avoid recursion); at runtime `wiener_tendencies!` loops over
#   `tendency_names`, `ocean_tendency_names`, `land_tendency_names`.
# - Atmosphere: one spectral pattern, one slot per layer of every variable, variables with layers
#   first, then single-layer ones (tendency_names order is vor, div, temp, pres, humid, the project
#   had pres last). Zero amplitudes keep their slots (realisation independent of other amplitudes).
#   All atmospheric noise is added to the spectral tendencies, which at that point are complete
#   and radius-scaled (hence amplitudes × R, vorticity and divergence × R²).
# - Ocean/land: one grid-point pattern, slots ONLY for nonzero amplitudes (the project forced the
#   SST only; extra slots would shift its random stream). Inconsistent with the atmosphere rule,
#   a leftover of the bitwise reproduction, revisit.
# - The pattern is drawn at the end of a step and applied in the next (stored in the variables, so
#   the first step of a new simulation is deterministic) because the project's random process was
#   updated at the end of the step. Not needed otherwise: drawing and applying within the same step
#   would make the pattern scratch memory.
#
# Placement in the time step (reviewed 2026-11-12, decided 2026-11-13): right before the implicit
# terms, i.e. the noise is on the right-hand side of the implicit solves like an explicit tendency
# (drift-implicit / stochastic theta method, see the docstring of `StochasticLeapfrog` for the
# scheme and the references). Where in the time step the noise is added matters in two places only,
# both implicit parts of the time stepping that act on the whole tendency:
# 1. implicit horizontal diffusion, tendency = (tendency + expl*x) * impl, damps the noise within
#    the step. Scalar leapfrog+RAW test of an OU process dx = -νx dt + a dW (project:
#    speedy/checks/leapfrog_noise_placement.jl), stationary variance / exact a²/2ν, noise before
#    vs after the implicit terms: 2Δtν = 0.1: 1.28 vs 1.54; 0.34 (T32's smallest scale,
#    impl = 0.744): 1.39 vs 2.49; 3: 0.75 vs 12.0. Added after (e.g. in update_prognostic!, the
#    split-step variant), the small-scale noise is no longer controlled by the diffusion.
# 2. the semi-implicit gravity-wave correction (centering 0.5) divides the noise in fast gravity
#    wave modes by (1 - iαω2Δt): energy input ×1/(1+(αω2Δt)²), ≈0.07 for external waves at T32's
#    smallest scales (estimate). Added after, the noise fully excites modes the scheme slows down.
# The dynamical core in between is indifferent: it computes its terms from the state and only
# accumulates into the tendencies (linear in what is already there: radius scaling, transform,
# sign and global mean of the pressure tendency). The pattern has no l = 0 component, so the
# pressure noise conserves mass also without the dynamical core's mass fixer.
# Wiener scaling: a/√Δt in the tendency gives the nominal diffusion a²σ²/2 for leapfrog+RAW
# (scalar test: 1.00; a/√(2Δt) would give 0.50): each of the two interleaved leapfrog sequences
# gets 2a√Δt⋅ξ every 2Δt, their mean (the physical mode) a√Δt⋅(ξⁿ + ξⁿ⁺¹), the difference (the
# computational mode) is damped by the filter. Ewald's semi-implicit leapfrog (see docstring) adds
# ΔWⁿ⁻¹ + ΔWⁿ per leapfrog step instead, which does not excite the computational mode: an option.
# With `model.dynamics = false` (primitive models) `random_process!` is not called: no noise,
# as before the move.
#
# Open: `wiener_tendencies!` loops over runtime Symbols (type-unstable: ~21 KB, ~80 μs per step at
# T32 L8 on CPU, ~1% of a step), an @generated driver like `update_prognostic_atmosphere!` would
# fix it; GPU untested; full test suite not run; CHANGELOG to mention random_process! -> update_random!.

"""Supertype of the noise components of stochastic time steppers, see `StochasticLeapfrog`."""
abstract type AbstractWienerNoise <: AbstractModelComponent end

"""Additive Wiener noise for the prognostic variables, used by `StochasticLeapfrog`.
Every prognostic variable `x` and every layer is forced by its own independent white noise
pattern `ξ` (independent in time, between variables and between layers) with a fixed
amplitude `a`, so the noise is additive: its covariance depends neither on the state nor on
any parameter of the model. With time step `Δt` the tendency receives `a/√Δt ⋅ taper ⋅ ξ`,
so the state receives an increment ∝ `a√Δt ⋅ ξ` per step, a Wiener increment whose variance
is independent of the time step.

The amplitudes are given per variable by the name of its tendency, e.g.
`AdditiveWienerNoise(spectral_grid; vorticity = 3e-8, η = 1e-3, soil_temperature = 1e-4)`,
variables not listed get no noise. Units are those of the variable per √s, e.g. 1/s^(3/2) for
vorticity and divergence, K/s^(1/2) for temperature. Which variables there are is inferred from
the tendencies of the model (`tendency_names`, `ocean_tendency_names`, `land_tendency_names`):

- atmospheric variables are forced with one `SpectralWhiteNoiseProcess` (spectrally flat by
  default) that has a pattern slot for every layer of every variable (variables with layers first,
  then single-layer variables, each in the order of `tendency_names`), also for zero amplitudes,
  so that switching one amplitude on or off does not change the noise realisation of the others.
  The noise is added to the spectral tendencies, so its covariance is diagonal in the spectral
  coefficients. The pattern has no global mean, the noise therefore conserves e.g. the mass
  (surface pressure). Variables with layers are tapered with `tapering(σ)`.
- ocean and land variables (e.g. sea surface temperature) with a nonzero amplitude are forced by a
  `GridWhiteNoiseProcess`, independently at every grid point and layer.

The amplitudes multiply a pattern of standard deviation `standard_deviation` (1/3 by default,
as for `SpectralAR1Process`). The defaults are the amplitudes used for T32 L8 PrimitiveWet
models in a 1500-year stochastic control run: climate-neutral, stable, humidity bounded by
clipping rather than by the climate (see the docstring of `StochasticLeapfrog`). They are not
universal, check them for other resolutions, time steps or models. The noise is drawn on
the model's architecture (see `architecture_rng`); `draw_on_cpu = true` draws it on the CPU also
on GPU (slower, but the same numbers as on CPU).

Fields are $(TYPEDFIELDS)"""
struct AdditiveWienerNoise{NF, A <: NamedTuple, F, VectorType, SP, GP} <: AbstractWienerNoise
    "[OPTION] Noise amplitude per variable, by the name of its tendency, zero for all others"
    amplitudes::A

    "[OPTION] Vertical tapering as a function of sigma, for atmospheric variables with layers"
    tapering::F

    "[DERIVED] Vertical tapering per layer"
    taper::VectorType

    "[OPTION] White noise process in spectral space for the atmospheric variables"
    spectral::SP

    "[OPTION] White noise process in grid space for the ocean and land variables"
    grid::GP
end

"""Default amplitudes of `AdditiveWienerNoise`, for T32 L8 PrimitiveWet models."""
const DEFAULT_WIENER_AMPLITUDES = (
    vorticity = 3.0e-8,                 # [1/s^(3/2)]
    divergence = 3.0e-8,                # [1/s^(3/2)]
    temperature = 3.0e-3,               # [K/s^(1/2)]
    humidity = 3.0e-8,                  # [kg/kg/s^(1/2)]
    sea_surface_temperature = 1.0e-3,   # [K/s^(1/2)]
)

"""$(TypedMethodSignatures(defaults = true))
Additive Wiener noise for the prognostic variables, all keyword arguments other than the ones
listed are noise amplitudes by variable name, added to (or replacing) `DEFAULT_WIENER_AMPLITUDES`.
`seed` seeds both white noise processes (they draw from independent streams), `wavenumber` and
`standard_deviation` are passed on to the `SpectralWhiteNoiseProcess`, `draw_on_cpu` to both."""
function AdditiveWienerNoise(
        SG::SpectralGrid;
        tapering = σ -> 1,
        seed = 0,
        wavenumber = 0,
        standard_deviation = 1 / 3,
        draw_on_cpu = false,
        amplitudes...,
    )
    NF = SG.NF
    amplitudes = map(NF, merge(DEFAULT_WIENER_AMPLITUDES, amplitudes))
    taper = on_architecture(SG.architecture, zeros(NF, SG.nlayers))
    spectral = SpectralWhiteNoiseProcess(SG; seed, wavenumber, standard_deviation, draw_on_cpu)
    grid = GridWhiteNoiseProcess(SG; seed, standard_deviation, draw_on_cpu)
    return AdditiveWienerNoise{NF, typeof(amplitudes), typeof(tapering), typeof(taper), typeof(spectral), typeof(grid)}(
        amplitudes, tapering, taper, spectral, grid,
    )
end

"""$(TYPEDSIGNATURES) Noise amplitude of the variable `name`, zero if not given."""
amplitude(noise::AdditiveWienerNoise{NF}, name::Symbol) where {NF} = get(noise.amplitudes, name, zero(NF))

"""$(TYPEDSIGNATURES)
The tendency variables of `model`, i.e. the definitions behind `tendency_names` and co., as they
are needed before the `Variables` exist. Collected from the model and its components except the
time stepper (its variables include the noise patterns, which are defined from these)."""
function tendency_variables(model::AbstractModel)
    components = (getproperty(model, name) for name in propertynames(model) if name != :time_stepping)
    all = (variables(model)..., Iterators.flatten(variables(c, model) for c in components)...)
    return unique(v -> (v.namespace, v.name), filter(v -> v isa TendencyVariable, collect(all)))
end

# number of pattern slots of a tendency (of one time step): its atmosphere, land or ocean layers, else 1
nslots(tendency) = size(tendency.data, 2)
nslots(v::TendencyVariable, model::AbstractModel) = nslots(get_step(allocate(v, model), 1))
has_layers(tendency) = ndims(tendency.data) == 2

function variables(noise::AdditiveWienerNoise, model::AbstractModel)
    tendencies = tendency_variables(model)
    atmosphere = filter(v -> v.namespace == Symbol(), tendencies)
    surface = filter(v -> v.namespace in (:ocean, :land) && !iszero(amplitude(noise, v.name)), tendencies)
    n = sum(v -> nslots(v, model), atmosphere)

    spectral = PrognosticVariable(
        :wiener_pattern, Spectral3D(; n), units = "1",
        desc = "Pattern of the Wiener increment of the next time step, one slot per variable and layer",
    )
    surface = isempty(surface) ? () : (
            PrognosticVariable(
                :wiener_pattern_surface, Grid3D(; n = sum(v -> nslots(v, model), surface)), units = "1",
                desc = "Grid-point pattern of the Wiener increment of ocean and land variables of the next time step",
            ),
        )
    return (spectral, surface...)
end

function initialize!(noise::AdditiveWienerNoise, model::AbstractModel)
    initialize!(noise.spectral, model)
    initialize!(noise.grid, model)

    # with the same seed both processes would draw the same normals, derive a distinct
    # (but reproducible) stream for the grid process
    if noise.spectral.seed == noise.grid.seed && noise.grid.seed != 0
        Random.seed!(noise.grid.random_number_generator, hash(noise.grid.seed, UInt(0x9e3779b97f4a7c15)))
    end

    # vertical tapering, in a kernel as a comprehension would scalar-index the vertical coordinates on GPU
    coord = model.geometry.vertical_coordinates
    launch!(
        architecture(noise.taper), LinearWorkOrder, (get_nlayers(coord),),
        wiener_taper_kernel!, noise.taper, noise.tapering, coord,
    )

    # amplitudes for variables this model does not time step are ignored, warn unless a default
    names = [v.name for v in tendency_variables(model) if v.namespace in (Symbol(), :ocean, :land)]
    ignored = filter(name -> !(name in names || haskey(DEFAULT_WIENER_AMPLITUDES, name)), keys(noise.amplitudes))
    isempty(ignored) || @warn "AdditiveWienerNoise: the model has no tendency for $ignored, amplitude ignored."
    return nothing
end

@kernel inbounds = true function wiener_taper_kernel!(taper, tapering, coordinate)
    k = @index(Global, Linear)
    taper[k] = tapering(sigma(k, coordinate))
end

"""Leapfrog time stepping with additive Wiener noise, i.e. the stochastic differential equation

    dx = F(x) dt + a dW

with the deterministic tendencies F of the model, white noise patterns of fixed amplitude `a`
(see `AdditiveWienerNoise`) and the semi-implicit leapfrog scheme of `Leapfrog` (Euler start,
Robert-Asselin-Williams filter, implicit diffusion and gravity waves).

The scheme is drift-implicit: the noise is explicit and added to the tendencies as `a/√Δt ⋅ ξ`
with a new pattern `ξ` (standard deviation `standard_deviation`) every time step, right before
the implicit terms, i.e. after the parameterizations and the dynamical core. Schematically,
with the linear terms `L` that are treated implicitly,

    xⁿ⁺¹ = xⁿ⁻¹ + 2Δt (1 - 2Δt α L)⁻¹ [F(xⁿ, xⁿ⁻¹) + a/√Δt ⋅ ξⁿ]

(`α = 1` for the horizontal diffusion, `α = centering` for the gravity waves of the semi-implicit
correction; for the first Euler step `Δt/2` instead of `2Δt`). The noise is on the right-hand
side of the implicit solves like every explicit tendency, so the implicit terms damp it within
the time step, as the exact solution `∫exp(L(t-s)) a dW(s)` of the linear problem does. This
is how stiff linear terms are treated in the standard schemes for stochastic differential
equations:

- the implicit Euler schemes of Kloeden and Platen (1992, chapter 12) and the stochastic
  theta method (Higham, 2000) evaluate (part of) the drift at the new time step and add the
  Wiener increment explicitly, `xⁿ⁺¹ = xⁿ + Δt[θ F(xⁿ⁺¹) + (1-θ) F(xⁿ)] + a ΔWⁿ`,
- the stochastic semi-implicit leapfrog scheme of Ewald and Témam (2005) and its weak version
  (Ewald), `xⁿ⁺¹ = xⁿ⁻¹ + 2Δt F₁(xⁿ) + 2Δt F₂(xⁿ⁺¹) + a (ΔWⁿ⁻¹ + ΔWⁿ)` with the implicit part `F₂`.

Here one increment `2a√Δt ⋅ ξⁿ` is added per leapfrog step instead of Ewald and Témam's sum of
two: the physical mode (the mean of the two interleaved leapfrog sequences) receives the
same `a√Δt (ξⁿ⁻¹ + ξⁿ)`, the computational mode that is excited in addition is damped by the
time filter. The alternative of adding the noise after the implicit terms is a valid scheme
too (split-step, Higham, Mao and Stuart, 2002) and converges to the same equation for `Δt → 0`,
but it leaves the noise in the stiff modes (the smallest scales, fast gravity waves) undamped
by the terms that are implicit because they are stiff. Ocean and land variables have no
implicit terms and are stepped with `EulerForward`, for them the scheme is the
Euler-Maruyama method.

The diffusion this scheme injects per variable is the nominal `a²σ²/2` (`σ` the
`standard_deviation`) only where neither the implicit terms nor the time filter matter, and
should be measured where it matters (e.g. 1.33× nominal for a T8 barotropic model).

The pattern `ξ` applied in a time step is drawn at the end of the previous time step and stored
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

"""$(TypedMethodSignatures(defaults = true))
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
        noise = nothing,
        kwargs...,
    )
    if isnothing(noise)
        noise = AdditiveWienerNoise(spectral_grid; kwargs...)
    elseif !isempty(kwargs)
        throw(ArgumentError("Provide either `noise` or its keyword arguments $(keys(kwargs)), not both."))
    end
    return Leapfrog(spectral_grid; Δt_at_T32, adjust_with_output, robert_filter, williams_filter, ocean, land, noise)
end

"""$(TYPEDSIGNATURES)
Add the Wiener increment of this time step (drawn at the end of the previous one) to the
tendencies, then draw the pattern for the next time step. Called in `time_step!` after the
dynamical core and before the implicit terms of the time stepping (`diffusion_and_implicit!`),
which therefore act on the noise like on every explicit tendency."""
function random_process!(vars::Variables, L::StochasticLeapfrog, model::AbstractModel)
    wiener_tendencies!(vars, L.noise, L, model)
    wiener_draw!(vars, L.noise)
    return nothing
end

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

    # ATMOSPHERE, slots of the spectral pattern: variables with layers first, then single-layer ones
    offset = 0
    for layers in (true, false), name in tendency_names(vars)
        tendency = get_tendency_step(getfield(vars.tendencies, name), L, noise)
        has_layers(tendency) == layers || continue

        # after the dynamical core the tendencies are scaled with the radius R (scale_tendencies!,
        # the time step is Δt/R), vorticity and divergence are stored scaled by R in addition
        # (scale_prognostic!), so need R once more
        a = amplitude(noise, name) * isdt * R * (name in (:vorticity, :divergence) ? R : one(R))
        wiener!(tendency, vars.prognostic.wiener_pattern, a, layers ? noise.taper : nothing, offset)
        offset += nslots(tendency)
    end

    # OCEAN, LAND, slots of the grid-point pattern for variables with a nonzero amplitude, not
    # radius-scaled: ocean and land variables are stepped with the unscaled time step
    offset = 0
    for (namespace, names) in ((:ocean, ocean_tendency_names(vars)), (:land, land_tendency_names(vars))), name in names
        a = amplitude(noise, name)
        iszero(a) && continue
        tendency = get_tendency_step(getfield(getfield(vars.tendencies, namespace), name), L, noise)
        wiener!(tendency, vars.prognostic.wiener_pattern_surface, a * isdt, nothing, offset)
        offset += nslots(tendency)
    end
    return nothing
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
