# NCycleLorenz for the shallow water model: bug fixes and implicit centering

> Status: **in progress**. Bugs fixed, implicit solver moved into the time stepper, per-substep implicit centering with an order condition (α_eff), 4-cycle variant A as new default; AB, ABBA left for future work.

Date of initial draft: 2026-11-12

Base revision: 069fb09f3d4fa539fbeb27ea4c352010d6737961

## Originating prompt

> Can you revisit the use of NCycleLorenz for the shallow water model? It's implemented and
> running but not default because when I last tested it required a shorter time step and it was
> more diffusive too. Have a look at the PRs that implemented this and check whether this looks
> all correct and matches the docs, particularly the implicit solver (play with the centering
> too). Then take the gravity wave example from the docs and compare to Leapfrog, I would like the
> NCycleLorenz to preserve the wave energy (variance of eta) well but also not require a (much)
> shorter time step than Leapfrog. Can you find a setup/formulation of NCycleLorenz that solves
> this example as well as Leapfrog?

## Revision log

> Can you open a PR with this? Bundle all your findings together! Don't bother making N=4 even
> better for now I want to have a look first

> Can you clarify how your changes compare to the Daisuke's version and Milan's version that I
> wrote down in the originating issue? / Yes map to the issue and document the exact mathematical
> formulation of the implicit too

- Mapped all formulations to the names in #917 and wrote down the exact implicit solve
  (here and in `docs/src/shallowwater.md`).
- Plan written after the investigation (exploratory work), bundling bug fixes, the analysis and
  the per-substep centering for N=3. N=4 centerings are only reported, not made default.

> [Review comment on #1402 on variant A vs B, then:] With this per-substep centering, is there
> a preference for variant A or B? / yes please switch to the more robust A centering but can you
> also implement the B centering so it's at least default when `variant = NCycleLorenzB()`

- The A centering (0.069, 1.028, 0.48) was from an optimisation without the orography
  constraint and is weakly unstable over orography for ω_LΔt ≳ 5 (1.0005 at β=0.1, ω_LΔt=5;
  1.03 at ω_LΔt=10). Replaced by the robust optimum (0.116, 1.039, 0.462).
- Added the robust B optimum (1.386, 0.468, -0.464) as default for `NCycleLorenzB`. Two decimals
  are not enough (rounded to (1.39, 0.47, -0.46) it amplifies with 1.0014).
- In the model B turned out more stable than A (table below), although linearly both have the
  same explicit limit.

> Can you update the docs example with the gravity waves into a 1km high Gaussian blob in the
> Pacific somewhere, maybe 5˚ large and run the simulation as long as it takes for the wave to be
> roughly 1/4 around the globe so it's like a ring around the planet / Can you rechoose the A/B
> centerings so that it's more stable and relax the gravity wave energy constraint somewhat. It's
> okay for NCycleLorenz to not preserve wave energy as much as Leapfrog but it being
> unpredictably unstable is a no go.

- Docs example: Gaussian blob, 5000 m (1000 m first, but then the ring is not visible in the
  layer thickness h), σ = 5˚ at 160˚W, 10˚N, 10 hours, plotting η and h.
- Centerings re-optimised for stability margin (see "Robust centerings" below).

> In the mean time, can you reorganize model.time_stepping and model.implicit. The latter is
> dependent on the former but currently exist as independent model components, I suggest to
> arrange them as model.time_stepping.implicit for all time steppers and then set implicit =
> nothing for Barotropic and the implicit solvers otherwise. Adapt the current function calls and
> barriers to this new order and label the PR as breaking

- Implicit solver moved into the time steppers (see "Implicit solver as part of the time stepper").

> Are the variables that ImplicitPrimitiveEquation requires actually still used? / Can you just
> remove the average temperature profile calculation? When the GlobalMeanReferenceTemperature
> needs it it could just calculate at that step and directly write into the vector in implicit?
> The callback is just an example, calculate the global surface temperature differently

- Removed `vars.dynamics.average_temperature_profile` and `temperature_average!` (computed every
  step in `transform!`). `GlobalMeanReferenceTemperature` computes the global-mean profile from the
  l=m=0 harmonic of the spectral temperature when the reference profile is recomputed
  (`global_mean_temperature`), `GlobalSurfaceTemperatureCallback` likewise for the lowest layer.

> Merge main in yes and try to find centerings for N=4 I am keen to know whether it's generally
> more/less stable or more/less diffusive. Can you also think whether there is a way to formalise
> the centerings mathematically. Like do you want them to be like Crank Nicolson on average or what
> do negative centerings actually mean. Intuitively only 0-1 makes sense to me. Is there something
> better than just trial and error for these? / Yes can you run a few example across resolutions
> with N=4 new A centering and check that it's always stable? / Yes please!! / And also can you
> write up some details about the mathematics of the centerings we discussed into the docs?

- Merged main. Derived the order condition / effective centering α_eff (see "Order condition
  and effective centering"), chose centerings by fixing α_eff and maximising stability.
- New defaults: `steps = 4` (was 3), 4-cycle A θ = (2.06, 0.93, 1.02, 0.09), 3-cycle B
  θ = (1.63, 0.55, -0.55); 3-cycle A unchanged (-1.03, 0.70, 1.05). `effective_centering`
  function. Mathematics in `docs/src/time_stepping.md` ("Implicit centering of the Lorenz N-cycle").

## Problem description

`NCycleLorenz` (#1031, #1035) is not the default for the shallow water model as it was found to
be more diffusive and to need a shorter time step than `Leapfrog`.

## Background

### Mapping to the formulations in #917

Originating issue: [#917](https://github.com/SpeedyWeather/SpeedyWeather.jl/issues/917). Notation per
substep k of the cycle, state ``u_i``, weight ``w = w_k``, explicit tendency ``F_E``, linear
gravity-wave terms ``L``, full tendency ``F = F_E + Lu``, accumulated tendency ``G``.

| Name in #917 | Formulation | Status |
|---|---|---|
| 0. Semi-implicit leapfrog | ``(u_{i+1}-u_{i-1})/2Δt = (1-2αΔtL)^{-1}[F(u_i) - L(u_i - u_{i-1})]`` | unchanged |
| 2. Daisuke's version (Hotta et al. 2016) | ``(u_{i+1}-u_i)/Δt = (1-αΔtL)^{-1}[wF_E(u_i) + (1-w)G + Lu_i]``, ``G \gets wF_E + (1-w)G`` | never on main, was described in the docs; unstable (Doppler shift, Coriolis) |
| 3. Milan's version | ``G \gets (1-wαΔtL)^{-1}[wF(u_i) + (1-w)G]``, ``u_{i+1} = u_i + ΔtG`` | on main; stable but damps gravity waves |
| **This PR** | ``G \gets (1-wθ_kΔtL)^{-1}[wF(u_i) + (1-w)G]``, ``u_{i+1} = u_i + ΔtG`` | Milan's version with centering ``θ_k`` per substep; ``θ_k = α`` recovers it exactly |

Equivalently, this PR evaluates the linear terms inside the Lorenz average at
``u_i + θ_k(u_{i+1} - u_i)``:
``G \gets w(F_E(u_i) + L(u_i + θ_k(u_{i+1}-u_i))) + (1-w)G``.

#### Exact implicit solve (shallow water)

``L``: ``∂D/∂t = -g∇²η``, ``∂η/∂t = -HD`` (D divergence, η interface displacement, H layer
thickness at rest, ∇² the eigenvalue ``-l(l+1)`` per harmonic of the radius-scaled equations).
With ``ξ = wθ_kΔt`` (divided by the radius scaling), per harmonic and substep:

```
R_D = w F_D + (1-w) G_D
R_η = w F_η + (1-w) G_η
G_D ← (R_D - ξ g ∇² R_η) / (1 - ξ² g H ∇²)
G_η ← R_η - ξ H G_D
```

with ``G_ζ ← w F_ζ + (1-w) G_ζ`` for vorticity. `implicit_correction!` stores
``F ← (G_new - (1-w)G_old)/w`` so that the shared `update_prognostic!`
(``G ← wF + (1-w)G``, ``u ← u + ΔtG``) gives the same result for all variables. Horizontal
diffusion ``F ← (F + ν∇²ⁿu)/(1 - Δtν∇²ⁿ)`` is applied before the implicit correction.

Daisuke's version as tested here: ``G ← wF_E + (1-w)G`` with
``F_E = F - Lu``, ``u ← u + Δt(1-αΔtL)^{-1}(G + Lu)``. The Lawson variant additionally
applies ``R = (1-αΔtL)^{-1}(1+(1-α)ΔtL)`` to state and memory, ``u ← R(u + ΔtG)``, ``G ← RG``.

### Bugs found

1. `ShallowWaterModel` + `NCycleLorenz` crashed on main (`BoundsError` in `implicit_correction!`):
   the spectral η tendency was allocated with the prognostic step count `ps` instead of the
   tendency step count `ts` (regression from #1150, `shallow_water.jl`).
2. After fixing 1., all shallow water tendencies were zero: the fused spectral intermediates
   (`kinetic_energy`, `uh`, `vh`, and `uT_anomaly`, `vT_anomaly`, `uq`, `vq` in the primitive
   models) share the step dimension of the fused tendency parent, i.e. 2 steps for
   `NCycleLorenz`. The batched transform writes into step 1 (`get_tendency_step`), but they were
   read with `get_step(var)`, which returns the *last* step (still zero). Now read with
   `get_tendency_step(var, time_stepping, DynamicalCore())`, as `u_tendency`, `v_tendency` already were.
3. Docs and code differ: the docs described the linear terms outside of the Lorenz average
   (`dx = (I - αΔtL)^-1 (G + Lx)` with `G` accumulating explicit tendencies only), the code puts
   them inside (`G = (I - wαΔtL)^-1 (wF + (1-w)G)`). Docs updated to the code.

### Linear stability analysis

Oscillation equation `x' = iω_E x + iω_L x` (explicit + implicit part), per-step amplification
over a full cycle (see "Methods" below to reproduce). Test cases: pure gravity wave
(ω_E = 0), Doppler shift (ω_E > 0), f-plane (Coriolis explicit, gravity waves implicit), and
"orography" ω_E = -β ω_L, β ∈ (0, 0.9): the reference layer thickness H used for L is larger
than the actual one (H - H_b), which semi-implicit leapfrog handles unconditionally.

| Scheme | Gravity wave damping/step at ωΔt=0.5, 1 | Doppler | f-plane | orography β |
|---|---|---|---|---|
| Leapfrog + RAW, α=0.5 | 0.9994, 0.997 | stable to ω_EΔt≈1 | stable | stable |
| Milan's version (on main), θ=0.5 (N=3) | 0.90, 0.75 | unstable (1.07) | stable | stable |
| Daisuke's version, α=0.5 | 1, 1 | unstable (1.2-1.5) | unstable | - |
| Lawson / integrating factor | 1, 1 | stable | unstable (1.04-1.5) | - |
| Milan's version, per-substep θ = (1, 1, 0), N=3 | exactly 1 (all-pass) | stable to 0.55 | stable | unstable (1.1-2.6) |
| Milan's version, θ = (0.069, 1.028, 0.48), A (first version of this PR) | 1, 0.9997 | stable to ≈0.3 | stable | unstable for ω_LΔt ≳ 5 (1.03 at 10) |
| Milan's version, θ = (0.116, 1.039, 0.462), A (second version of this PR) | ≥ 0.997 for ωΔt ≤ 1.5 | stable to 0.4 | stable | stable |
| Milan's version, θ = (1.386, 0.468, -0.464), B (second version of this PR) | ≥ 0.998 for ωΔt ≤ 1.5 | stable to 0.4 | stable | stable |
| Milan's version, θ = 0.5, B (main) | ≥ 0.72 | unstable (1.12) | stable | stable |

Variant A vs B: the weights are the same in a different order, A = (1, 3/2, 3), B = (1, 3, 3/2),
so the centering has to be optimised per variant (A's centering applied to B is unstable). The
explicit amplification of A and B is identical for linear constant-coefficient problems (limit
ω_EΔt < 0.57 for both). The A optimum is well defined (all optimiser runs converge to it), the B
optimum less so (other runs land at e.g. (10.6, -5.1, 10.3)) and has a negative centering on the
last substep. Phase speed at ωΔt = 0.6 relative to exact: A 0.84, B 0.82.

The explicit Lorenz N-cycle itself is limited to ω_EΔt < √3/3 ≈ 0.58 (N=3) and ≈ 0.71 (N=4)
per tendency evaluation, compared to 1 for leapfrog. The robust optimum for N=4 found so far,
θ ≈ (2.31, 0.86, 0.86, 0.08), damps gravity waves more (≈0.994/step) and is not tested in the model.

The hyperdiffusion is applied implicitly with Δt before the Lorenz weighting; for N=3 the
weight w=3 makes stiff diffusion unstable (|A| ≈ 1.2 for νΔt ≳ 10). Using wΔt makes it stable,
but in the experiments this was not the limiting factor (not changed here).

### Experiments: docs gravity wave example (T128, 2 days, var(η)/var(η₀))

| Δt at T32 | Leapfrog | main (θ=0.5) | A, θ = (0.069, 1.028, 0.48) | A, θ = (0.116, 1.039, 0.462) | B, θ = (1.386, 0.468, -0.464) |
|---|---|---|---|---|---|
| 10 min | 0.36 | 0.011 | 0.43 | 0.43 | 0.39 |
| 15 min | | | | 0.41 | 0.36 |
| 20 min | 0.38 | 0.011 | 0.31 | 0.38 | 0.32 |
| 25 min | | | | 0.24 | 0.25 |
| 30 min | 0.35 | 0.011 | 0.30 | unstable | 0.25 |
| 40 min | 0.29 | 0.011 | unstable | | 0.19 |
| 50 min | 0.23 | 0.011 | | | 0.17 |
| 60 min | 0.18 | 0.011 | | | 0.15 |

Galewsky jet (`ZonalJet`, 20 days, T85), largest stable Δt_T32 tested:
leapfrog 50 min (60 unstable), B 40 min (45 unstable), A 25 min (30 unstable), main 25 min,
A with (0.069, 1.028, 0.48) 20 min. Daisuke's version and the Lawson variant were unstable
beyond Δt_T32 ≈ 10-20 min in the gravity wave example.

Why B is more stable than A in the model while linearly equal is not understood yet.

### Methods

Analysis scripts are not committed; this section describes them so they can be recreated.

**Linear stability analysis.** Units Δt = 1. Explicit part ``E`` and implicit part ``L`` are
matrices (scalars for the 1-variable cases), ``I`` identity. The state over a substep is the pair
``(x, G)``, so each substep is a linear map on the stacked vector ``[x; G]``; the cycle matrix ``T``
is the product of the N substep maps (in order k = 0..N-1), and the amplification per step is
``ρ(T)^{1/N}`` (spectral radius). Since ``w_0 = 1`` erases ``G`` at the start of a cycle, for a pure
gravity wave the cycle reduces to a scalar map ``x ↦ R(iω_L)x`` (the ``[1,1]`` element of ``T``).

- Milan's version with centering ``θ_k``, weight ``w``: ``S = (I - wθ_kL)^{-1}``,
  ``δ = S(w(E+L)x + (1-w)G)``, ``x ← x + δ``, ``G ← δ``, i.e. substep matrix
  ``[[I + Sw(E+L), S(1-w)], [Sw(E+L), S(1-w)]]``.
- Daisuke's version: ``S = (I - αL)^{-1}``, ``G ← wEx + (1-w)G``, ``x ← x + S(G + Lx)``.
- Lawson variant: ``R = S(I + (1-α)L)``, ``G̃ = wEx + (1-w)G``, ``x ← R(x + G̃)``, ``G ← RG̃``.
- Leapfrog reference: state ``(x_{n-1}, x_n)``,
  ``(I - 2αL)x_{n+1} = (I + 2(1-α)L)x_{n-1} + 2Ex_n``, then the RAW filter with
  ``d = ν/2(x_{n-1} - 2x_n + x_{n+1})``, ``x_n ← x_n + α_W d``, ``x_{n+1} ← x_{n+1} - (1-α_W)d``,
  ``ν = 0.1``, ``α_W = 0.53``; for the "physical mode" take the eigenvalue closest to ``exp(i(ω_E+ω_L))``.
- Weights: A ``w_k = N/(N-k)``, B ``w_k = N/k``, ``w_0 = 1``.

Test cases (each maximised over ``ω_LΔt`` on a grid of 30 points in [0, 3] and 40 in [3, 40]):
pure gravity wave ``E = 0, L = iω_L``; Doppler ``E = iω_E``, ``ω_EΔt ∈ {0.1, 0.2, 0.3, 0.4}``;
f-plane with ``x = (ζ, D, η)``, ``E = [[0, -f, 0], [f, 0, 0], [0, 0, 0]]``,
``L = [[0, 0, 0], [0, 0, ω_L], [0, -ω_L, 0]]``, ``fΔt ∈ {0.1, 0.35}``; orography ``E = -βL``,
``β ∈ {0.1, 0.2, 0.3, 0.5, 0.7, 0.9}``. Explicit limit: largest ``ω_EΔt`` with ``L = 0``, ``θ = 0`` and
amplification ≤ 1.

Hyperdiffusion: scalar ``x' = -νx`` applied as in the code, tendency ``-νx/(1 + νΔt·c)`` then Lorenz
weighting, with ``c = 1`` (current) or ``c = w`` (weighted).

**Optimisation of the centering.** Objective: maximise the minimum gravity-wave amplification
per step over ``ω_LΔt ∈ [0.05, 1.5]`` (15 points), subject to the maximum amplification over all
test cases (incl. pure gravity waves up to ``ω_LΔt = 40``) being ≤ 1 + 1e-6, implemented as a
penalty of 1000 × excess. Random search: 28 starts ``θ ~ 0.6 + 0.8·N(0, 1)`` (componentwise),
2000 iterations of Gaussian perturbations with step 0.4 halved every 400 iterations, accepting
improvements only. Rank the starts by objective and check that the best ones agree.
The A optimum was reproduced by all top starts; for B the top starts disagree.
Afterwards round and re-check (two decimals are not enough, three are).

**Model experiments.** Experimental formulations were implemented by defining a new
`AbstractImplicit` subtype and its `implicit_correction!(vars, implicit, ::AbstractNCycleLorenz,
::ShallowWater)` in a script, modifying the prognostic variables in place where needed.
- Gravity waves: the docs example (T128, 1 layer, `RandomWaves(lmin=10, lmax=30)` with amplitude
  2000 m, `EarthOrography(smoothing=false)`, `ImplicitShallowWater(centering=0.5)`), 2 days,
  `adjust_with_output=false`, `Random.seed!(1234)`. Metric: ``Σ_{l>0} c_m |η_{lm}|²`` with
  ``c_m = 1`` for m = 0 and 2 otherwise, relative to the initial value, from a callback; NaN or
  growth > 10 counts as unstable.
- Galewsky jet: `ShallowWaterModel` with `ZonalJet` initial conditions and default orography, T85,
  20 days; stable = no NaNs detected.

### Order condition and effective centering

Method: symbolic expansion (sympy) of one cycle of the oscillation equation
``x' = (λ_E + λ_I)x`` with ``a = λ_EΔt``, ``b = λ_IΔt``, per substep
``G ← (w(a+b)x + (1-w)G)/(1 - wθ_k b)``, ``x ← x + G``, starting from ``G = 0`` (``w_0 = 1``), in powers
of ``a, b`` up to total order 3, compared with ``exp(N(a+b))``. Result: all first-order terms and the
``a²`` term vanish (explicit Lorenz cycle is exact to order N), the ``ab`` and ``b²`` terms are equal
and linear in θ, ``κ c⋅θ``:

| N | Variant | c⋅θ | κ |
|---|---|---|---|
| 3 | A | θ₀ - θ₁ + 2θ₂ | 3/2 |
| 3 | B | θ₁ + θ₂ | 3/2 |
| 4 | A | θ₁ - θ₂ + θ₃ | 4 |
| 4 | B | θ₁ + θ₂ + θ₃ | 4/3 |

Comparing with a θ-method over the cycle (``R = (1+(1-α)Nb)/(1-αNb)``, error ``N²(α-1/2)b²``) gives
``α_eff = 1/2 + κ c⋅θ/N²`` (`effective_centering`). α_eff of earlier centerings: θ = 1/2 on every
substep 2/3 (N=3), 5/8 (N=4); A (0.069, 1.028, 0.48) and (0.116, 1.039, 0.462) 0.500;
A (-1.03, 0.70, 1.05) 0.562; B (1.386, 0.468, -0.464) 0.501; B (1.222, 1.142, -0.905) 0.539.
α_eff predicts the damping in the model: second order (α_eff = 1/2) keeps most gravity wave
energy, α_eff = 0.52 already damps considerably more over 2 days.

Design procedure: fix α_eff (one linear constraint, solved for the last θ with nonzero c), maximise
over the remaining N-1 centerings the explicit stability limit ω_E,max (largest value on a ladder
0.200, 0.225, ... up to the explicit limit of the cycle, 0.575 for N=3 and 0.705 for N=4, found by
bisection) for which all test cases of "Robust centerings" below are stable, tie-break by the
minimum gravity wave amplification for ω_LΔt ∈ [1, 4]. Bounded |θ| ≤ 3 or θ ∈ [0, 1] ("box").
Random search, 24 starts (some seeded with known centerings), 400 iterations, step 0.3 halved every
100. Results (amplification per step, minimum over ωΔt ≤ 1 and ωΔt ∈ [1, 4]):

| N | Variant | α_eff | box | θ | ω_E,max | amp ≤ 1 | amp 1-4 |
|---|---|---|---|---|---|---|---|
| 3 | A | 0.50 | no | (0.155, 2.131, 0.988) | 0.575 | 0.83 | 0.72 |
| 3 | A | 0.52 | yes | (0.912, 0.99, 0.099) | 0.575 | 0.977 | 0.86 |
| 3 | A | 0.55 | yes | (0.08, 0.999, 0.61) | 0.575 | 0.956 | 0.93 |
| 3 | B | 0.50 | no | (1.631, 0.55, -0.55) | 0.575 | 0.989 | 0.97 |
| 3 | B | ≤ 0.55 | yes | - | 0 | | |
| 4 | A | 0.50 | no | (2.056, 0.931, 1.02, 0.089) | 0.70 | 0.994 | 0.88 |
| 4 | A | 0.52 | no | (2.276, 0.899, 0.905, 0.087) | 0.70 | 0.996 | 0.90 |
| 4 | A | 0.50 | yes | (0.992, 0.65, 0.707, 0.057) | 0.325 | 0.96 | 0.74 |
| 4 | A | 0.60 | yes | (0.985, 0.979, 0.641, 0.062) | 0.625 | 0.84 | 0.77 |
| 4 | B | 0.60 | no | (2.977, -0.884, 1.713, 0.372) | 0.70 | 0.96 | 0.97 |

The N=4 A optimum is reproducible (independent starts: θ ≈ (2, 0.9, 0.9-1.0, 0.09)), robust to
θ₀ ∈ [2.0, 2.1] and rounding to (2.06, 0.93, 1.02, 0.09), but not to (2.06, 0.90, 1.00, 0.10)
(unstable in the linear tests).

Model sweep (method as in "Methods", η variance relative to initial, Δt_T32 = 20/30/40/50/60 min):

| Scheme | Random waves T128 | Blob T128 | Jet T85 |
|---|---|---|---|
| Leapfrog | 0.38/0.35/0.29/0.23/0.18 | 0.46/0.46/0.45/0.44/0.36 | stable to 50 |
| **4A (2.06, 0.93, 1.02, 0.09)** | 0.39/0.32/0.22/0.05/0.09 | 0.47/0.50/0.49/0.24/0.39 | stable to 40 |
| 4A (2.276, 0.899, 0.905, 0.087), α_eff 0.52 | 0.02/0.03/0.06/0.02/0.06 | 0.17/0.14/0.14/0.09/0.11 | stable to 40 |
| 4A box (0.985, 0.979, 0.641, 0.062) | 0.01 everywhere | 0.05/0.03/0.03/0.02/0.02 | stable to 40 |
| 4B (2.977, -0.884, 1.713, 0.372) | 0.01 everywhere | 0.04-0.01 | stable to 40 |
| **3B (1.631, 0.55, -0.55)** | 0.26/0.10/0.06/0.05/0.03 | 0.45/0.36/0.29/0.32/0.25 | stable to 40 |
| 3A box (0.912, 0.99, 0.099) | 0.02-0.01 | 0.20/0.15/0.13/0.11/0.10 | stable to 20 |
| 3A box (0.08, 0.999, 0.61) | 0.01 everywhere | 0.10-0.04 | stable to 40 |

Resolutions, 4A (2.06, 0.93, 1.02, 0.09) vs leapfrog (gravity wave cases 2 days, jet 20 days):

| Case | T32 | T64 | T128 | T256 |
|---|---|---|---|---|
| Random waves, 4A | 0.08/0.02/0.02/0.01/0.01 | 0.22/0.10/0.08/0.05/0.02 | 0.39/0.32/0.22/0.05/0.09 | 0.46/0.42/0.42/0.27/0.34 |
| Random waves, LF | 0.07/0.08/0.06/0.05/0.05 | 0.25/0.17/0.15/0.11/0.11 | 0.38/0.35/0.29/0.23/0.18 | 0.42/0.41/0.37/0.31/0.35 |
| Blob, 4A | 0.28/0.16/0.14/0.04/0.09 | 0.49/0.39/0.28/0.26/0.16 | 0.48/0.50/0.49/0.24/0.39 | 0.48/0.48/0.48/0.43/0.50 |
| Blob, LF | 0.34/0.27/0.26/0.24/0.25 | 0.45/0.36/0.37/0.37/0.29 | 0.46/0.46/0.45/0.44/0.36 | 0.48/0.47/0.46/0.46/0.46 |
| Jet, 4A stable up to | ≥ 60 | 40 | 30 | 30 |
| Jet, LF stable up to | ≥ 60 | ≥ 60 | 40 | 40 |

All 40 gravity wave runs with 4A were stable; the jet has a clean stability threshold at every
resolution (stable below, unstable above), ≈ 0.75x leapfrog's time step, matching the ratio of
the explicit limits (0.70 vs 1).

### Robust centerings (previous defaults)

The previous defaults, optimised for gravity wave energy, turned out unpredictably unstable in
the model: A with (0.116, 1.039, 0.462) blew up in the Gaussian blob example at Δt_T32 = 30 and
40 min but not at 50, in the random waves example at 30 and 40 but not at 50 and 60.

Re-optimised with stability first (method as in "Methods" with these changes):
- Extended test set, each maximised over ω_LΔt ∈ [0, 40]: pure gravity waves; orography
  β ∈ {0.1, 0.3, 0.5, 0.7, 0.9, 1}; Doppler ω_EΔt; Doppler combined with orography
  ``E = iω_E - βiω_L``, β ∈ {0.3, 0.7}; f-plane with ``fΔt = ω_EΔt``; with ω_EΔt (5 values) up to
  a level ω_E,max.
- Objective: the largest ω_E,max on the ladder 0.300, 0.325, ..., 0.575 for which all cases are
  stable (≤ 1 + 1e-6), subject to gravity wave amplification ≥ d_min for ω_LΔt ∈ [0.05, 1] and
  ≤ 1 for ω_LΔt ∈ [0.05, 40]. Tie-break among equal ω_E,max by high-frequency damping. 28
  starts, 600 iterations, step 0.4 halved every 150.
- d_min ∈ {0.99, 0.98, 0.95} all reach the full explicit limit ω_E,max = 0.575; optima are well
  defined (top starts agree). A tie-break favouring *little* damping over ω_LΔt ≤ 10 found
  optima with large centerings (e.g. A (7.2, -4.1, 10.1)) that damped even more in the model,
  so the linear gravity wave damping does not predict the energy loss in the model; only the
  stability transfers.

Model sweep (variance of η relative to initial, X = blew up; gravity waves T128 2 days with
random waves or the Gaussian blob, Galewsky jet T85 20 days), Δt_T32 = 20/30/40/50/60 min:

| Scheme | Random waves | Gaussian blob | Galewsky jet |
|---|---|---|---|
| Leapfrog | 0.38/0.35/0.29/0.23/0.18 | 0.46/0.46/0.45/0.44/0.36 | stable to 50 |
| A (0.116, 1.039, 0.462) (previous) | 0.38/X/X/0.16/0.18 | 0.47/8.19/7.46/0.40/0.34 | stable to 20 |
| B (1.386, 0.468, -0.464) (previous) | 0.32/0.25/0.19/0.17/0.15 | 0.46/0.43/0.36/0.41/0.32 | stable to 40 |
| **A (-1.03, 0.70, 1.05)** | 0.03/0.17/0.12/0.07/0.04 | 0.10/0.10/0.11/0.15/0.13 | stable to 40 |
| **B (1.222, 1.142, -0.905)** | 0.18/0.11/0.08/0.04/0.02 | 0.18/0.22/0.22/0.23/0.18 | stable to 40 |

The previous B centering was not unstable in these tests but only reaches ω_E,max = 0.425
in the linear analysis, so it was replaced as well. Rounding: (1.2, 1.14, -0.91) is linearly
unstable, use three decimals for B.

### Implicit solver as part of the time stepper

The implicit solver depends on the time stepper (e.g. ξ = 2αΔt for leapfrog, ξ = wθ_kΔt for the
Lorenz N-cycle), so it moves from `model.implicit` to `model.time_stepping.implicit`:

- `Leapfrog`, `NCycleLorenz` and `EulerForward` get a field `implicit`. Default for `Leapfrog` and
  `NCycleLorenz` is the placeholder `DefaultImplicit()`, for `EulerForward` (ocean/land) `nothing`.
- `default_implicit(spectral_grid, M)`: `nothing` for `Barotropic`, `ImplicitShallowWater` for
  `ShallowWater`, `ImplicitPrimitiveEquation` for `PrimitiveEquation`.
- The generic model constructor `M(spectral_grid; kwargs...)` replaces a `DefaultImplicit()` in
  `time_stepping` by `default_implicit(spectral_grid, M)` (`resolve_implicit`, rebuilding the
  time stepper with `with_implicit`), so `Leapfrog(spectral_grid)` passed to any model works as
  before. Passing `implicit` to a model constructor throws an `ArgumentError` pointing to the
  time stepper. The models' default `time_stepping` is constructed with
  `implicit = default_implicit(spectral_grid, <model type>)`.
- `variables(time_stepping, model)` collects the variables of the implicit solver (e.g. the
  average temperature profile of `ImplicitPrimitiveEquation`), previously collected as a model
  component.
- All uses of `model.implicit` now read `model.time_stepping.implicit`; function barriers that
  took `implicit` as an argument are unchanged.
- Breaking: models no longer have an `implicit` field (one type parameter less) and the
  constructor keyword `implicit` is removed.

## Summary of changes

- Fix η tendency allocation for the shallow water model (`ts` steps).
- Read fused spectral intermediates with `get_tendency_step` (shallow water and primitive models).
- `NCycleLorenz` option `implicit_centering`: per-substep centering, defaults
  `(2.06, 0.93, 1.02, 0.09)` for `steps=4` with `NCycleLorenzA`, `(-1.03, 0.70, 1.05)` and
  `(1.63, 0.55, -0.55)` for `steps=3` with `NCycleLorenzA`, `NCycleLorenzB`, `nothing`
  (= `implicit.centering` on every substep) otherwise. Used by the shallow water implicit solver.
- `NCycleLorenz` default `steps = 4` (was 3).
- `effective_centering(variant, θ)`: α_eff of the cycle from the order condition.
- Docs: Lorenz N-cycle algorithm as implemented, section on implicit centering, implicit solver
  section in `time_integration.md`, gravity wave example from a Gaussian blob.
- Breaking: implicit solver as part of the time stepper, `model.time_stepping.implicit`.

## Testing and verification

- Unit tests: default/explicit `implicit_centering` per substep; `effective_centering` (θ = 1/2 gives
  2/3 for N=3 and 5/8 for N=4, the second-order defaults 1/2, and the expansion of a cycle of the
  oscillation equation agrees to second order with a θ-method with centering α_eff); shallow water
  gravity waves with the 4-cycle A and 3-cycle A and B at Δt_T32 = 30 and 60 min run without NaNs
  and do not increase the η variance.
- Model experiments as described in "Methods" below.

## Documentation changes

`docs/src/time_stepping.md`, Lorenz N-cycle section incl. "Implicit centering of the Lorenz N-cycle"
(order condition, α_eff, signs of the centerings, design procedure, defaults);
`docs/src/shallowwater.md` semi-implicit Lorenz N-cycle; `docs/src/time_integration.md`.

## Known limitations

- The 4-cycle A default needs ≈ 0.75x leapfrog's time step where advection limits (Galewsky jet),
  the explicit limit of the cycle. It preserves gravity wave energy like leapfrog at T128/T256 for
  Δt_T32 ≲ 40 min but damps more at T32/T64 and longer time steps; the η variance is not monotonic
  in Δt (dips at 50 min) although stable.
- The 4-cycle A centering is sensitive in some directions (stable for θ₀ ∈ [2.0, 2.1], not for
  (2.06, 0.90, 1.00, 0.10)), so it is pinned to two decimals.
- No robust second-order centering for the 3-cycle variant A (default slightly off-centred,
  α_eff = 0.56).
- AB and ABBA: the centering tuple has length N and is indexed by `step_counter mod N`, so A and B
  cycles cannot get different centerings; they use `implicit.centering` on every substep.
- The implicit primitive equation solver does not use the per-substep centering.
- Implicit hyperdiffusion uses Δt instead of wΔt.

## Future work

- Consider `NCycleLorenzB` as default variant; understand why it is more stable than A in the model.
- Centering indexed over the full AB/ABBA sequence.
- Third-order conditions for the centerings (the a²b, ab², b³ terms are quadratic in θ), weighted implicit diffusion, primitive equations.
