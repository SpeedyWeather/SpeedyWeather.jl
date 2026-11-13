# NCycleLorenz for the shallow water model: bug fixes and implicit centering

> Status: **in progress**. Bugs fixed and per-substep implicit centering implemented for the 3-cycle variants A and B; N=4, AB, ABBA left for review.

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

Daisuke's version as tested here (`impl.jl`, ``μ = 0``): ``G ← wF_E + (1-w)G`` with
``F_E = F - Lu``, ``u ← u + Δt(1-αΔtL)^{-1}(G + Lu)``. The Lawson variant (``μ = 1``) additionally
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
over a full cycle (scripts in `ncycle-lorenz-shallow-water/`). Test cases: pure gravity wave
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
| **Milan's version, θ = (0.116, 1.039, 0.462), A (this PR)** | ≥ 0.997 for ωΔt ≤ 1.5 | stable to 0.4 | stable | stable |
| **Milan's version, θ = (1.386, 0.468, -0.464), B (this PR)** | ≥ 0.998 for ωΔt ≤ 1.5 | stable to 0.4 | stable | stable |
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

| Δt at T32 | Leapfrog | main (θ=0.5) | A, θ = (0.069, 1.028, 0.48) | **A, θ = (0.116, 1.039, 0.462)** | **B, θ = (1.386, 0.468, -0.464)** |
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

## Summary of changes

- Fix η tendency allocation for the shallow water model (`ts` steps).
- Read fused spectral intermediates with `get_tendency_step` (shallow water and primitive models).
- `NCycleLorenz` option `implicit_centering`: per-substep centering, default
  `(0.116, 1.039, 0.462)` for `steps=3` with `NCycleLorenzA`, `(1.386, 0.468, -0.464)` with
  `NCycleLorenzB`, `nothing` (= `implicit.centering`
  on every substep) otherwise. Used by the shallow water implicit solver.
- Docs: Lorenz N-cycle algorithm as implemented, section on implicit centering.

## Testing and verification

- Unit tests: default/explicit `implicit_centering` per substep; shallow water gravity waves with
  `NCycleLorenz` (variants A and B) run without NaNs and keep at least half the η variance of
  leapfrog after 1 day.
- Experiment scripts in `ncycle-lorenz-shallow-water/` (`one.jl`, `one_theta.jl`, `jet1.jl`).

## Documentation changes

`docs/src/time_stepping.md`, Lorenz N-cycle section.

## Known limitations

- Variant A still needs a shorter time step than leapfrog (≈0.5x for the gravity wave example and
  the Galewsky jet at T85). Variant B gets to ≈0.8x (jet) and ≈1x (gravity waves) but its
  centering is a less well-defined optimum. The default variant is unchanged (A).
- AB and ABBA: the centering tuple has length N and is indexed by `step_counter mod N`, so A and B
  cycles cannot get different centerings; they use `implicit.centering` on every substep.
- The implicit primitive equation solver does not use the per-substep centering.
- Implicit hyperdiffusion uses Δt instead of wΔt.

## Future work

- Consider `NCycleLorenzB` as default variant; understand why it is more stable than A in the model.
- Centering indexed over the full AB/ABBA sequence.
- N=4 centering (higher explicit limit ≈0.71), weighted implicit diffusion, primitive equations.
