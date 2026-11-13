# NCycleLorenz for the shallow water model: bug fixes and implicit centering

> Status: **in progress**. Bugs fixed and per-substep implicit centering implemented for the 3-cycle; N=4 left for review.

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

- Plan written after the investigation (exploratory work), bundling bug fixes, the analysis and
  the per-substep centering for N=3. N=4 centerings are only reported, not made default.

## Problem description

`NCycleLorenz` (#1031, #1035) is not the default for the shallow water model as it was found to
be more diffusive and to need a shorter time step than `Leapfrog`.

## Background

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
| N-cycle as in code, θ=0.5 (N=3) | 0.90, 0.75 | unstable (1.07) | stable | stable |
| L outside Lorenz average (docs), α=0.5 | 1, 1 | unstable (1.2-1.5) | unstable | - |
| Lawson / integrating factor | 1, 1 | stable | unstable (1.04-1.5) | - |
| Code, per-substep θ = (1, 1, 0), N=3 | exactly 1 (all-pass) | stable to 0.55 | stable | unstable (1.1-2.6) |
| **Code, θ = (0.069, 1.028, 0.48), N=3** | 1, 0.9997 | stable to ≈0.3 | stable | stable |

The explicit Lorenz N-cycle itself is limited to ω_EΔt < √3/3 ≈ 0.58 (N=3) and ≈ 0.71 (N=4)
per tendency evaluation, compared to 1 for leapfrog. The robust optimum for N=4 found so far,
θ ≈ (2.31, 0.86, 0.86, 0.08), damps gravity waves more (≈0.994/step) and is not tested in the model.

The hyperdiffusion is applied implicitly with Δt before the Lorenz weighting; for N=3 the
weight w=3 makes stiff diffusion unstable (|A| ≈ 1.2 for νΔt ≳ 10). Using wΔt makes it stable,
but in the experiments this was not the limiting factor (not changed here).

### Experiments: docs gravity wave example (T128, 2 days, var(η)/var(η₀))

| Δt at T32 | Leapfrog | NCycle N=3 (old, θ=0.5) | NCycle N=3, θ = (0.069, 1.028, 0.48) |
|---|---|---|---|
| 10 min | 0.36 | 0.011 | 0.43 |
| 20 min | 0.38 | 0.011 | 0.31 |
| 30 min | 0.35 | 0.011 | 0.30 |
| 40 min | 0.29 | 0.011 | unstable |
| 60 min | 0.18 | 0.011 | - |

Galewsky jet (`ZonalJet`, 20 days): T85 stable up to Δt_T32 = 20 min (old 25 min, leapfrog ≥ 40 min),
T31 up to 30 min (leapfrog ≥ 40 min). The L-outside-average and Lawson variants were unstable
beyond Δt_T32 ≈ 10-20 min in the gravity wave example.

## Summary of changes

- Fix η tendency allocation for the shallow water model (`ts` steps).
- Read fused spectral intermediates with `get_tendency_step` (shallow water and primitive models).
- `NCycleLorenz` option `implicit_centering`: per-substep centering, default
  `(0.069, 1.028, 0.48)` for `steps=3` with `NCycleLorenzA`, `nothing` (= `implicit.centering`
  on every substep) otherwise. Used by the shallow water implicit solver.
- Docs: Lorenz N-cycle algorithm as implemented, section on implicit centering.

## Testing and verification

- Unit tests: default/explicit `implicit_centering` per substep; shallow water gravity waves with
  `NCycleLorenz` run without NaNs and keep at least half the η variance of leapfrog after 1 day.
- Experiment scripts in `ncycle-lorenz-shallow-water/` (`one.jl`, `one_theta.jl`, `jet1.jl`).

## Documentation changes

`docs/src/time_stepping.md`, Lorenz N-cycle section.

## Known limitations

- NCycleLorenz still needs a shorter time step than leapfrog (≈0.75x in the gravity wave
  example, ≈0.5x for the Galewsky jet at T85), from the explicit stability limit of the 3-cycle.
- The implicit primitive equation solver does not use the per-substep centering.
- Implicit hyperdiffusion uses Δt instead of wΔt.

## Future work

- N=4 centering (higher explicit limit ≈0.71), weighted implicit diffusion, primitive equations.
