# StochasticLeapfrog: additive Wiener noise, upstreamed from FokkerPlanckSensitivities.jl

> Status: **in progress**. Implemented and unit-tested on `mg/gfdt-upstream` (draft PR #1397). Still open: rebase onto `main`, full test suite, GPU, docs build, cleanup before review.

Date of initial draft: 2026-11-13 (this status document; work on the branch started 2026-11-11)

Base revision: `24fdb6847a3bdd49f7e23e94a5b64de06dd89d48` (branch point). HEAD when written: `6892775024a22f609652757eecfb960d8ad380ac`.

## Originating prompt

> The prompts that started this work came from earlier sessions and were not recorded. This
> document was written in response to: "First inspect this branch to understand what we are doing
> here in this branch, then make a markdown to save this understand for yourself for later, so we
> can continue"

## Revision log

> - 2026-11-11: First `StochasticLeapfrog` commit. Random numbers drawn on the model's
>   architecture (`architecture_rng`). Random-number buffer in the AR1 process.
> - 2026-11-12: `SpectralAR1Process` gets `n_patterns`. Renamed `random_process!(vars, process)` →
>   `update_random!`. Amplitudes keyed by tendency name, with the variables inferred from the
>   model's tendencies instead of a hard-coded list. Up to `52c666f4` the code reproduced the origin
>   project **bitwise**, with the noise added before the dynamical core.
> - 2026-11-13, **user decision**: the noise moves to right before the implicit terms. The call
>   stays **explicit in `time_step!`** (the user did not want it hidden inside a
>   `diffusion_and_implicit!` method). All atmospheric noise is now spectral. Results match the
>   origin project only to rounding. Docs and references added.
> - 2026-11-13, later: `wiener_tendencies!` made type-stable (`772a26d1`), with a JET guard in
>   `test/dynamics/dispatch.jl`. Noise-factor formulas and the `StochasticLeapfrog` docstring
>   switched to typeset LaTeX. The `# IMPLEMENTATION NOTES` block was removed from the source
>   (`1066382e`). The PrimitiveDry cases and the "noise before the implicit terms" testset were
>   removed from `test/dynamics/stochastic_leapfrog.jl` to shorten the tests (`c8946a7b`).

## Problem description

The sister project `/Users/max/Nextcloud/FokkerPlanckSensitivities.jl` computes climate
sensitivities with the adjoint Fokker–Planck approach (Thuburn 2005). It needs SpeedyWeather to
integrate the SDE `dx = F(x) dt + a dW`, where the noise sets the diffusion `K = a²σ²/2` of the
Fokker–Planck equation. That project had its own noise code (a `SpectralWhiteNoiseProcess` and a
`GridWhiteNoiseProcess` under a `CombinedNoiseProcess` that drives an additive forcing). This
branch upstreams that capability into SpeedyWeather as a time stepper. "gfdt" in the branch name
refers to this response-theory / fluctuation-dissipation use.

## Background

- Before this branch: `AbstractRandomProcess` had only `SpectralAR1Process`, used as
  `model.random_process` by SPPT (`parameterizations/stochastic_physics.jl`). It was advanced with
  `random_process!(vars, process)` at the end of `update_prognostic!`.
- In the origin project, `speedy/checks/upstream_wiener_equivalence.jl` is the equivalence gate.
  It runs in the `speedy-upstream` env, which points at this checkout via `[sources]`
  (`julia +1.11 --project=speedy-upstream`). `speedy/checks/leapfrog_noise_placement.jl` is the
  scalar OU placement test. The origin project's log is `docs/dev/PROGRESS.md` (entries for
  2026-11-12 and 2026-11-13). Its pinned envs `speedy/` and `speedy-p4*` still reproduce old runs;
  the new upstream code does not reproduce them bitwise.

## Summary of changes

All changes are in `SpeedyWeather/` plus the docs. The other submodules are untouched.

### `src/dynamics/random_process.jl`: random process interface revised

- `update_random!(pattern, process)` draws the next pattern into any array; AR1 reads its previous
  state from that array. This lets any component own a process and its own pattern.
  `update_random!(vars, process)` advances `model.random_process`, by default in
  `vars.prognostic.random_pattern`.
- `architecture_rng(arch, seed; draw_on_cpu)` returns `Random.Xoshiro` on CPU and
  `GPUArrays.RNG` on GPU. `draw_on_cpu = true` always uses Xoshiro (the same numbers on every
  architecture). The GPUArrays compat was raised to `11.5.1` for this.
- `complex_normals!(r, rng)` draws in place and scales by `2/√2`, written that way to keep the
  random numbers of earlier versions bit for bit.
- `SpectralAR1Process` gains `n_patterns` and a `random_numbers` buffer, so a draw no longer
  allocates. `random_pattern_variables(n)` gives 2D variables for `n ≤ 1` and `Spectral3D/Grid3D(; n)`
  otherwise.
- New `SpectralWhiteNoiseProcess`: memoryless and spectrally flat by default (`wavenumber = 0`).
  No noise on l = 0 or l = lmax. The m = 0 coefficients are real and scaled by √2, so the pattern
  is homogeneous and band-limited. No clamp by default.
- New `GridWhiteNoiseProcess`: independent normals at every grid point, drawn directly in grid space.

### `src/time_stepping/steppers/leapfrog.jl`

- `Leapfrog` gains a 7th type parameter `N` and a field `noise` (default `nothing`).
- `initialize!(L)` calls `initialize!(L.noise, model)`, and `variables(L, model) = variables(L.noise, model)`.

### `src/time_stepping/steppers/stochastic_leapfrog.jl` (new)

Until `1066382e` the file started with a long `# IMPLEMENTATION NOTES` block, the most detailed
design record. It is summarised in this document; the full text is in the history:
`git show 772a26d1:SpeedyWeather/src/time_stepping/steppers/stochastic_leapfrog.jl`.

- `StochasticLeapfrog` is a **type alias**: a `Leapfrog` whose `noise` is an `AbstractWienerNoise`.
  It is not a new stepper type. The constructor `StochasticLeapfrog(spectral_grid; seed, amplitudes...)`
  passes `Leapfrog` kwargs on to `Leapfrog` and everything else to `AdditiveWienerNoise`. Passing
  both `noise` and noise kwargs raises an `ArgumentError`.
- `AdditiveWienerNoise`: `amplitudes` is a NamedTuple keyed by tendency name, merged over
  `DEFAULT_WIENER_AMPLITUDES` (vor/div 3e-8, temperature 3e-3, humidity 3e-8, SST 1e-3). The
  defaults are tuned for T32 L8 PrimitiveWet from a 1500-year control run in the origin project.
  Other fields: `tapering(σ)` → `taper` per layer, `spectral::SpectralWhiteNoiseProcess`,
  `grid::GridWhiteNoiseProcess`. When both processes get the same nonzero seed, the grid one is
  reseeded with a derived hash.
- `tendency_variables(model)` collects the `TendencyVariable`s of the model and all its components
  except `time_stepping` (excluded to avoid recursion). It exists because `Variables` do not exist
  yet when `variables(noise, model)` is called.
- Variables: `vars.prognostic.wiener_pattern` (`Spectral3D`) has one slot per layer of every
  atmospheric tendency. Variables with layers come first, then single-layer ones, each in
  `tendency_names` order. Zero-amplitude variables keep their slots, so the realisation of one
  variable does not depend on the other amplitudes. `vars.prognostic.wiener_pattern_surface`
  (`Grid3D`) has slots **only for ocean/land variables with a nonzero amplitude**.
- `random_process!(vars, L::StochasticLeapfrog, model)` calls `wiener_tendencies!` (adds the stored
  pattern) and then `wiener_draw!` (draws the pattern for the next step).
- Scaling: the tendency receives `a/√Δt · taper · ξ`, with `Δt = L.Δt` in seconds. Atmospheric
  tendencies are radius-scaled at that point, so the amplitude is multiplied by `R`, and by `R²`
  for vorticity and divergence. Ocean/land amplitudes are not scaled.
- `wiener_tendencies!` dispatches to two `@generated` drivers, `wiener_atmosphere!` and
  `wiener_surface!`. They are unrolled over the names from the `Variables` type, like
  `update_prognostic_atmosphere!`. The per-variable work is in `wiener_slots!` and
  `wiener_surface_slots!`.

### `src/time_stepping/time_integration.jl`

- Both `time_step!` methods call `random_process!(vars, time_stepping, model)` between
  `dynamics_tendencies!` and `diffusion_and_implicit!`. The fallback for any
  `AbstractTimeStepper` is a no-op.
- With `model.dynamics = false` the hook is never called, so there is no noise.
- `update_prognostic!` now calls `update_random!(vars, model.random_process)`.
- **Naming trap:** `random_process!` now means the 3-argument stepper hook. The old 2-argument
  function is `update_random!`.

### Docs, CHANGELOG, tests

- `docs/src/time_stepping.md`: new section `[Stochastic leapfrog](@id stochastic_leapfrog)` covering
  the scheme, the references and an `@docs` block. `docs/src/time_integration.md`: usage example.
  `docs/src/references.bib`: `EwaldTemam2005`, `EwaldWeakSchemes`, `Higham2000`,
  `HighamMaoStuart2002`, `KloedenPlaten1992`.
- `CHANGELOG.md` entry for #1397. It already mentions the `random_process!` → `update_random!`
  rename, so the matching "Open" item in the implementation notes is done.
- Tests: `test/dynamics/stochastic_leapfrog.jl` (new, 7 testsets),
  `test/parameterizations/white_noise_process.jl` (new),
  `test/parameterizations/random_process.jl` (rename plus an `n_patterns` testset),
  `test/parameterizations/stochastic_physics.jl` (rename).

## Key design decisions and why

1. **Placement right before the implicit terms** (drift-implicit, stochastic-theta-like). Only the
   implicit terms care where the noise enters:
   - Implicit diffusion damps the noise within the step. In the scalar OU test (stationary variance
     divided by the exact value, noise added before vs after the implicit terms):
     2Δtν = 0.1 gives 1.28 vs 1.54, 0.34 gives 1.39 vs 2.49, and 3 gives 0.75 vs 12.0.
   - The semi-implicit correction damps the noise in fast gravity waves.

   The dynamical core is indifferent because it only accumulates into the tendencies. Adding the
   noise after the implicit terms (split-step, e.g. in `update_prognostic!`) was the user's earlier
   proposal and was rejected: it leaves stiff modes undamped.
2. **One increment `2a√Δt ξⁿ` per leapfrog step.** It gives the nominal diffusion `a²σ²/2` for
   leapfrog+RAW (1.00 in the scalar test; `a/√(2Δt)` would give 0.50). The computational mode is
   excited and the RAW filter damps it. Ewald's `ΔWⁿ⁻¹ + ΔWⁿ` would avoid exciting it; that is a
   possible option.
3. **The pattern is drawn at the end of a step and applied in the next**, and stored in the
   prognostic variables. As a result the first step of a new simulation is deterministic, and
   `run!` continues with the stored pattern. This was needed only for the bitwise reproduction.
4. **Mass is conserved** because the pattern has no l = 0 component, so no explicit mass fixer is
   needed. A test covers this.

## Testing and verification

As of 2026-11-13, according to the origin project's PROGRESS.md (**not re-run in the session that
wrote this document**):

- Upstream tests passed 718/718: stochastic leapfrog, random processes, white noise, SPPT, time
  stepping.
- The equivalence gate passed. After 5 steps the state matched to 1e-3 in Float32 (max 1.35e-4)
  and 6.3e-11 in Float64, so the difference is rounding. Noise patterns were bitwise equal; zero
  amplitudes gave a bitwise-equal `Leapfrog`; barotropic T8 Float64 matched to 2.9e-15 after 100 days.
- Not run: the full test suite, GPU, and the docs build.

Type-stability fix (2026-11-13, measured on this Mac, T32 L8 PrimitiveWet, every noise path active):

| | before | after |
|---|---|---|
| trajectory after 20 steps | reference | bitwise identical |
| allocations per `wiener_tendencies!` call | 22.4 KB | 3.5 KB (KernelAbstractions launches, ~560 B each) |
| time per call (noisy) | ~53 μs | ~27 μs |
| JET runtime-dispatch reports in the `StochasticLeapfrog` `time_step!` | 46 | 0 |

After the fix the stochastic tests passed 492/492 (stochastic leapfrog, white noise, random
process, SPPT; before the test shortening in `c8946a7b`). `dispatch.jl` passed 4/4, including
the new `StochasticLeapfrog` testset.

Commands, run from the repository root. `SpeedyWeather/test/Manifest.toml` is stale (July, no
SpeedyWeather entry), so the JET test below stacks the test env behind the package env instead of
resolving it:

```bash
julia +1.12 --project=SpeedyWeather --check-bounds=yes -e 'using Test, SpeedyWeather; include("SpeedyWeather/test/dynamics/stochastic_leapfrog.jl")'
julia +1.12 --project=SpeedyWeather --check-bounds=yes -e 'using Test, SpeedyWeather; include("SpeedyWeather/test/parameterizations/white_noise_process.jl"); include("SpeedyWeather/test/parameterizations/random_process.jl"); include("SpeedyWeather/test/parameterizations/stochastic_physics.jl")'
JULIA_LOAD_PATH="@:$PWD/SpeedyWeather/test:@stdlib" julia +1.12 --project=SpeedyWeather --check-bounds=yes -e 'using Test, SpeedyWeather; include("SpeedyWeather/test/dynamics/dispatch.jl")'
```

## Open items (before un-drafting #1397)

- [ ] **Rebase or merge `origin/main`**, which is 6 commits ahead (two-band shortwave #1262, NaN
      gradient fix #1401, precipitable water #1370, …). Expect conflicts in `CHANGELOG.md` and
      `SpeedyWeather/Project.toml`. `main` is already at `version = "0.25.0-DEV"`, which also covers
      the version bump required for a new public type; keep it. The branch currently reads `0.24.0`.
- [ ] Push: local commits are ahead of `origin/mg/gfdt-upstream`.
- [ ] Run the full `Pkg.test("SpeedyWeather")`, GPU tests, and the docs build. GPU could be tried
      locally with Metal (`SpeedyWeather/test/GPU/MetalGPU`).
- [x] Type instability in `wiener_tendencies!`: fixed with `@generated` drivers (`772a26d1`), see
      Testing and verification.
- [ ] Ocean/land slots exist only for nonzero amplitudes, unlike the atmosphere rule. This is a
      leftover of the bitwise reproduction, which was dropped on 2026-11-13, so it can now be made
      consistent.
- [ ] Draw-at-end-of-step plus storing the pattern in prognostics was also only needed for bitwise
      reproduction. Drawing and applying in the same step would make the pattern scratch memory.
      Decide whether to change this (it changes first-step determinism and restart semantics).
- [ ] `wiener_surface_slots!` passes `L` (the Leapfrog) to `get_tendency_step` for ocean/land
      tendencies instead of `time_stepper(L, Val(:ocean/:land))`. Harmless today (every stepper
      falls back to tendency step 1) but inconsistent with `ocean.jl` and `soil_temperature.jl`.
- [x] `# IMPLEMENTATION NOTES` block removed from `stochastic_leapfrog.jl` (`1066382e`).
- [ ] Write the PR description (the body is currently empty).
- [ ] Maybe: an Ewald-style `ΔWⁿ⁻¹ + ΔWⁿ` increment option.

## Known limitations

- The default amplitudes are for T32 L8 PrimitiveWet only.
- The effective diffusion matches the nominal `a²σ²/2` only where neither the implicit terms nor
  the filter matter (1.33× nominal for a T8 barotropic model). It should be measured for the
  scales of interest.
- `model.dynamics = false` means no noise.
- `GridWhiteNoiseProcess` is uniform in grid coordinates, not homogeneous on reduced grids.

## Future work

- Multiplicative or state-dependent noise could be another `AbstractWienerNoise` subtype.
- A stochastic NCycleLorenz would need its own `random_process!` method.
