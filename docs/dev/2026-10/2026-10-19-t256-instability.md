# T256L8 PrimitiveWetModel instability in long runs

> Status: **in progress**. The instability is confirmed on GPU (Float32 and Float64) and on CPU. The decisive
> factor is the **shape** of the semi-implicit reference temperature profile, which is frozen at the initial
> conditions: a reference with at least the model's top-level static stability (isothermal 300 K, or the
> current profile refreshed periodically) removes the cold spots at T256 and T128. A linear stability analysis
> shows the stale reference is *not* linearly unstable, so the mechanism involves the flow over the Antarctic
> orography. In a 2-year, 2-seed ensemble both cures survived all runs, while both default runs blew up.

Date of initial draft: 2026-10-19

Base revision: `e410eed0` (branch `mg/instability`, identical to `main`)

## Originating prompt

> In another investigation I saw an issue that the T256L8 model on GPU is instable. I need to investigate
> this. First confirm it, both on GPU and CPU. Then we need to investigate the cause

## Revision log

> while we wait, write down the current progress of the investigation in a markdown file

- 2026-10-19: initial draft, written while the variant experiments were running.
- 2026-10-19 (update): CPU confirmation, results of all variant experiments, T128 with the cures, linear
  stability analysis of the semi-implicit scheme (finding 5), 2-year ensemble started. The hypothesis changed
  from "reference too cold" (Simmons et al., 1978) to "reference with the wrong stratification at the top":
  a reference shifted by +30 K but with the same shape still blows up, and the linear analysis rules out the
  classical semi-implicit instability.
- 2026-10-19 (update): 2-year ensemble result (finding 6): both cures stable in both seeds, default 8/8 unstable.

## Problem description

The default `PrimitiveWetModel` at T256L8 (Float32, `OctahedralGaussianGrid`, Δt = 300 s, all default
components) produces NaNs after 7 to 12 months of simulated time. At T128L8 the same configuration runs
stably for at least 10 years.

## Background

The instability was first seen in the Reactant DFT work (`docs/dev/2026-10/2026-10-16-reactant-spectral-transform.md`
on `mg/reactant-dft`, scripts in `/p/tmp/maxgelbr/code/Speedy2/claude_reactant_gpu_2026-10-13/dft_impl_2026-10-18/`).
Long-run stability controls there (H100, 2 years each) found:

| run | Slurm | last finite | first NaN |
|---|---|---|---|
| CUDA.jl Float32, `mg/reactant-dft` | 2568562 | day 180 | day 210 |
| CUDA.jl Float32, `main` at `e410eed0` | 2582219 | day 270 | day 300 |
| Reactant DFT, `HIGHEST` precision | 2568563 | | after day 330 |
| Reactant DFT, TF32 | 2568499 | | after day 360 |
| CUDA.jl Float32 **T128L8**, 10 years | 2569913 | stable | |

The global means stayed steady in all runs. The minimum temperature fell steadily over the southern
winter, from ~190 K to ~136 K, before the blow-up. Reactant, TF32 and the DFT were therefore not the cause.

## Setup of this investigation

Scripts, logs and output are in `/p/tmp/maxgelbr/code/claude_t256_instability_2026-10-19/`, outside the repository:

- `stability_run.jl`: a long-run driver configured through environment variables (`ARCH`, `TRUNC`, `NF`, `DAYS`,
  `VARIANT`, `DT_T32_MIN`, `RESTART`, `IMPLICIT_PROFILE`, `IMPLICIT_UPDATE_DAYS`, `IMPLICIT_FIXED`). It steps with
  `time_step!(simulation)` after a single `initialize!(simulation; period)`, like `run!`. Every `DIAG_HOURS` it writes
  diagnostics with the *location* (lat, lon, layer) of the extremes: T_min, T_max, max |u|, |v|, |ζ|, |D|,
  min p_s, and the number of grid points below 170 K. Output goes to stdout and to a CSV. It snapshots the
  prognostic variables (JLD2) every 30 days, and on a blow-up it dumps the last 3 daily states as restart points.
- `analyze_snapshot.jl`: loads snapshots into a CPU model and reports column static stability (potential
  temperature differences between adjacent layers), the column at the coldest point, and polar-cap statistics.
- `linear_stability.jl`: linear stability analysis of the semi-implicit leapfrog with a reference profile that
  differs from the actual one (finding 5).
- `summarize.py`, `status.sh`: compact views of the CSVs and logs.
- Further driver options: `IMPLICIT_FIXED` (`iso300` or `plusN`), `PERTURB_SEED` (0.01 K white noise on the
  initial temperature, for independent realizations), and per-layer leapfrog noise in the CSV (`noiseD_k*`,
  `noiseT_k*`: rms over spectral coefficients of the difference between the two leapfrog time levels).
- Node `cso14c65` does not mount `/home`, so the Julia depot is missing there; the CPU jobs exclude it.
- `job_gpu1.sh` (one run per GPU, `gpupreempt`) and `job_cpu*.sh` (`medium`) are the Slurm scripts. The Julia
  environment is `env/` (CUDA.jl 6.4.2, `CUDA_Runtime_jll` pinned to 13.0, SpeedyWeather from this checkout).

Restarts need care. A run started from a snapshot rebuilds the semi-implicit operators from the restart state
(see below), so it is not the same model as the uninterrupted run. `IMPLICIT_PROFILE=initial` keeps the
reference profile of the initial conditions, as in an uninterrupted run.

## Findings so far

### 1. Confirmed on GPU, independent of floating-point precision

| run (H100, from day 0) | Slurm | result |
|---|---|---|
| T256L8 Float32, default | 2582904 | **NaN at day 280** (5 Oct) |
| T256L8 Float64, default | 2582907 | **NaN at day 358** (22 Dec) |
| T128L8 Float32, default | 2582906 | stable for 400 days |

A default run with 0.01 K white noise added to the initial temperature (`ens_default_s2`, seed 2) blew up at
day 181. Including the earlier runs, the T256 default configuration has blown up in all seven attempts, between
day 181 and day 358. Float64 delays but does not prevent it, so Float32 rounding is not the cause.

A restart from the day-270 snapshot of the Float32 run, with the initial reference profile, blew up at day
279.75 (the original at day 280). The snapshot restart reproduces the original run.

### 2. Confirmed on CPU

- On day 1 the CPU and GPU states agree to 3–4 digits (T_min 182.4 K at 35.3°N on both, max |ζ| 7.04e-4 vs
  7.05e-4 s⁻¹), so both architectures integrate the same model.
- **The CPU restart from the GPU's day-270 snapshot (initial reference profile, `R270_init_cpu`, Slurm 2585789)
  blew up at day 279.75**, exactly like the GPU restart from the same snapshot (day 279.75) and the original
  GPU run (day 280). The instability is not GPU-specific.
- The CPU restart from day 240 and the fresh CPU run (Slurm 2582905) are still running. CPU throughput is
  ~197 s per simulated day (~1.2 SYPD with 32 threads, barely faster than single-threaded), so the fresh CPU
  run reaches the blow-up window only after ~15 h.

### 3. Signature: cold spots and unstable stratification over Antarctica

The minimum temperature is in the **top layer (k1, σ = 0.0625) at 75–88°S**, mostly over the slopes of the East
Antarctic plateau (surface height 2–3.5 km). Between days 150 and 400:

| run | lowest T_min | max. number of points < 170 K | max \|u\| |
|---|---|---|---|
| T256 Float32 default | 125.8 K (day 276, 80°S, k1) | 1949 | 109 m/s |
| T256 Float64 default | 141.7 K (day 291, 83°S, k1) | 1620 | 126 m/s |
| T128 Float32 default (stable) | 142.4 K (day 177, 75°S, k1) | 375 | 104 m/s |

The cold spots persist into November–December, when Antarctica is sunlit, so they are not simply polar-night
radiative cooling. The `OneBandLongwave` scheme is self-limiting (emission ∝ T⁴).

The snapshot analysis (`analyze_snapshot.jl`) shows:

- The top layers are **not** convectively unstable. θ_k1 − θ_k2 stays above +30 K everywhere, even over the
  coldest spots, so a superadiabatic top layer is ruled out.
- The **mid-troposphere around Antarctica becomes strongly superadiabatic** at T256. On day 277,
  θ_k4 − θ_k5 = −41 K at 62°S 83°E (T = 195 K above 244 K, about 29 K/km), with 1000–4400 unstable columns
  per layer pair from k3 to k7. At T128 in the same season, the minimum is −18 to −25 K with 130–730 unstable
  columns. The model is hydrostatic and has no dry convective adjustment: Betts–Miller needs moisture, and
  `BulkRichardsonDiffusion` acts only in the boundary layer.
- The maxima of |v| and |ζ| sit over and around Antarctica nearly all year, at T128 as well. Examples are
  |v| ≈ 80 m/s at 82°S, σ = 0.44 in late December, and |ζ| ≈ 5–7e-4 s⁻¹, about 4× the polar Coriolis
  parameter. Some of this is spurious small-scale flow over the steep orography.

### 4. Lead hypothesis: the semi-implicit reference temperature profile is stale

`ImplicitPrimitiveEquation` builds its operators from `vars.dynamics.average_temperature_profile` whenever Δt
changes, which only happens in the 3-step leapfrog spin-up (`reinitialize!`). In a run started from the
default initial conditions, the reference profile is therefore the **initial-condition profile for the whole
run**. The model climate drifts far from it:

| layer (σ) | k1 (0.06) | k2 (0.19) | k3 (0.31) | k4 (0.44) | k5 (0.56) | k6 (0.69) | k7 (0.81) | k8 (0.94) |
|---|---|---|---|---|---|---|---|---|
| reference, frozen at t = 0 [K] | 200.0 | 225.6 | 243.1 | 255.2 | 264.7 | 272.5 | 279.2 | 285.1 |
| global mean, day ~300 [K] | 231.8 | 226.1 | 228.4 | 239.9 | 251.0 | 260.2 | 268.5 | 276.9 |
| actual − reference [K] | **+32** | +1 | −15 | −15 | −14 | −12 | −11 | −8 |

At the top, the reference is 32 K too cold and has the wrong sign of the lapse rate: the reference cools with
height, while the model has a stratospheric inversion. Semi-implicit schemes are known to become unstable
when the actual temperatures exceed the reference temperature (Simmons, Hoskins and Burridge, 1978). This is
why operational models use a warm isothermal reference (~300 K).

Restarts from the same snapshots, differing only in the reference profile:

| experiment | reference profile | result |
|---|---|---|
| restart from day 240 (`R240_init`) | initial conditions (as the original) | NaN at day 293 |
| restart from day 240 (`R240_curr`) | rebuilt from the day-240 state | stable to day 300 (end of run), T_min ≥ ~170 K |
| restart from day 270 (`N270_initial`) | initial conditions | NaN at day 279.6 |
| restart from day 270 (`N270_current`) | rebuilt from the day-270 state | stable to day 280 (end of run) |

Fresh T256L8 runs from the default initial conditions, one realization each (Slurm 2585790–2585796 and
earlier), statistics over days 150–365 ("days n>0": number of daily diagnostics with any point below 170 K):

| run | change | outcome | lowest T_min | max. points < 170 K | days n>0 |
|---|---|---|---|---|---|
| default F32 | — | NaN day 280 | 125.8 K (k1) | 1949 | 63 |
| default F64 | `NF = Float64` | NaN day 358 | 141.7 K (k1) | 1620 | 119 |
| `orog_T128` | `EarthOrography(; smoothing_fraction = 0.5)` | **NaN day 202** | 155.8 K (k1) | 228 | 8 |
| `plus30` | reference = initial profile + 30 K | **NaN day 286** | 148.3 K (k1) | 1138 | 90 |
| `dt_half` | Δt = 150 s | stable 365 d | 128.3 K (k1) | 735 | 109 |
| `diff_2h` | hyperdiffusion `time_scale = Hour(2)` | stable 365 d | 152.8 K (k3) | 684 | 71 |
| `update30` | reference refreshed every 30 days | stable 365 d | 164.2 K (k1) | 78 | 33 |
| `update1` | reference refreshed daily | stable 400 d | 165.5 K (k3) | 18 | 4 |
| `iso300` | reference isothermal 300 K | stable 365 d | **186.1 K** (k2) | **0** | 0 |
| `no_orog` | `NoOrography` | stable 365 d | **186.1 K** (k3) | **0** | 0 |

The same at T128L8 (400 days each), where the default is stable but shows the same cold spots:

| run | lowest T_min | max. points < 170 K | days n>0 |
|---|---|---|---|
| T128 default | 142.4 K (k1) | 375 | 127 |
| T128 `update30` | 166.3 K (k3) | 54 | 8 |
| T128 `iso300` | **192.3 K** (k2) | **0** | 0 |

Interpretation:

- Survival alone is a weak test, as default runs blew up anywhere between day 181 and 358 (`dt_half` and
  `diff_2h` survived with cold spots as severe as the default). The cold-spot statistics are the stronger
  signal.
- **The shape of the reference matters, not its temperature.** `plus30` matches the top-layer mean
  temperature but keeps the tropospheric lapse rate between k1 and k2, and still blows up. The isothermal
  300 K reference and the refreshed (current) profile both have at least the model's static stability at the
  top, and both remove or strongly reduce the cold spots. In σ coordinates the stability at the k1/k2
  interface, κT/σ − ∂T/∂σ, is ≈280 K for the frozen reference, ≈350 K for `plus30`, ≈570 K for the model
  climate and ≈690 K for `iso300`.
- **The Antarctic orography is the forcing, but through its large scale.** Without orography there are no
  cold spots. Smoothing the T256 orography to roughly T128 does not help, so the sharper T256 slopes are not
  why T256 blows up and T128 does not.
- T128 has the same cold spots, at lower amplitude, and the same cures remove them. T256 is the same
  phenomenon across a stability threshold rather than a separate problem.
- `iso300` has stronger subtropical jets (|u| up to 129 m/s at k2–k3, ~40°N/S, against ~90–110 m/s
  otherwise); global means are unchanged (T_top 232.4 vs 231.5 K, T_sfc 275.2 vs 276.0 K). This needs a look
  before choosing it as the fix.

### 5. The stale reference is not linearly unstable

`linear_stability.jl` analyses the semi-implicit leapfrog exactly as implemented, linearized about a resting,
horizontally uniform atmosphere, per spherical-harmonic degree: the reference operators (R, U_ref, L_ref, W)
implicitly (backward, α = 1), the residuals (L_act − L_ref)·D and (U_act − U_ref)·lnpₛ explicitly at time i,
the linear terms at i−1 as in the model, implicit hyperdiffusion, and the Robert–Asselin–Williams filter. All
operators come from the model's own `ImplicitPrimitiveEquation` and `HyperDiffusion`. The largest eigenvalue
modulus of the two-step map is the growth per time step.

- The analysis reproduces the known instability: with T_ref = 100 K or 80 K (below T_act/2, the classical
  condition for the centred scheme), up to |λ| = 1.4–2.6 per step, nearly all degrees unstable.
- For every realistic reference (the frozen initial profile, `plus30`, `iso300`, the climate profile, even
  150 K), and for the global-mean climate, a cold-spot column (129–218 K) and a warm polar-stratosphere column
  as actual profiles, the scheme is **stable** at T128 and T256, for α = 1 and α = 0.5, with and without
  diffusion and filter (max |λ| = 1, the neutral balanced mode).

So the stale reference does not cause the classical semi-implicit instability. Its effect must come through
what the linear model leaves out: the strong flow over the Antarctic plateau (`no_orog` cures it), horizontal
variation of the basic state, or differences in how strongly the backward-implicit scheme damps particular
vertical modes.

A first look at the noise contradicts a simple "the stale reference makes more noise" picture. From the same
snapshot, 10-day mean leapfrog time-level difference, stale divided by refreshed reference:

| restart from | D k1 | D k2 | D k3 | D k4 | D k8 | T k1 | T k4 | T k8 |
|---|---|---|---|---|---|---|---|---|
| day 240 | 0.35 | 0.42 | 0.66 | 0.76 | 0.94 | 0.72 | 0.74 | 0.95 |
| day 270 | 0.61 | 0.67 | 0.79 | 0.98 | 1.21 | 0.90 | 0.87 | 1.12 |

With the stale reference, the time-step-scale divergence aloft is *smaller*, not larger. One reading, not yet
tested: the stale reference suppresses the divergent (gravity-wave) adjustment in the upper layers, so cold
anomalies over the plateau are not removed by adiabatic descent and can accumulate.

The linear model does not support that reading either (`linear_modes.jl`, T256, actual = climate profile).
The least-damped oscillatory mode trapped in the top two layers (> 50% of |D|² in k1–k2) is damped almost
equally with every reference: per day 0.33 (stale) vs 0.34 (climate) vs 0.30 (`iso300`) at l = 10,
4.0e-5 vs 4.9e-5 vs 2.0e-5 at l = 30, and by more than 15 orders of magnitude per day for l ≥ 60. So neither
linear growth nor differential linear damping separates the stale from the refreshed reference. The mechanism
must be nonlinear or tied to the horizontally varying state over the plateau. Examples are the explicit
residual R_d(T − T_ref)∇ln pₛ over the steep plateau edges, and the vertical advection residual
(ΔT_local − ΔT_ref)·σ̇ under the large σ-velocities of flow across the slopes. The cold-spot columns
(e.g. 129, 160, 168, 181, … K on day 277) cool with height at all levels, like the stale reference
(200, 226, 243, … K), unlike the model's global mean, which has a top-level inversion. A temperature tendency
budget at the cold spots, run from the same snapshot with the stale and the refreshed reference, would show
which term differs.

### Ruled out or unlikely so far

- **Float32 precision**: Float64 also blows up.
- **GPU-specific bug**: the CPU restart blows up at the same time as the GPU restart (day 279.75).
- **Reactant, TF32, DFT**: the plain CUDA.jl and FFT/Legendre setup blows up too (earlier investigation).
- **A convectively unstable top layer**: θ_k1 − θ_k2 > +30 K everywhere.
- **A recent regression in the implicit solver re-initialization**: #1297 only moved the `reinitialize!` call.
  The reference profile has been frozen at the initial conditions in `run!` before as well.
- **A reference that is too cold** (classical semi-implicit instability): `plus30` still blows up, and the
  linear analysis is stable for all realistic references.
- **The sharper T256 orography**: `orog_T128` still blows up.
- **The time step**: `dt_half` has cold spots as severe as the default.

### 7. Released versions: the instability is not a regression

> bisect previous released version whether this issue already existed before
>
> exclude 0.18.1 from the tests it's too old and very different in internal structure

Each release is tested from its own git worktree (`versions/wt_<tag>`) with its own Julia environment
(`versions/env_<tag>`, CUDA.jl 6.4.x, all five packages from the worktree via `[sources]`). The test runs the
default `PrimitiveWetModel` at T256L8 (`truncation = 256`, or `trunc = 255` before v0.22) with a single `run!`
over 730 days, as a user would. Chunked `run!` calls would rebuild the semi-implicit operators every chunk,
and the stale reference profile is a lead suspect. A callback prints daily diagnostics and aborts at the first
non-finite value (`versions/release_run.jl`, `versions/job_release.sh`, logs in `versions/logs/`). The compute
nodes have no internet, so the assets of the older versions were prefetched on the login node first.

| release | date | result | signature before the blow-up |
|---|---|---|---|
| v0.23.0 | 2026-10-08 | **NaN at day 297** | |
| v0.22.1 | 2026-08-21 | **NaN at day 155** | k1 cold spots at 77–86°S, T_min 156 K, ~200 points < 170 K, \|v\| ≈ 80 m/s |
| v0.21.1 | 2026-06-15 | **NaN at day 45** | *different*: no cold spots (T_min ≈ 188 K over Greenland, 0 points < 170 K), max \|u\| creeping to ~100 m/s, then NaN within a day |
| v0.20.3 | 2026-06-03 | **NaN at day 69** | abrupt like v0.21.1: no point < 170 K, T_min ≈ 189 K over Greenland, max \|u\| 96 → 102 m/s in the last 3 days |
| v0.19.0 | 2026-04-23 | **NaN at day 54** | abrupt like v0.21.1: no point < 170 K on any day, max \|u\| ≈ 95 m/s; in the last 3 days T_min fell 185 → 174.5 K at 76–80°S, k1 |

v0.19.0 and v0.20.3 run ~10× slower on the GPU (~1.4 simulated days per minute), so within the 6 h limit
they cover ~500 days, which spans the first-year window of all blow-ups so far (day 155–358).

v0.21.1 fails through an older, known instability. The v0.22.0 changelog lists #1176: "Fix `SpeedLimitDrag`
pointing along the quadrant diagonal …; default `drag` raised 4e-7 -> 3e-6 so the limiter actually binds,
fixes previous high-resolution instabilities". This earlier failure masks whether the cold-spot instability
existed before v0.22. To unmask it, v0.19.0, v0.20.3 and v0.21.1 are rerun with only the drag coefficient
raised to the post-#1176 default, `SpeedLimitDrag(spectral_grid; drag = 3e-6)` (option `DRAG=3e-6`, 12 h
limit, tags `rel_<version>_drag3e-6`). The direction fix of #1176 is not included.

| release + `drag = 3e-6` | result | signature |
|---|---|---|
| v0.21.1 | survived the first southern winter, running (day 401 at 14:04) | **full cold-spot signature** over days 150–365: lowest T_min 134.4 K (day 262, 81.7°S, k1), up to 999 points < 170 K, 125 days with cold points. That is as strong as v0.23.0 over the same days (128.5 K, 1028 points, 115 days), which blew up at day 297 |
| v0.20.3 | survived the first southern winter, running (day 295 at 15:01) | **full cold-spot signature**: lowest T_min 143.9 K (day 288, 77.9°S, k2), up to 647 points < 170 K, 56 days with cold points |
| v0.19.0 | **NaN at day 118** | abrupt, as without the drag change: no point < 170 K, max \|u\| ≈ 100 m/s; the higher drag coefficient alone does not remove the early v0.19 instability |

**Conclusion.** The cold-spot instability is not a regression. It was already present in v0.21.1
(June 2026), and is developing in v0.20.3, but there it was hidden behind the older high-resolution
instability that #1176 fixed in v0.22.0. Surviving one winter with the full signature (v0.21.1 + drag) is
luck, like `dt_half` and `diff_2h` (finding 4). No release from v0.19.0 on is free of a T256L8 instability,
so there is nothing to bisect. This fits the cause: the reference profile has been frozen at the initial
conditions within a `run!` in all these versions. v0.19.0 still fails early, from an abrupt instability that
the higher drag coefficient alone does not remove (possibly the drag-direction bug, also fixed in #1176 but not
ported here). That one is not investigated further.

### 6. The cures hold in a 2-year ensemble

2-year (730-day) T256L8 runs, each from the default initial conditions plus 0.01 K white noise in temperature
(`PERTURB_SEED` 1 and 2), on one H100 each (Slurm 2585888–2585893). Statistics from day 150 to the end:

| run | outcome | lowest T_min | max. points < 170 K | days n>0 | max \|u\| |
|---|---|---|---|---|---|
| default, seed 1 | **NaN day 257** | 141.9 K (k1) | 1305 | 27 | 106 m/s |
| default, seed 2 | **NaN day 181** | 166.2 K (k1) | 44 | 3 | 105 m/s |
| `update30`, seed 1 | stable 730 d | 166.8 K (k3) | 15 | 9 | 108 m/s |
| `update30`, seed 2 | stable 730 d | 161.5 K (k3) | 43 | 11 | 118 m/s |
| `iso300`, seed 1 | stable 730 d | 182.7 K (k2) | 0 | 0 | 125 m/s |
| `iso300`, seed 2 | stable 730 d | 179.5 K (k3) | 0 | 0 | 128 m/s |

The default has now blown up in 8 of 8 runs, the latest at day 358. Both cures ran 2 years in both
realizations, through two southern winters. `iso300` removes the cold spots completely, but again shows
subtropical jets of 125–128 m/s. `update30` keeps the jets close to the default (108–118 m/s) and leaves a
few mild cold points (≥ 161 K, in k3, not k1).

## Next steps

1. Collect the fresh CPU run (Slurm 2582905), which reaches the blow-up window ~13 h after this update.
2. Pin down the mechanism. The linear analysis rules out the classical instability, so look at what it
   leaves out: the damping of individual vertical modes by the backward-implicit scheme with each reference,
   and the gravity-wave noise over the Antarctic plateau. The leapfrog time-level difference (`noiseD_k*`,
   `noiseT_k*` in the CSVs of the `N240_*` and `N270_*` restarts) compares the noise with the stale and the
   refreshed reference from the same state.
3. Write the fix plan for sign-off. Candidates are a periodic re-initialization of the implicit operators
   from the current global-mean profile, a fixed isothermal reference profile, or a climatological
   reference profile instead of the initial-condition one. `iso300` changes the jets and needs checking.
   Each candidate needs checking for differentiability (`reinitialize!` is inactive for Enzyme) and Reactant.
   **Drafted** (periodic re-initialization, default daily): plan
   `docs/dev/2026-10/2026-10-19-implicit-reference-update.md` and code on the local branch
   `mg/implicit-reference-update` (worktree `/p/tmp/maxgelbr/code/claude_t256_instability_2026-10-19/fix_worktree`),
   uncommitted. With the draft, the default T256L8 model is stable for 730 days in 3 of 3 runs, and T128L8
   loses its cold spots.
