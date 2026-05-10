# Dynamical Core GPU Optimization Plan

## Goal

Make the primitive-equation dynamical core (`dynamics_tendencies!(vars, lf, model::PrimitiveEquation)`) faster on GPU while not regressing on CPU. The two main levers identified by the user:

1. **Reduce or fuse `transform!` calls** — they are already optimized internally, but each call has launch overhead and GPU benefits from batching multiple transforms into one fewer-but-bigger call.
2. **Other GPU-friendly optimizations** — kernel fusion, reduced memory traffic, fewer kernel launches.

The model state and high-level algorithm are described in [`docs/src/primitiveequation.md`](docs/src/primitiveequation.md). The function being analyzed is at [`SpeedyWeather/src/dynamics/tendencies.jl:48-110`](SpeedyWeather/src/dynamics/tendencies.jl#L48-L110).

## General rule: keep variable-level readability via views

Whenever we fuse multiple variables into a single backing buffer (a stacked `Field`, a wider `LowerTriangularArray`, etc.) for transform batching or kernel fusion, the rest of the code must keep accessing them through **named `view`s** — not raw `[:, k:k+nlayers]`-style slices sprinkled across kernels.

Concretely:
- Allocate one fused parent buffer in `vars.scratch` (or wherever appropriate).
- Expose each logical variable as a named `view` of that parent: `temp_tend_grid = view(parent, :, 1:nlayers)`, `uT_grid = view(parent, :, nlayers+1:2nlayers)`, etc.
- Pass those named views into the kernels and the spectral post-ops (`curl!`, `divergence!`, `∇²!`) just like today's code passes the standalone fields.
- Only the single batched `transform!` call references the parent buffer directly.

This way the GPU gets one fat transform but the source code still reads as if `temp_tend_grid`, `uT_grid`, `vT_grid` are independent arrays. Apply this rule throughout Phases B and C.

---

## 1. Inventory: transforms per call to `dynamics_tendencies!` (PrimitiveWet, no tracers)

Counting `transform!` calls inside the primitive-equation `dynamics_tendencies!`:

| # | Step | File:line | Direction | Layers | Notes |
|---|------|-----------|-----------|--------|-------|
| T1 | `pressure_gradient_flux!` | [tendencies.jl:130](SpeedyWeather/src/dynamics/tendencies.jl#L130) | spec→grid | 1 (2D) | `dpres_dx_spec → dpres_dx`, unscale_coslat |
| T2 | `pressure_gradient_flux!` | [tendencies.jl:131](SpeedyWeather/src/dynamics/tendencies.jl#L131) | spec→grid | 1 (2D) | `dpres_dy_spec → dpres_dy`, unscale_coslat |
| T3 | `surface_pressure_tendency!` | [tendencies.jl:332](SpeedyWeather/src/dynamics/tendencies.jl#L332) | grid→spec | 1 (2D) | `pres_tend_grid → ūv̄∇lnpₛ` |
| T4 | `vordiv_tendencies!` | [tendencies.jl:476](SpeedyWeather/src/dynamics/tendencies.jl#L476) | grid→spec | nlayers | `u_tend_grid → u_tend` |
| T5 | `vordiv_tendencies!` | [tendencies.jl:477](SpeedyWeather/src/dynamics/tendencies.jl#L477) | grid→spec | nlayers | `v_tend_grid → v_tend` |
| T6 | `temperature_tendency!` | [tendencies.jl:610](SpeedyWeather/src/dynamics/tendencies.jl#L610) | grid→spec | nlayers | `temp_tend_grid → temp_tend` |
| T7 | `flux_divergence!` (T) | [tendencies.jl:774](SpeedyWeather/src/dynamics/tendencies.jl#L774) | grid→spec | nlayers | `uT_grid → uT` |
| T8 | `flux_divergence!` (T) | [tendencies.jl:775](SpeedyWeather/src/dynamics/tendencies.jl#L775) | grid→spec | nlayers | `vT_grid → vT` |
| T9 | `horizontal_advection!` (q) | [tendencies.jl:718](SpeedyWeather/src/dynamics/tendencies.jl#L718) | grid→spec | nlayers | `q_tend_grid → q_tend` |
| T10 | `flux_divergence!` (q) | [tendencies.jl:774](SpeedyWeather/src/dynamics/tendencies.jl#L774) | grid→spec | nlayers | `uq_grid → uq` |
| T11 | `flux_divergence!` (q) | [tendencies.jl:775](SpeedyWeather/src/dynamics/tendencies.jl#L775) | grid→spec | nlayers | `vq_grid → vq` |
| T12 | `bernoulli_potential!` | [tendencies.jl:961](SpeedyWeather/src/dynamics/tendencies.jl#L961) | grid→spec | nlayers | `KE_grid → KE` |

**Totals per dynamics_tendencies!**
- spec→grid: 2 (both 2D, single layer)
- grid→spec: 1 × 2D + 9 × 3D = 10
- + 3 more grid→spec per active tracer (T*div, u*T, v*T flavour — ~3 per tracer)

In addition, the surrounding driver `time_step!` does another spec→grid for the new state ([transform.jl:158-175](SpeedyWeather/src/time_stepping/transform.jl#L158-L175)) — these are already needed and not in scope here.

The Legendre transform is already batched in the vertical (one `transform!` of an `(lm × N)` array runs all `N` slices at once; see [`SpeedyTransforms/src/legendre.jl`](SpeedyTransforms/src/legendre.jl)). So conceptually, fusion = stacking variables into a single LTA's "vertical" axis and doing one bigger transform.

---

## 2. Optimization opportunities

Sorted by expected GPU impact, with implementation cost.

### 2.1 [HIGH] Fuse the three grid→spec transforms inside `temperature_tendency!`

Currently for temperature we do, in order:
1. Build `temp_tend_grid` (T'D + adiabatic, in a kernel)
2. `transform!(temp_tend, temp_tend_grid, S)`         ← T6
3. Build `uT_grid`, `vT_grid` (in a kernel inside `flux_divergence!`)
4. `transform!(uT, uT_grid, S)`                        ← T7
5. `transform!(vT, vT_grid, S)`                        ← T8
6. `divergence!(temp_tend, uT, vT, S, add=true, flipsign=true)` (spectral)

**Optimization**: pack `temp_tend_grid`, `uT_grid`, `vT_grid` into a single 3-batch `Field` (or call the transform on a `view` of the 4D buffer of shape `(npoints × nlayers × 3)`), do **one** fused grid→spec, then `divergence!` on the appropriate spectral slices and `add!` the temp slice into `temp_tend`.

Trades 3 transforms → 1 transform of size `3·nlayers`. The Legendre matmul cost is roughly the same (it scales linearly), but FFT plan setup, scratch reuse, and kernel launches are amortized once across the batch. On GPU that's typically a 2–3× wallclock improvement on this stage.

**Files touched**: [`tendencies.jl`](SpeedyWeather/src/dynamics/tendencies.jl) — `temperature_tendency!`, `flux_divergence!`. Likely needs a small helper that allocates a 3-layer scratch field, or extends `vars.scratch.grid` with a "stacked" buffer.

**Risk**: requires the FFT/Legendre plans in `SpectralTransform` to support `nlayers` larger than the model's `nlayers`. Currently `S.nlayers` is allocated to the model's `nlayers`; the scratch arrays (`g_north`, `g_south`) and column scratch are sized to `S.nlayers`. We need to grow them to `≥ 3·nlayers` (or whatever the largest fused batch is) at construction. This is a one-time setting; the change is local to `SpectralTransform` constructor.

### 2.2 [HIGH] Fuse the three grid→spec transforms inside `humidity_tendency!`

Identical structure to 2.1 but for humidity (T9, T10, T11). Same approach, same buffer-resize requirement. Apply once 2.1 is shaken out.

Combined with 2.1 this turns 6 separate `transform!` calls into 2 (or 1 if 2.4 is also done).

### 2.3 [HIGH] Fuse `vordiv_tendencies!` u_tend + v_tend transforms

T4 and T5 are independent, both grid→spec, identical layout. Pack them in a 2-batch field, call one `transform!`, then run `curl!` and `divergence!` on the relevant slice pair.

This is a small but cheap win: 2 transforms → 1.

### 2.4 [MEDIUM-HIGH] Single mega-batch of all grid→spec transforms in `dynamics_tendencies!`

Detailed in §6 below — this is now the leading candidate and gets its own section.

### 2.5 [HIGH] Fuse the two spec→grid transforms in `pressure_gradient_flux!`

T1 and T2 are two single-layer (2D) spec→grid transforms with `unscale_coslat=true`. Currently the most expensive overhead-per-flop calls in the whole tendency function (very small payload, full plan/launch cost).

Pack `dpres_dx_spec` and `dpres_dy_spec` into a 2-layer LTA, do one transform, take views.

### 2.6 [MEDIUM] Kernel fusion: combine flux preparation kernels

Currently `flux_divergence!` launches `_flux_divergence_kernel!` to compute `uA_grid`, `vA_grid`. For temperature we also have a separate kernel for `T'*D + adiabatic`. These all read `u`, `v`, `coslat⁻¹`, `whichring`, and a tracer `A`. Fuse them into one kernel:

```julia
@kernel function _temp_grid_prep_kernel!(temp_tend_grid, uT_grid, vT_grid,
                                         u, v, T, divergence, humid, ...)
    ij, k = @index(Global, NTuple)
    j = whichring[ij]; coslat⁻¹j = coslat⁻¹[j]
    Tᵥ = ...; T_anom = T[ij,k]
    # adiabatic + T'*D contribution
    temp_tend_grid[ij,k] += T_anom * div[ij,k] + κ*Tᵥ*(...)
    # flux prep
    Tcoslat = T_anom * coslat⁻¹j
    uT_grid[ij,k] = u[ij,k] * Tcoslat
    vT_grid[ij,k] = v[ij,k] * Tcoslat
end
```

Reduces 2 kernel launches to 1 and re-uses the cached load of `u, v, T` per thread. Same idea for humidity and tracers. Same idea for u_tend/v_tend + KE preparation: the bernoulli kernel reads the same `u, v` already loaded by `vordiv_tendencies` — they could share a kernel.

### 2.7 [MEDIUM] Avoid the temperature anomaly add/subtract on grid

Lines [tendencies.jl:71](SpeedyWeather/src/dynamics/tendencies.jl#L71) and [tendencies.jl:107](SpeedyWeather/src/dynamics/tendencies.jl#L107) do:

```julia
vars.grid.temperature.data .-= implicit.temp_profile'
... (lots of work that uses temperature as anomaly) ...
vars.grid.temperature.data .+= implicit.temp_profile'
```

This means **two full passes over the temperature field** (`npoints × nlayers` reads+writes each) per timestep just to subtract/add a small per-layer profile. On GPU this is bandwidth-bound and avoidable.

Better: pass `temp_profile` into the kernels that need the anomaly (they already get `Tₖ` for virtual temperature) and do `temp_grid[ij,k] - Tₖ[k]` on the fly. The two kernels that consume the anomaly are `_vordiv_tendencies_kernel!` ([tendencies.jl:484-513](SpeedyWeather/src/dynamics/tendencies.jl#L484-L513)) and `_temperature_tendency_kernel!` ([tendencies.jl:617-650](SpeedyWeather/src/dynamics/tendencies.jl#L617-L650)) plus the flux-prep step (uses `temp` as `T'`, the anomaly form).

Drop both `.-=` / `.+=` calls. Saves 2 full grid passes over T.

### 2.8 [MEDIUM] Avoid re-reading u_mean*∇lnp in `vertical_velocity!`

[tendencies.jl:358-360](SpeedyWeather/src/dynamics/tendencies.jl#L358-L360) computes `ūv̄∇lnp = u_mean*dpres_dx + v_mean*dpres_dy` on the grid — but the same quantity is already computed earlier in `surface_pressure_tendency!` ([tendencies.jl:329](SpeedyWeather/src/dynamics/tendencies.jl#L329)) accumulated into `pres_tend_grid`.

Either:
- store `ūv̄∇lnp` in a tiny dedicated scratch and reuse it in both places, or
- fold the computation of `ūv̄∇lnp` (and the `u_mean*dpres_dx + v_mean*dpres_dy` math) into the final iteration of the `_vertical_integration_kernel!` (we already have `u_mean`, `v_mean` accumulated by then).

Either way, avoid the redundant pointwise broadcast in `vertical_velocity!`.

### 2.9 [MEDIUM] Fuse `linear_virtual_temperature!` into the spectral `geopotential!` kernel

`linear_virtual_temperature!` ([virtual_temperature.jl:28-47](SpeedyWeather/src/dynamics/virtual_temperature.jl#L28-L47)) does one full LTA broadcast `Tᵥ = T + (Tₖ μ) q`, just to feed it into the spectral `geopotential!` kernel ([tendencies.jl:74](SpeedyWeather/src/dynamics/tendencies.jl#L74), kernel at [geopotential.jl:153-170](SpeedyWeather/src/dynamics/geopotential.jl#L153-L170)).

Pass `T_spec`, `q_spec`, `Tₖ`, `μ` directly into `geopotential_spectral_kernel!` and compute `Tᵥ` on the fly per `(lm, k)`. Saves a full LTA-sized read+write pass on GPU.

### 2.10 [MEDIUM] Fold `linear_pressure_gradient!` into `bernoulli_potential!`

[tendencies.jl:86](SpeedyWeather/src/dynamics/tendencies.jl#L86): `Φ.data .+= R_dry .* Tₖ' .* lnpₛ.data` (LTA broadcast across nlayers).
[tendencies.jl:101](SpeedyWeather/src/dynamics/tendencies.jl#L101): `bernoulli_potential!` does `bernoulli += geopot` (LTA add).

Both are LTA passes at the same time in the algorithm. The linear pressure gradient term `R_dry*Tₖ*lnpₛ` is added under the same `−∇²` as the geopotential and the kinetic energy, by design (see `primitiveequation.md` §"Pressure gradient"). So just compute `bernoulli += geopot + R_dry*Tₖ*lnpₛ` in one kernel/broadcast — saves a full Φ pass.

### 2.11 [LOW] `surface_pressure_tendency` mass-conservation gate

[tendencies.jl:337](SpeedyWeather/src/dynamics/tendencies.jl#L337): `pres_tend.data[1:1] .= 0`. This is a length-1 slice operation; on GPU this triggers a minimum-size kernel launch. Replace with a scalar `getindex/setindex!` (still launches but cheaper) or fold into a host-side small write. Verifying this is actually a measurable cost requires a profile.

### 2.12 [LOW] `reset_tendencies!` allocations

[tendencies.jl:1007-1024](SpeedyWeather/src/dynamics/tendencies.jl#L1007-L1024) walks the tendency NamedTuple recursively. Each `fill!` is a kernel on GPU. Currently this fires for every leaf. Two ideas: (a) pre-collect the array list once at construction time, (b) check whether some "tendency" arrays are written by an `=` rather than `+=` and don't need zeroing. Probably a small win; profile first.

### 2.13 [LOW] Pre-allocated temp profile broadcast

The TODO in [tendencies.jl:70](SpeedyWeather/src/dynamics/tendencies.jl#L70) and [virtual_temperature.jl:44](SpeedyWeather/src/dynamics/virtual_temperature.jl#L44) flag broadcast issues with LTA + per-layer vector. Once fixed (in LowerTriangularArrays), several of the "use `.data` to break broadcast" lines can be cleaned up. Mostly hygiene/maintainability rather than perf, but enables more KA-style fused operations.

### 2.14 [LOW] Remove redundant `vars.scratch.grid.a .= 0` for dry models

[tendencies.jl:457](SpeedyWeather/src/dynamics/tendencies.jl#L457) and [tendencies.jl:591](SpeedyWeather/src/dynamics/tendencies.jl#L591) zero the scratch only when humidity is missing. Two zero-fills per timestep on a 3D grid. We could share a single "zero scratch" per timestep, or read a sentinel zero in the kernels for dry models. Mostly relevant for `PrimitiveDry`.

---

## 3. Proposed implementation order (incremental, each PR-able)

Each step keeps the test suite passing and is independently measurable.

**Phase A — easy fusion wins (no API change to `SpectralTransform`)**

A1. **2.5** Fuse `dpres_dx`, `dpres_dy` spec→grid transforms (single 2-layer batched transform). Local change in `pressure_gradient_flux!`.

A2. **2.7** Drop the temperature anomaly add/subtract; fold `-Tₖ` into the consumer kernels (`_vordiv_tendencies_kernel!`, `_temperature_tendency_kernel!`, flux prep).

A3. **2.10** Fold `linear_pressure_gradient!` into the bernoulli combination.

A4. **2.8** Reuse `u_mean*dpres_dx + v_mean*dpres_dy` between `surface_pressure_tendency!` and `vertical_velocity!`.

A5. **2.9** Fold linear virtual temperature into the spectral geopotential kernel.

**Phase B — multi-variable transform batching (requires growing scratch)**

This phase requires `SpectralTransform` to allocate scratch (`g_north`, `g_south`, column scratch) sized to `K · nlayers`, where K is the largest fused batch (≈9 for PrimitiveWet, 6 for PrimitiveDry).

B1. **2.3** Fuse `u_tend`/`v_tend` in `vordiv_tendencies!` into a 2-batch transform.

B2. **2.1** Fuse temperature: `temp_tend`, `uT`, `vT` into a 3-batch transform.

B3. **2.2** Fuse humidity: `q_tend`, `uq`, `vq` into a 3-batch transform.

B4. **2.4** Mega-batch: pack everything (u_tend, v_tend, temp_tend, uT, vT, q_tend, uq, vq, KE) into one 9-batch transform. Re-arrange `dynamics_tendencies!` so all grid-space prep runs first, then a single transform, then all spectral ops.

**Phase C — kernel fusion (after batching is stable)**

C1. **2.6** Fuse temp tendency kernel + temp flux prep into a single kernel; same for humidity and tracers; same for vordiv + bernoulli (KE).

C2. **2.13** Clean up LTA broadcast hacks once LTA broadcast is fixed.

C3. **2.11**, **2.12**, **2.14** based on profiling results.

---

## 4. Practical notes

### 4.1 Scratch sizing for batched transforms

In [`SpeedyTransforms/src/spectral_transform.jl`](SpeedyTransforms/src/spectral_transform.jl) the scratch memory (`ScratchMemory`, `ColumnScratchMemory`) is allocated with `S.nlayers`. For Phase B we need to over-allocate. Options:

- Add a `transform_batch::Int` kwarg to `SpectralGrid` / `SpectralTransform` that defaults to 1 and is set ≥ 9 when used by the primitive-equation dynamical core.
- Or auto-size by inspecting the model: when constructing the model, set `transform_batch = max(needed_for_dycore, 1)`.

Memory cost for batch=9: scratch arrays (north/south) currently take ~`O(nfreq_max × nlayers × nlat_half × sizeof(Complex{NF}))`. At T31 / 8 layers / 16 nlat_half / Float32 that's ~few MB; at 9× still negligible (~tens of MB). Not a problem.

### 4.2 Fused buffers as views vs. separate fields

Two implementation styles for the fused grid-space buffer:

- **(a) Stacked field**: one `Field` of size `(npoints × (K·nlayers))` and named `view(...)`s into it.
- **(b) Tuple of fields with stacked LTA**: keep grid-space `Field`s separate (one per variable), but allocate the spectral output as one large LTA. The `transform!` API operates on a single `Field`–`LTA` pair, so style (a) is simpler from the API perspective.

Style (a) simplifies the transform call but requires that each kernel writing into a "named slice" use the proper view. In practice, kernels already accept `AbstractArray`, so passing a view is fine; the only worry is performance regression from view indexing. KernelAbstractions usually handles this well, but worth a microbenchmark.

### 4.3 Correctness checks

- Run `SpeedyWeather/test/dynamics/` after each step.
- Bit-reproducibility: prefer keeping the order of floating-point reductions identical to the current implementation where possible. Some fusions (e.g. 2.10) change associativity — record this in the changelog and verify that the existing tolerance-based tests still pass.
- For Phase B, carefully verify that all `vars.dynamics.*` and `vars.scratch.*` consumers see fresh data when the algorithm reorders.

### 4.4 Benchmarking

Use the existing benchmark suite ([`SpeedyWeather/benchmark/`](SpeedyWeather/benchmark/), see `manual_benchmarking.jl`) before/after each phase. Compare GPU and CPU separately; the goal is GPU speedup, but CPU should not regress more than ~10%. The CLAUDE.md note: "deviations of +/- 20% are normal" — so target ≥20% sustained speedup on GPU per phase to call it a real win, and watch CPU within 20% of baseline.

Current PrimitiveDry / PrimitiveWet baselines should be recorded in `benchmark/README.md` before starting Phase A.

### 4.5 Things explicitly NOT in scope

- Changing the spectral transform algorithm itself (FFT/Legendre).
- Vertical advection scheme changes (WENO etc.) — already kernelized.
- Implicit correction (`implicit_correction!`) — separate from `dynamics_tendencies!`.
- Parameterizations.
- Output / NetCDF.

---

## 5. Open questions — answered

**1. Does the fused 4D → 4D transform work as-is?** **Yes, by construction.** The transform pipeline is layer-count-agnostic:

- `transform!(field, coeffs, scratch, S)` does no per-layer dispatch — it calls `_fourier!` and `_legendre!`.
- `_fourier!` reads `nlayers = size(field, 2)` ([fourier.jl:36, 67, 133, 171](SpeedyTransforms/src/fourier.jl)).
- `_legendre!` reads `nlayers = axes(specs, 2)` ([legendre.jl:51, 138](SpeedyTransforms/src/legendre.jl#L51)).
- `ismatching` allows `nlayers <= S.nlayers` ([spectral_transform.jl:297, 307](SpeedyTransforms/src/spectral_transform.jl#L297)).

So passing a `(lm × K·model_nlayers)` LTA / `(npoints × K·model_nlayers)` Field through `transform!` Just Works, **provided** `S.nlayers ≥ K·model_nlayers` at construction. That's exactly what S1' (the `transform_batch` kwarg) handles. No additional plumbing needed.

A microbench is still worthwhile to confirm performance scaling is roughly linear in `K`, but no design risk remains.

**2. Is the spec→grid 2D fusion (item 2.5) worth it?** **Cannot answer from code alone — pure benchmark question.** Without measuring on actual hardware, the rough analytical picture:

- Two single-layer transforms ≈ 2× per-launch overhead, 2× work.
- Fusing into one 2-layer transform ≈ 1× per-launch overhead, 2× work.
- On GPU: launch overhead is typically tens of μs while a T31 single-layer transform itself is also tens of μs → ~30–50 % stage-local speedup is plausible.
- On CPU: launch overhead is negligible; gain ≈ 0 %.

Whether the absolute saving matters depends on this stage's profile share. Action: **defer 2.5 until Phase A is profiled.** If `pressure_gradient_flux!` doesn't show up as a hot spot, skip it.

**3. Where else is `vars.scratch.grid.a/b` used?** Mapped out across the whole codebase via grep:

| File | Use | Phase |
|------|-----|-------|
| [`parameterizations/vertical_diffusion.jl:114, 175, 264, 276, 278`](SpeedyWeather/src/parameterizations/vertical_diffusion.jl#L114) | dry static energy, K, Ri, humidity stand-in | parameterizations |
| [`parameterizations/convection.jl:58-59, 296`](SpeedyWeather/src/parameterizations/convection.jl#L58) | reference T/q profiles | parameterizations |
| [`parameterizations/radiation/{short,long}wave_transmissivity.jl`](SpeedyWeather/src/parameterizations/radiation/shortwave_transmissivity.jl) | transmissivity scratch `t` | parameterizations |
| [`dynamics/tendencies.jl:457-458, 591-592`](SpeedyWeather/src/dynamics/tendencies.jl#L457) | zero-fill stand-in for missing humidity in dry model | dycore (`vordiv_tendencies!`, `temperature_tendency!`) |
| [`dynamics/tendencies.jl:763-764`](SpeedyWeather/src/dynamics/tendencies.jl#L763) | `uA_grid`, `vA_grid` in `flux_divergence!` | dycore |
| [`dynamics/tendencies.jl:916, 944`](SpeedyWeather/src/dynamics/tendencies.jl#L916) | `bernoulli_grid` | dycore (and ShallowWater) |
| [`dynamics/geopotential.jl:67-68`](SpeedyWeather/src/dynamics/geopotential.jl#L67) | zero-fill stand-in for missing humidity in spectral geopotential! | dycore |

Pattern: every consumer follows write-before-read within a single pass. No persistent state lives in `scratch.grid.a/b`. Parameterizations and the dycore use them strictly sequentially within a timestep — no race.

**Implication for the plan**: after Phase B the dycore's uses (rows 4–7 of the table) move onto the fused parent buffer and stop touching `scratch.grid.a/b`. Parameterizations keep using them. We don't delete the scratch buffers; we just stop scribbling on them from the dycore. Easy.

**One subtlety to be careful of:** in S2 we replace `vars.tendencies.grid.{u,v,T,q}` with views of the fused parent buffer. **Parameterizations write to those names** (e.g. [`stochastic_physics.jl:47-49`](SpeedyWeather/src/parameterizations/stochastic_physics.jl#L47), [`convection.jl:46-47`](SpeedyWeather/src/parameterizations/convection.jl#L46), [`vertical_diffusion.jl:95-97`](SpeedyWeather/src/parameterizations/vertical_diffusion.jl#L95), surface fluxes, radiation, large-scale condensation). They need the views to behave identically to standalone Fields under broadcast (`.+=`, `[ij,k] += ...`) and KA kernels. As §6.4 noted, this is already exercised in production and should be transparent — but it is the single highest-value smoke test in S2.

**4. CPU performance regression risk for Phase B.** Concrete picture from looking at `_legendre!`:

- The hot inner loop is `_fused_oddeven_matvec!` ([legendre.jl:78](SpeedyTransforms/src/legendre.jl#L78)) — a per-`(m, j)` matrix-vector multiply where the "vector" dimension is `nlayers`. **Bigger nlayers actually helps SIMD vectorization**, since the inner loop `for k in 1:nlayers` is the SIMD axis.
- Working-set scaling: `g_north`, `g_south` are `(nfreq_max × nlayers × nlat_half) × Complex{NF}`. At T31/L8/K=9 (Float32): ≈ 5 MB per buffer × 2 = 10 MB. Comfortably in L2/L3 cache. At T127/L32/K=9: ≈ 50× larger, ~500 MB → blows out L3 on most CPUs.
- So the regression risk is **resolution-dependent**, not just K-dependent. For the typical SpeedyWeather demo scales (T31–T63, L8–L31), K=9 stays cache-friendly. For higher resolutions where tests run longer, we may need a fallback.

**Decision**: don't pre-commit to a CPU-fallback path. Run the benchmark suite on T31/T63 first; only add a CPU dispatch (split parent into chunks of size ≤ original `nlayers`) if regression > 20 %. The fallback is mechanically trivial (a `for chunk in eachslice(parent, ...) transform!(...) end` loop) — easy to bolt on if needed.

---

## 6. Mega-batch grid→spec transform: detailed design

This is the leading optimization candidate and warrants more careful analysis.

### 6.1 What goes into the fused buffer

For PrimitiveWet, the grid-space arrays that need to be transformed grid→spec inside one call to `dynamics_tendencies!` are (each `(npoints × nlayers)` unless marked):

| Slot | What | Currently lives in |
|------|------|--------------------|
| 1 | `u_tend_grid` | `vars.tendencies.grid.u` (named TendencyVariable) |
| 2 | `v_tend_grid` | `vars.tendencies.grid.v` (named TendencyVariable) |
| 3 | `temp_tend_grid` | `vars.tendencies.grid.temperature` (named TendencyVariable) |
| 4 | `humid_tend_grid` | `vars.tendencies.grid.humidity` (named TendencyVariable, **PrimitiveWet only**) |
| 5 | `uT_grid` | `vars.scratch.grid.a` (re-used scratch, see [tendencies.jl:763](SpeedyWeather/src/dynamics/tendencies.jl#L763)) |
| 6 | `vT_grid` | `vars.scratch.grid.b` (re-used scratch) |
| 7 | `uq_grid` | `vars.scratch.grid.a` (re-used **after** uT/vT have been transformed) |
| 8 | `vq_grid` | `vars.scratch.grid.b` |
| 9 | `KE_grid` = ½(u²+v²) | `vars.scratch.grid.a` (re-used **after** flux_divergence) |
| (10) | `pres_tend_grid` (2D) | `vars.tendencies.grid.pressure` (Grid2D) |

PrimitiveDry: drop slots 4, 7, 8 → **6 multi-layer slots**. PrimitiveWet, no tracers: **9 multi-layer slots**. Each active tracer adds 3 slots (tracer_tend, u·tracer, v·tracer).

The pressure tendency is single-layer. Two options: pack as the last slice of the parent (parent shape `(npoints × (K·nlayers + 1))`), or keep separate. Packing it in saves one tiny transform; the asymmetry is minor.

### 6.2 Where the fused parent buffer lives

User's preference: **no new struct type** — just a regular `Field` allocated as a sibling field, "next to the old grid tendencies".

I suggest a tighter spelling: allocate the parent as a `ScratchVariable` with custom dimensions (a Grid3D-like with `K·nlayers` layers, where K is determined per-model). It is reachable from `vars.tendencies.grid.fused` (or `vars.scratch.grid.tendency_buffer` — bikeshed naming).

The four named tendency entries `vars.tendencies.grid.u/v/temperature/humidity` are **replaced by `field_view`s into the parent**. Auxiliary slots (`uT`, `vT`, `uq`, `vq`, `KE`) are exposed as `field_view`s built locally inside `dynamics_tendencies!` (or as additional `ScratchVariable` views — TBD).

### 6.3 Critical question: how do the named tendency fields become views of the parent?

This is the heart of the design and the biggest risk. Two approaches:

**Approach (i) — patch up after `Variables(model)`**: let the variable system allocate `vars.tendencies.grid.u` etc. normally as standalone Fields, then in a post-construction step:

1. Allocate the parent buffer (size `npoints × K·nlayers`).
2. Build a new NamedTuple where `u`, `v`, `temperature`, `humidity` are replaced by `field_view(parent, :, slot_range)`.
3. Discard the original standalone Field allocations (let GC clean them).
4. Return the modified `Variables`.

Pros: zero changes to the variable system. Just a small "fuse_tendencies!(vars, model)" step run in the model's `Simulation` constructor right after `Variables(model)`.

Cons: the standalone arrays are briefly allocated and thrown away. Tiny one-time waste. Also, replacing fields in an immutable `NamedTuple` needs careful handling — `merge` does this cleanly.

**Approach (ii) — extend the variable allocation system**: add a marker (e.g. a `fused = true` kwarg or a `FusedTendencyVariable` subtype) so that `allocate(...)` knows to allocate them all from a shared parent.

Pros: cleaner semantics, the parent is a first-class concept.

Cons: bigger surface change to the variable system, slightly intrusive.

**Recommendation**: Approach (i). It's local, easily reverted, and keeps the variable system simple. We're carving out an optimization for one specific group of variables, not building a general feature.

### 6.4 Field views: do they work as drop-ins?

`RingGrids` exposes `field_view(field, :, i, args...)` ([field.jl:447](RingGrids/src/field.jl#L447)) which returns a `Field` wrapping `view(field.data, ...)`. So a `field_view(parent, :, k1:k2)` is a `Field{T, 2, SubArray{...}, ...}` — it implements the same `AbstractField{T, 2}` interface as a standalone Field.

This means:
- `vars.tendencies.grid.u .+= ...` (broadcast in parameterizations) — works, same broadcast style as Field.
- KA kernels reading/writing `vars.tendencies.grid.u[ij, k]` — works, indexes through the SubArray.
- `transform!(spec, vars.tendencies.grid.u, S)` — works, but we are NOT going to call this on the named view; we'll call it on the parent.
- `architecture(field_view)` — likely works since it dispatches on the array type and the SubArray retains the parent's GPU array type.

Risk to verify by microbenchmark: contiguous SubArray indexing in KA kernels. The slice is `view(parent, :, k1:k2)` which is contiguous (stride-1 along the leading axis, stride = `npoints` along the second axis — same as a plain matrix). KA kernels should handle this with no overhead, but on Metal/AMDGPU there have historically been edge cases. Add a smoke test on each backend before committing.

### 6.5 Spectral side: does the same trick work for an LTA?

After the mega transform we have a spectral parent `LowerTriangularArray` of shape `(lm × K·nlayers)`. We need named views into this for `curl!`, `divergence!`, `∇²!`.

Need to verify: does `LowerTriangularArrays` support `view(::LTA, :, k1:k2)` returning an LTA-flavored object that the spectral operators accept?

```bash
grep -rn "Base.view\|function view" LowerTriangularArrays/src/
```

If yes, plug-and-play. If not, we need either:
- Add LTA-view support, OR
- Pass a contiguous LTA slice via `LowerTriangularArray(view(parent.data, ...), parent.spectrum)` (similar to `field_view`).

The sub-LTA needs to behave correctly for `curl!`, `divergence!`, `∇²!` which are element-wise spectral kernels parameterized by `(lm, k)` index. Should be straightforward.

**Action item before B4**: confirm LTA views work with `curl!` / `divergence!` / `∇²!` (or add a 3-line `lta_view` helper).

### 6.6 Order-of-operations rewrite

Current order (annotated with what's grid-space (G) vs spectral (S) vs transform (T)):

```
forcing!, drag!                        # writes to vars.tendencies.grid.{u,v,T,q} (G)
pressure_gradient_flux!                # spec-grad + 2 spec→grid (T) + grid mul
linear_virtual_temperature!            # S broadcast: Tᵥ = T + (Tₖμ)q
T_grid -= temp_profile                 # G broadcast
geopotential!(vars, geopot, orog)     # S kernel (Φ from Tᵥ)
vertical_integration!                  # G kernel (u_mean, v_mean, sums)
                                       # S kernel (D_mean)
surface_pressure_tendency!             # G kernel (pres_tend_grid += u_mean*∇x..)
                                       # 1 grid→spec (T)
                                       # S add (pres_tend -= ūv̄∇lnp + D_mean)
vertical_velocity!                     # G kernel
linear_pressure_gradient!              # S broadcast: Φ += R*Tₖ*lnp
vertical_advection!                    # G kernel (modifies u/v/T/q tendencies)
vordiv_tendencies!                     # G kernel + 2 grid→spec (T) + S curl/div
temperature_tendency!                  # G kernel + 1 grid→spec (T)
                                       # flux_divergence!: G kernel + 2 grid→spec (T) + S div
humidity_tendency!                     # horizontal_advection!: G kernel + 1 grid→spec (T)
                                       # flux_divergence!: G kernel + 2 grid→spec (T) + S div
bernoulli_potential!                   # G kernel + 1 grid→spec (T) + S add + S ∇²
tracer_advection!                      # per tracer: ~3 grid→spec
T_grid += temp_profile                 # G broadcast
```

Reordered for mega-batch:

```
PHASE 0 — pre-grid-space prep (unchanged in spirit):
forcing!, drag!
pressure_gradient_flux!                # ideally batched 2-layer spec→grid (item 2.5)
linear_virtual_temperature!            # ideally folded into geopotential kernel (item 2.9)
T_grid -= temp_profile                 # ideally dropped (item 2.7)
geopotential!(vars, geopot, orog)
vertical_integration!                  # ideally also computes pres_tend_grid += ūv̄∇lnp on-the-fly (item 2.8)
surface_pressure_tendency_grid!        # the +=ūv̄∇lnp grid kernel (split out of current surface_pressure_tendency!)
vertical_velocity!                     # reuses ūv̄∇lnp from above
linear_pressure_gradient!              # may be folded into bernoulli combination (item 2.10)
vertical_advection!

PHASE 1 — grid-space tendency prep, all writing into named views of the parent:
fused_grid_prep!:                      # one-or-few kernel(s) writing parent slots:
   slot 1,2: vordiv contribution      → u_tend_grid, v_tend_grid
   slot 3:   T'D + adiabatic           → temp_tend_grid
   slot 4:   q*D                       → humid_tend_grid (PrimitiveWet)
   slot 5,6: uT, vT                    → uT_grid, vT_grid
   slot 7,8: uq, vq                    → uq_grid, vq_grid (PrimitiveWet)
   slot 9:   ½(u²+v²)                  → KE_grid
   slot10 (2D): pres_tend_grid (already prepared in PHASE 0; just shares the transform)

PHASE 2 — single mega-batch transform:
transform!(parent_spec, parent_grid, scratch_memory, S)

PHASE 3 — spectral post-ops on named views:
curl!(vor_tend, u_tend_spec, v_tend_spec, S, add=true)
divergence!(div_tend, u_tend_spec, v_tend_spec, S, add=true)
add!(temp_tend, temp_tend_spec)                  # spectral fold-in of temp_tend slot
divergence!(temp_tend, uT_spec, vT_spec, S, add=true, flipsign=true)
add!(humid_tend, humid_tend_spec)
divergence!(humid_tend, uq_spec, vq_spec, S, add=true, flipsign=true)
bernoulli_combine!(bernoulli_spec, KE_spec, geopot_spec, R*Tₖ*lnpₛ_spec)
∇²!(div_tend, bernoulli_spec, S, add=true, flipsign=true)
pres_tend -= pres_tend_grid_spec + div_mean      # the existing - sign convention
pres_tend.data[1:1] = 0

PHASE 4 — tracers (per active tracer):
either (a) keep separate transforms per tracer (3 each), or
(b) extend the parent buffer to include tracer slots (preferred when there are tracers).

PHASE 5 — restoration:
T_grid += temp_profile                 # ideally dropped (item 2.7)
```

### 6.7 Interaction with existing scratch and "shared" arrays

`vars.scratch.a`, `vars.scratch.b` (spectral 3D) are **also** used outside `dynamics_tendencies!`:
- by `transform!(vars, lf, model::PrimitiveEquation)` for U,V spectral velocities ([transform.jl:135-136](SpeedyWeather/src/time_stepping/transform.jl#L135-L136))
- by `implicit_correction!` ([implicit.jl:369-370](SpeedyWeather/src/time_stepping/implicit.jl#L369-L370))

So we must **not** repurpose them. The new fused parent is an additional allocation, not a replacement of `scratch.a/b`. (`scratch.a/b` will become unused inside `dynamics_tendencies!` after the refactor; their continued use elsewhere is fine.)

`vars.scratch.grid.a`, `vars.scratch.grid.b` (grid 3D) are used inside `dynamics_tendencies!`:
- as `uA_grid`, `vA_grid` in `flux_divergence!` (used for T, q, tracers — overwritten between calls)
- as `bernoulli_grid` in `bernoulli_potential!`
- as a zero-fill for missing humidity in PrimitiveDry ([tendencies.jl:457-458](SpeedyWeather/src/dynamics/tendencies.jl#L457-L458))

After the refactor these scratch buffers are no longer needed by the dycore (the fused parent owns those slots). Verify no other code path writes to them, then they can be removed (or left for now and removed in a cleanup pass).

### 6.8 `reset_tendencies!` interaction

[tendencies.jl:1007-1024](SpeedyWeather/src/dynamics/tendencies.jl#L1007-L1024) walks `vars.tendencies` and `fill!`s every leaf array. With the named tendencies being views of the parent, this would currently:
- `fill!(view of parent for u)` → K1·nlayers slots zeroed
- `fill!(view of parent for v)` → another K1·nlayers slots zeroed
- ... etc.

This is N kernel launches on GPU (N = number of named tendency entries). Replacement: if the parent is also reachable in `vars.tendencies.grid` (e.g. as `vars.tendencies.grid.fused`), `_reset_tendencies!` will fill the parent in one launch — and the named-view fills become redundant zero-fills of already-zero memory.

Two cleaner options:
- (a) Skip the views in `reset_tendencies!` (their slots are zeroed via the parent fill).
- (b) Skip the parent in `reset_tendencies!` (use only the named views — but auxiliary slots like `uT_grid` are write-before-read so don't need zeroing anyway).

Either works. (b) is simpler — keeps `reset_tendencies!` walking the same set of named entries it does today; the parent is **not** in `vars.tendencies.grid` but in `vars.scratch.grid` (or similar) so it's not visited.

This argues for: **put the parent buffer in `vars.scratch.grid.fused`**, not in `vars.tendencies.grid.*`. That places it "next to" the existing scratch grid arrays (which is in spirit "next to the old grid tendencies"). The named views (`u`, `v`, `T`, `q`) live in `vars.tendencies.grid` as today; they're just views into the parent under the hood.

I'll go with this layout in the implementation plan unless you prefer otherwise. Concretely:

```
vars.scratch.grid.fused           ::Field, size (npoints × K·nlayers [+ 1])  ← parent
vars.tendencies.grid.u            = field_view(vars.scratch.grid.fused, :, slot_u)
vars.tendencies.grid.v            = field_view(vars.scratch.grid.fused, :, slot_v)
vars.tendencies.grid.temperature  = field_view(vars.scratch.grid.fused, :, slot_T)
vars.tendencies.grid.humidity     = field_view(vars.scratch.grid.fused, :, slot_q)   # wet only

# auxiliary, declared as locals inside dynamics_tendencies!:
uT_grid   = field_view(vars.scratch.grid.fused, :, slot_uT)
vT_grid   = field_view(vars.scratch.grid.fused, :, slot_vT)
uq_grid   = field_view(vars.scratch.grid.fused, :, slot_uq)   # wet only
vq_grid   = field_view(vars.scratch.grid.fused, :, slot_vq)   # wet only
KE_grid   = field_view(vars.scratch.grid.fused, :, slot_KE)
```

The slot offsets are model-dependent constants computed once at construction.

### 6.9 Sizing the SpectralTransform scratch

The current `SpectralTransform` scratch (`g_north`, `g_south`, column scratch) is sized for `S.nlayers`. The check is `nlayers <= S.nlayers` ([spectral_transform.jl:307](SpeedyTransforms/src/spectral_transform.jl#L307)). We need `S.nlayers ≥ K·model.nlayers` (or `+1` if pressure is packed in).

Where `S.nlayers` is set: when `SpectralTransform` is constructed for the model, currently with the model's `nlayers`. We need to thread through a "transform batch factor" so that the model constructs its `SpectralTransform` with `nlayers = K·model_nlayers + 1`.

Options:
- New kwarg on `SpectralGrid` or `SpectralTransform` (e.g. `transform_batch_layers::Int`). Defaults to 1 (back-compat). Set automatically by primitive-equation models.
- Or: auto-detect by walking the model's tendency variable definitions and computing the largest plausible batch.

I prefer the explicit kwarg with a sensible default, set by the primitive-equation model constructor.

Memory cost: `g_north`, `g_south` are `(nfreq_max × nlayers × nlat_half) × Complex{NF}`. At T31, nlat_half=24, nfreq_max ≈ 48, NF=Float32, nlayers=8, K=9: `48·72·24·8 = ~660K` complex × 8 bytes = ~5 MB per buffer × 2 buffers = 10 MB. Negligible. At T127 (nlat_half ~96, nfreq_max ~192, nlayers=32) and K=9: ~50× larger, ~500 MB. Still fine on GPU but worth budgeting; can quantify exactly during implementation.

### 6.10 `parameterization_tendencies_only!` — does the same refactor break it?

[tendencies.jl:519-555](SpeedyWeather/src/dynamics/tendencies.jl#L519-L555) is used when `dynamics=false`. It transforms u, v, temperature, humidity grid tendencies (parameterization output) into spectral and applies curl/div for vor/div tendencies.

With u, v, T, q grid tendencies now living in views of the fused parent: this path can do its own batched transform of just those 3–4 named slots (or, simpler, just keep the existing per-variable transforms — the fused parent buffer happens to be sized for the dycore's needs but the views can still be transformed individually).

No correctness issue — just an opportunity to also batch-transform here. Lower priority.

### 6.11 Order/numerical reproducibility

The reorder doesn't change which arithmetic operations happen — only the order of `transform!` calls. Within `transform!` itself the math is deterministic per layer. So bit-reproducibility is preserved between the old per-variable transforms and the fused transform, **assuming** `_legendre!` and `_fourier!` use the same per-layer reductions whether `nlayers = K·model_nlayers` or `K` separate calls of `model_nlayers`. A spot check would confirm this — the inner Legendre loop is over `lm` and the matvec is per-layer-independent ([legendre.jl:62-94](SpeedyTransforms/src/legendre.jl#L62-L94)), so it should be bitwise identical.

### 6.12 Failure modes / things that could surprise us

1. **Field views in KA kernels on Metal** — historically flaky for non-trivial views. Mitigation: contiguous slice (`:, k1:k2`) only. Smoke-test early.
2. **LTA views** — may not exist; might need to be added. Check first.
3. **Adapt rules** — `Adapt.@adapt_structure` recurses into fields. The fused parent + views need to adapt cleanly when moving Variables to GPU. Likely fine since `field_view` returns a regular `Field` with a `SubArray` data — `Adapt` treats it like any AbstractArray. Verify.
4. **Differentiability / Enzyme** — the existing `_reset_tendencies!` was carefully written to avoid Union-typed iteration ([tendencies.jl:1013](SpeedyWeather/src/dynamics/tendencies.jl#L1013)). Our changes mostly preserve type stability if we avoid runtime branching on `haskey(...)` for humidity. Plan: dispatch on model type (`PrimitiveDry` vs `PrimitiveWet`) for the layout choice, so the slot indices are compile-time constants per method.
5. **Tests under `--check-bounds=yes`** — `field_view` uses `@boundscheck` semantics consistent with regular Field; should be a no-op concern.
6. **CPU regression** — fused transform has a much larger inner working set (`g_north`, `g_south`, column scratch). Cache-blocked CPU Legendre may slow down at K≥9. Plan: have a CPU dispatch fallback that sends the parent through `transform!` slot-by-slot if K·nlayers exceeds some threshold. Decide threshold by benchmarking.

### 6.13 Step-by-step implementation order for §6

Each step below is a self-contained PR.

S1. **LTA-view helper** (only if LTA views don't already work) — add `lta_view(::LowerTriangularArray, :, k1:k2)` analogous to `field_view`. Add a tiny test.

S2. **Transform-batch kwarg** — add `transform_batch_layers` (or similar) to `SpectralGrid` / `SpectralTransform`. Default 1. Verify memory growth.

S3. **Allocate the fused parent + rewire named tendencies as views** — implement `fuse_tendencies!(vars, model)` (Approach (i) from §6.3). Run after `Variables(model)`. Verify that all existing tests pass with the rewire applied but no transform-batching yet (transforms stay per-variable). This catches issues with views as drop-ins early.

S4. **Wire the auxiliary slots into `dynamics_tendencies!`** — replace `vars.scratch.grid.a/b` reuse with the fused parent slots. `dynamics_tendencies!` declares its `uT_grid`, `vT_grid`, etc. as `field_view`s of the parent. Transforms still per-variable — this is just preparing the data layout. Verify tests.

S5. **Replace the inner transforms with one mega-batch** — pack the 9 (or 6) variables into the parent (already done in S4), call `transform!` once on the parent, run the spectral post-ops on the named spectral-side views. This is the actual perf win step. Benchmark.

S6. **Pack the 2D pressure tendency into the parent** — extend parent by 1 layer, drop the standalone `transform!(ūv̄∇lnpₛ, ...)` call. Benchmark.

S7. **Tracer integration** — extend the parent layout to include tracer slots when active tracers exist. Or: if the tracer count is low and dynamic, keep tracer transforms separate. Decide based on workload.

### 6.14 Open questions — answered

**(Q1) LTA views.** **Already exists.** `lta_view` is defined at [`LowerTriangularArrays/src/lower_triangular_array.jl:726-728`](LowerTriangularArrays/src/lower_triangular_array.jl#L726-L728), as the LTA analogue of `field_view`:

```julia
lta_view(L::LowerTriangularArray, c::Colon, i, args...) =
    LowerTriangularArray(view(L.data, c, i, args...), L.spectrum)
```

It returns a real `LowerTriangularArray` whose `.data` is a `SubArray`. It is already used in production by [`SpeedyWeather/src/variables/variables.jl:307`](SpeedyWeather/src/variables/variables.jl#L307) for leapfrog-step views (`get_step` in the prognostic state). So the spectral operators (`curl!`, `divergence!`, `∇²!`) demonstrably already accept LTA views as drop-ins. No new helper needed.

**(Q2) Where `SpectralTransform` picks up `nlayers`.** The convenience constructor at [`SpeedyWeather/src/dynamics/spectral_grid.jl:264-273`](SpeedyWeather/src/dynamics/spectral_grid.jl#L264-L273) pulls `nlayers` directly from the `SpectralGrid`:

```julia
function (::Type{S})(spectral_grid::SpectralGrid; one_more_degree=true, kwargs...) where {S<:AbstractSpectralTransform}
    (; NF, spectrum, grid, nlayers) = spectral_grid
    ...
    return S(spectrum, grid; NF, nlayers, kwargs...)
end
```

The inner `SpectralTransform(spectrum, grid; NF, nlayers, ...)` constructor at [`SpeedyTransforms/src/spectral_transform.jl:90`](SpeedyTransforms/src/spectral_transform.jl#L90) accepts `nlayers` as a keyword. The model itself constructs the transform via `SpectralTransform(spectral_grid)` — see [`primitive_dry.jl:102`](SpeedyWeather/src/models/primitive_dry.jl#L102), [`primitive_wet.jl:107`](SpeedyWeather/src/models/primitive_wet.jl#L107), etc.

**Implication for the plan.** Two clean ways to give the transform a larger scratch:
- (a) Add a `transform_batch::Int = 1` field to `SpectralGrid`. The convenience constructor multiplies: `S(spectrum, grid; NF, nlayers = nlayers * transform_batch, kwargs...)`. Default 1 keeps back-compat. Primitive-equation models set it ≥ K.
- (b) Bypass the convenience: have the primitive-equation model build its transform via the inner constructor with an explicit `nlayers` override.

I recommend (a) — it is one extra integer in `SpectralGrid` and zero ambiguity at call sites.

**Note on the kwargs collision**: the convenience constructor destructures `nlayers` from `spectral_grid` and passes it explicitly, so passing `nlayers` again through `kwargs...` would error. Either (a) or (b) above sidesteps this; option (a) does it inside the constructor.

**(Q3) `Adapt.adapt_structure` on a Field with `SubArray` data.** The relevant rules are:

```julia
# RingGrids/src/field.jl:524
Adapt.adapt_structure(to, field::AbstractField) = Adapt.adapt(to, field.data)

# LowerTriangularArrays/src/lower_triangular_array.jl:793
Adapt.adapt_structure(to, L::LowerTriangularArray) = Adapt.adapt(to, L.data)
```

Note: both rules unwrap to the underlying data array — they do **not** rebuild a Field/LTA on the target side. That's intentional: KA kernels receive the raw array view, which works because the kernels we have all index by `[ij, k]` style which is just `AbstractArray` indexing. Existing production code already does this with views (e.g. [`drag.jl:96-100`](SpeedyWeather/src/dynamics/drag.jl#L96-L100), [`particle_advection.jl:117-118`](SpeedyWeather/src/dynamics/particle_advection.jl#L117-L118), output writers in [`output/variables/`](SpeedyWeather/src/output/variables/)) — passing `field_view(...)` results into kernels works today on CPU and CUDA.

For `on_architecture` (which actually moves the data across devices, not just adapt for a kernel call), the path is different — defined at [`field.jl:529-...`](RingGrids/src/field.jl#L529). Our fused parent buffer is allocated once on the target architecture and the views are constructed locally on that architecture, so we never need `on_architecture` to adapt a `SubArray`-backed Field. The pattern: allocate parent with `on_architecture(arch, zeros(NF, grid, K*nlayers))`, then `field_view` on it.

**Verdict:** safe in principle. Smoke-test on each GPU backend (CUDA at least) is still warranted but the existing `field_view`/`lta_view` usage means we're not breaking new ground.

**(Q4) Mass-conservation gate `pres_tend.data[1:1] .= 0`.** Trivially compatible. If `pres_tend = lta_view(parent_spec, :, slot_p)`, then:
- `pres_tend.data` is a 2D `SubArray` (the slice at layer `slot_p` of `parent_spec.data`).
- `pres_tend.data[1:1]` is a length-1 view of that SubArray's first row.
- `.= 0` writes a single scalar zero into the first harmonic at that layer.

The implementation stays one line. (As a nano-optimization, we could replace the `[1:1] .= 0` with a scalar `[1] = 0` to avoid the broadcast wrapper — that's item 2.11 from earlier.)

---

### 6.15 Net effect on Phase B / S0

All four open questions resolved positively:
- Q1: `lta_view` ready, in use elsewhere.
- Q2: clear path via a `transform_batch::Int` kwarg on `SpectralGrid`.
- Q3: `field_view`/`lta_view` are already exercised on GPU code paths in production; smoke-test rather than design risk.
- Q4: trivial.

So **S0 collapses to one preliminary**: add `transform_batch` to `SpectralGrid` and thread it through the convenience constructor. Then proceed straight to S1, which becomes effectively a no-op (no new helpers needed) and we can start with S3 (`fuse_tendencies!`).

**Revised S-sequence:**
- ~~S1 (lta_view helper)~~ — not needed, exists.
- **S1':** Add `transform_batch::Int = 1` to `SpectralGrid`; multiply into `nlayers` in the SpectralTransform convenience constructor. Default unchanged → no behavior change. Add a test that an `S` constructed with `transform_batch=K` accepts an LTA/Field of `K·model_nlayers` layers.
- **S2 (✅ DONE — see §6.16):** Declarative fusion via a `fuse::Symbol` field on variable definitions. Allocator builds a shared parent and swaps named entries to views. Cross-type fusion supported. Production tendencies `vars.tendencies.grid.{u,v,T,q}` are now views of `vars.tendencies.grid.tend_grid`.
- **S3:** Re-wire `dynamics_tendencies!` auxiliary slots (`uT/vT/uq/vq/KE`) to use parent views instead of `scratch.grid.{a,b}` reuse — by declaring them as `ScratchVariable(..., fuse=:tend_grid)` so they share the same parent automatically.
- **S4:** Replace per-variable transforms with one `transform!(parent_spec, parent_grid, ...)` and run spectral post-ops on slot views. **This is the perf win step.** Benchmark.
- **S5:** Pack 2D pressure tendency into the parent (one extra layer). Drop the standalone 2D `transform!`.
- **S6:** Tracer integration (extend layout with active tracer slots, or keep tracers on the side).

---

### 6.16 Implemented: declarative fuse mechanism (S2 ✅)

A general fusion mechanism is now in place. **Any** variable definition can opt in via a `fuse::Symbol` keyword; the allocator handles the rest.

#### How it works

Add `fuse=:foo` to a variable definition. All variables sharing the same `(namespace, fuse)` pair — across any combination of variable types (`PrognosticVariable`, `GridVariable`, `TendencyVariable`, `DynamicsVariable`, `ParameterizationVariable`, `ParticleVariable`, `ScratchVariable`) — are allocated as views of a single shared parent buffer.

The parent itself lives at one canonical location: `vars.scratch.fused.<fuse_symbol>`. Views appear at their natural locations (`vars.tendencies.grid.u`, etc.) and reach the same memory.

Example:
```julia
# barotropic.jl
TendencyVariable(:u, Grid3D(), namespace=:grid, fuse=:tend_grid, ...)
TendencyVariable(:v, Grid3D(), namespace=:grid, fuse=:tend_grid, ...)
# primitive_dry.jl
TendencyVariable(:temperature, Grid3D(), namespace=:grid, fuse=:tend_grid, ...)
# primitive_wet.jl
TendencyVariable(:humidity, Grid3D(), namespace=:grid, fuse=:tend_grid, ...)
```

After `Variables(model)` for `PrimitiveWetModel` with `nlayers=4`:
```
vars.scratch.fused.tend_grid      # parent Field, size (npoints × 16) — canonical home
vars.tendencies.grid.u            # field_view of parent, slots 1–4
vars.tendencies.grid.v            # field_view of parent, slots 5–8
vars.tendencies.grid.temperature  # field_view of parent, slots 9–12
vars.tendencies.grid.humidity     # field_view of parent, slots 13–16
```

Add a `ScratchVariable(:foo, Grid3D, namespace=:grid, fuse=:tend_grid)` and the parent grows to 20 layers; the scratch view appears at `vars.scratch.grid.foo` (slots 17–20). The parent itself stays at `vars.scratch.fused.tend_grid` only — no duplicate references in other groups.

#### Key implementation details

- **Where**:
  - [`SpeedyWeather/src/variables/variables.jl`](SpeedyWeather/src/variables/variables.jl) — added `fuse::Symbol = Symbol()` to the `@kwdef`-generated struct (uniform across all 7 variable types). `build_fuse_parents` runs once before per-type allocation, walking all variables across all types and allocating one parent per `(namespace, fuse)`. `_allocate_namespace` emits views for fused members; standalone members allocate normally via `zero(v, model)`. After per-type allocation, `build_fused_namespace` collapses the fuse-parents Dict into a `(; fuse_sym = parent, ...)` NamedTuple and merges it into `scratch` as `scratch.fused`.
  - [`SpeedyWeather/src/variables/dimensions.jl`](SpeedyWeather/src/variables/dimensions.jl) — `fused_slots(dims, model)` and `allocate_fused(vars, model)` methods for `Grid2D`, `Grid3D`, `Spectral2D`, `Spectral3D`. Fallback errors clearly for unsupported dim types. All members in a fuse group must share the same dim type (validated at allocation).
  - [`SpeedyWeather/src/dynamics/tendencies.jl`](SpeedyWeather/src/dynamics/tendencies.jl) — `_reset_tendency!` is the original simple `fill!`. When a tendency entry is a fuse view (e.g. `vars.tendencies.grid.u`), `fill!` writes through the view and zeroes only that slot of the parent, leaving non-tendency slots in the same parent untouched.

- **Cross-type sharing**: when fuse-group members live in multiple variable types (e.g. one `TendencyVariable` and one `ScratchVariable` with the same fuse symbol), they all become views of the same canonical parent at `vars.scratch.fused.<fuse_symbol>`. The parent is referenced exactly once.

- **Reset semantics**: `reset_tendencies!` walks `vars.tendencies` and `fill!`s each leaf as before. Fuse views fill their slot of the parent directly. Pure-scratch fuse groups (no tendency members) are never touched by `reset_tendencies!` — their state persists across timesteps.

- **Fuse-symbol uniqueness**: enforced globally at construction. Two `(namespace, fuse)` pairs sharing the same fuse symbol error with a clear message. Fuse symbols must be unique across the whole model (since they're flat keys in `vars.scratch.fused`).

- **Default behavior**: `fuse = Symbol()` (empty) → no fusion, exact previous behavior. All existing variable declarations continue to work unchanged.

#### What was fused for the dycore use case

| Model | Fused tendency-grid variables | Parent size |
|-------|-------------------------------|-------------|
| Barotropic | `u`, `v` | `(npoints × 2·nlayers)` |
| PrimitiveDry | `u`, `v`, `temperature` | `(npoints × 3·nlayers)` |
| PrimitiveWet | `u`, `v`, `temperature`, `humidity` | `(npoints × 4·nlayers)` |

Left intentionally **standalone** (not fused):
- `vorticity` Grid3D `:grid` tendency — declared but not transformed by `dynamics_tendencies!` for primitive equations.
- `divergence` Grid3D `:grid` tendency — same.
- `pressure` Grid2D `:grid` tendency — different dim type; can't fuse with Grid3D under the same-dim-type constraint.

#### Verified

- `PrimitiveDryModel` and `PrimitiveWetModel` build with the expected fused parent at `vars.scratch.fused.tend_grid` and slot views at `vars.tendencies.grid.{u, v, T[, q]}`. Parent does NOT appear under fuse symbol in any other group.
- Round-trip writes: `vars.scratch.fused.tend_grid .= 99` then `reset_tendencies!` zeroes the tendency slots (1..3·nlayers / 1..4·nlayers) — i.e. fill! through views correctly writes through to the parent at the right offsets.
- Cross-group sharing: `ScratchVariable` with `fuse=:tend_grid` shows up at `vars.scratch.<namespace>.<member>` as a view of the same canonical parent.
- 1-hour `run!` completes for both PrimitiveDry and PrimitiveWet — parameterizations writing to view-backed `vars.tendencies.grid.{u, v, T, q}`, dycore reading them, per-variable `transform!` calls on view-backed Fields, and `reset_tendencies!` all work transparently.

#### Limitations / known follow-ups

- ~~**GPU smoke-test still pending**~~ — **VERIFIED on CUDA (§6.17)**. All mechanisms work correctly on GPU.
- **Cross-architecture transfer** of an existing `Variables` (i.e. `on_architecture(arch, vars)` moving from CPU to GPU) is not fuse-aware — would adapt parent and views independently and break sharing. Not exercised in the typical flow (Variables are constructed once on the target arch), so deferred unless tests trip it.
- **Mixed dim types in one fuse group** (e.g. Grid2D pressure + Grid3D tendencies) not supported. Could be added by sizing the parent as Grid3D with Grid2D members occupying 1 layer each.

---

## 7. Status

**Infrastructure**
- [x] **S2: Declarative fuse mechanism** (§6.16) — `fuse::Symbol` field on all variable types; cross-type sharing; reset-aware. `vars.tendencies.grid.{u,v,T,q}` are now views of the shared parent `vars.tendencies.grid.tend_grid`. Verified by 1h `run!` of PrimitiveDry + PrimitiveWet.
- [ ] **S0/S1': `transform_batch::Int` on `SpectralGrid`** — needed before S4 to size the SpectralTransform scratch for `K·nlayers`.
- [ ] **GPU smoke-test of view-backed tendencies** — confirm the existing fuse machinery survives a `run!` on CUDA/Metal/AMDGPU.

**Phase A — easy local wins** (do not depend on fusion)
- [ ] A1 (2.5): fuse 2D spec→grid for ∇lnpₛ
- [ ] A2 (2.7): drop temperature anomaly add/subtract
- [ ] A3 (2.10): fold linear pressure gradient into bernoulli
- [ ] A4 (2.8): reuse ūv̄∇lnp
- [ ] A5 (2.9): fold linear_virtual_temperature into spectral geopotential

**Phase B — multi-variable transform batching** (builds on S2)
- [ ] B1 (2.3): fuse u_tend/v_tend (already declared as views of `tend_grid`; need to actually batch the transform)
- [ ] B2 (2.1): fuse temperature tendency block — declare `uT, vT` as ScratchVariables with `fuse=:tend_grid`, then mega-transform their slots together with `temp_tend`
- [ ] B3 (2.2): fuse humidity tendency block — same pattern with `uq, vq, humid_tend`
- [ ] B4 (2.4): mega-batch grid→spec — single `transform!(parent_spec, parent_grid, ...)` call
- [ ] S5: pack 2D pressure tendency into the parent
- [ ] S6: tracer integration

**Phase C — kernel fusion**
- [ ] C1 (2.6): kernel fusion of grid-space prep kernels
- [ ] C2 (2.13): LTA broadcast cleanup
- [ ] C3 (2.11/2.12/2.14): small remaining items, profile-driven
