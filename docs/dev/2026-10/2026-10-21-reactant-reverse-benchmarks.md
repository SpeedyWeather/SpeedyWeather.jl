# Reactant reverse-mode (AD) benchmarks, and a timestep multiplier for benchmark201 only

> Status: **in progress**. Implemented and verified locally (Phase 0 confirmed budget 32 and the
> default physics); `--hi-res`, `--SPEEDY` and `reactant-tpu` added on request, `--SPEEDY` not with
> Reactant yet, `reactant-tpu` untested. Outstanding: the full runs on the H100 and a TPU, and the PR
> number in the CHANGELOG line.

Date of initial draft: 2026-10-21

Base revision: f5ecd475691a5d2df15c0304589588a163765734 (branch `mg/benchmark-reverse`)

## Originating prompt

> - plan a new benchmark set for the main benchmark that is supposed to benchmark the performance of the
>   reverse model / AD. It shold time the time in seconds it takes for 50, 100, 200, 2000 time steps of
>   the model forward and for a gradient / reverse over that amount of steps with Reactant. This should
>   only benchmark the `SpectralTransform` variant at different resolutions. It should be opt-in and not
>   always executed (as it might take too long).
> - Also make the when given a multiplier in the command line to extend the simulation time that's
>   benchmark, that multiplier should only be applied to the main benchmark benchmark201

## Revision log

> 1) Four truncations is fine for now
> 2) yes go for Binomial
> 3) yes that's good, but can we achieve that it's not overwritten once it's launched without --reverse
> 4) yes the default model
>
> Just adjust the plan, don't implement yet

All four open decisions were settled as recommended (see "Decisions" at the end). Section 4 now spells
out the guarantee that a run without `--reverse` never overwrites the reverse results. It also adds a
check at startup that the results store is readable: today an unreadable store is silently replaced,
which would wipe the reverse results along with everything else. The testing step for the carry-over
is now more specific.

> okay go ahead

Phase 0 ran (results in "Protocol" below) and confirmed every assumption, so the implementation went
ahead as planned with budget 32. The benchmark environment's manifest was stale (Reactant 0.2.292 and
ProgressMeter 1.11, below the compat bounds), so it was regenerated with Julia 1.12.7, giving
Reactant 0.2.296. Enzyme was added to `benchmark/Project.toml`.

> add a further --hi-res flag that adds T342 and T512 to the :benchmark201
>
> edit benchmark scripts to use a 1-based notation throughout so T256 instead of T255 as we do in the
> model code
>
> Don't run CPU setups in the hi-res configurations

Added the `--hi-res` flag (section 5). It applies on GPU architectures only and adds T342 and T512 at
8, 16 and 24 layers. The benchmark scripts and this document now use 1-based truncations throughout:
the earlier T31, T63, T127 and T255 are now T32, T64, T128 and T256. The only computed label that
changed is the heading of the individual-dynamics-functions tables, which subtracted 1 from the
truncation. Results already stored in `benchmark_results.json` and README.md keep the old heading
until each architecture is rerun.

> okay I reconsidered, the reverse mode results should not be kept when the next benchmark is
> launched without --reverse
>
> add another flag --SPEEDY that signals a configuraiton that is more similar to original SPEEDY.f90.
> So a configuration with seasonal prescribed climatoglies in ocean and land, also quickly review if
> what other options might be best for that

The carry-over of section 4 was removed: every run now replaces everything stored for its label, as
before this PR, and the reverse suite is a section like any other. The check at startup that the
results store is readable stays. Added the `--SPEEDY` flag (section 6) with prescribed seasonal
climatologies of sea surface temperature, land surface temperature and soil moisture; section 6 also
reviews the other options.

> add a reactant-tpu option to run the benchmark on a tpu with reactant

Added the architecture `reactant-tpu` (section 7, 2026-10-22). It is treated like a GPU for the
resolution choices, and its machine info lists the TPU devices. It is untested on a TPU: none is
available here.

## Problem description

This plan makes two changes to `SpeedyWeather/benchmark/manual_benchmarking.jl`:

1. **A new opt-in suite that measures the cost of a gradient with Reactant.** It times, in seconds,
   50, 100, 200 and 2000 time steps forward, and the reverse-mode gradient over the same steps. It
   uses the `PrimitiveWetModel` with the `SpectralTransform` at several resolutions. The suite is too
   slow to run every time, so it runs only when requested.
2. **The timestep multiplier applies to `:benchmark201` only.** Today the optional second argument
   lengthens the timed runs of every `BenchmarkSuite`. It should lengthen only the main benchmark,
   `:benchmark201`, the PrimitiveWet resolution sweep.

## Background

**How the benchmarks run.**
- [manual_benchmarking.jl](../../../SpeedyWeather/benchmark/manual_benchmarking.jl) takes
  `arch [multiplier]` and loads the backend packages.
- It then includes three files:
  - `benchmark_suite.jl`: the suite types and the CPU/GPU hooks.
  - `reactant_runner.jl`: the Reactant hooks, only for `reactant-*`.
  - `define_benchmarks.jl`: the `benchmarks` dictionary.
- It runs every suite and writes one record per architecture into `assets/benchmark_results.json`.
  The record has three fields:
  - `meta`
  - `markdown`: all suite tables as one blob.
  - `overview`: the structured `:benchmark201` numbers.
- `README.md` and the docs page (`docs/generate_benchmarks_page.jl`) are regenerated from that JSON.
  The docs page prints `record["markdown"]` verbatim.
- A run replaces its architecture's whole record.

**The multiplier.**
- It is parsed at [manual_benchmarking.jl:86](../../../SpeedyWeather/benchmark/manual_benchmarking.jl#L86).
- It is set on every suite that has a `timestep_multiplier` field
  ([manual_benchmarking.jl:102](../../../SpeedyWeather/benchmark/manual_benchmarking.jl#L102)).
  That is every suite except `:benchmark700`.
- `BenchmarkSuite` applies it in `n_timesteps`
  ([benchmark_suite.jl:91](../../../SpeedyWeather/benchmark/benchmark_suite.jl#L91)).

**Reactant AD today.**
- `time_steps!(simulation, n; checkpointing, mincut)`
  ([SpeedyWeatherReactantExt.jl:115](../../../SpeedyWeather/ext/SpeedyWeatherReactantExt.jl#L115))
  runs n time steps as one traced loop, a `stablehlo.while` for n > 4.
- Reverse mode is tested in
  [test/reactant/differentation.jl](../../../SpeedyWeather/test/reactant/differentation.jl):
  - The loss is `sum(abs2, vorticity)` after n steps.
  - The gradient is `Enzyme.autodiff(Reverse, loss_n_steps, Active, Duplicated(vars, dvars),
    Const(model), Const(n), Const(options))`.
  - It is compiled with `reactant_compile_options(; raise_first = true)`.

What those tests showed with Reactant 0.2.294/0.2.295:
- With `checkpointing = Reactant.Binomial(2)`, the PrimitiveWet gradient at n = 6 compiles in about
  3.5 minutes in CI.
- `Reactant.Periodic(2)` does not finish compiling within an hour.
- `mincut = false` does not finish compiling within 30 minutes.
- `checkpointing = false` stores every step, so memory grows linearly with n.
- CI's PrimitiveWet gradient uses `convection = nothing` and a constant longwave transmissivity
  ([setup.jl](../../../SpeedyWeather/test/reactant/setup.jl)). That setup was chosen for the
  CPU-vs-Reactant comparison. As a result, nobody has yet compiled the gradient of the **default**
  PrimitiveWet model, with Betts-Miller convection.

**Recomputation with `Binomial`.** `Binomial(budget)` is revolve with `budget` checkpoints.
- n steps need r recomputation sweeps, where `binomial(budget + r, r) ≥ n`.
- Each step is recomputed up to r times.
- With budget = 32:
  - n = 50, 100 and 200 need r = 2, since `binomial(34, 2) = 561`.
  - n = 2000 needs r = 3, since `binomial(35, 3) = 6545`.

**Expected cost on reactant-gpu (H100).** The numbers come from `README.md`, for the PrimitiveWet
model at L8 with the `SpectralTransform`.

| Truncation | Forward time per step | 2000 steps forward |
| --- | --- | --- |
| T32 | 0.23 ms | 0.5 s |
| T64 | 0.36 ms | 0.7 s |
| T128 | 0.73 ms | 1.5 s |
| T256 | 2.3 ms | 4.5 s |

Compilation dominates the cost:
- Each forward program takes 16–45 s to compile (README).
- Each gradient program takes minutes.
- Every (resolution, n) pair is a separate program, because n is baked into it.

## Summary of changes

### 1. The multiplier applies to `:benchmark201` only

- In the driver loop
  ([manual_benchmarking.jl:98-106](../../../SpeedyWeather/benchmark/manual_benchmarking.jl#L98-L106)),
  iterate `for (key, suite) in benchmarks`.
- Set `suite.timestep_multiplier = TIMESTEP_MULTIPLIER` only when `key == :benchmark201`. Every other
  suite keeps the default of `1`.
- `n_timesteps` itself is unchanged.
- Update three texts:
  - the usage header of the driver;
  - its `@info` line, which should say that the multiplier applies to benchmark201 only;
  - the comment on the `timestep_multiplier` field
    ([benchmark_suite.jl:76](../../../SpeedyWeather/benchmark/benchmark_suite.jl#L76)).
- The new reverse suite has fixed step counts and no multiplier.

### 2. Opt-in with the `--reverse` flag

```bash
julia --project=. manual_benchmarking.jl reactant-gpu --reverse      # all suites + the reverse-mode suite
julia --project=. manual_benchmarking.jl reactant-gpu 10 --reverse   # and benchmark201 runs 10× longer
```

- Arguments that start with `--` are flags. All other arguments stay positional (`arch`, then
  `multiplier`), so existing invocations still work.
- An unknown flag is an error.
- `--reverse` with an architecture other than `reactant-*` fails before anything runs, with the
  message "the reverse-mode benchmarks need reactant-cpu or reactant-gpu".
- With `--reverse`, the driver also:
  - runs `using Enzyme` together with the other conditional backend loads. This loads
    `SpeedyWeatherEnzymeExt` and its `make_zero(::Variables)`.
  - includes the new `reverse_benchmark_suite.jl` after `reactant_runner.jl`.
- `Enzyme` is added to `benchmark/Project.toml`.

### 3. New suite `BenchmarkSuiteReverse`, as `:benchmark800`

The new file `benchmark/reverse_benchmark_suite.jl` defines:
- `BenchmarkSuiteReverse <: AbstractBenchmarkSuite`;
- its `run_benchmark!` method;
- its `write_results` method.

The existing `run_benchmark_suite!` therefore works unchanged: it catches errors per configuration and
calls `free_memory!` between configurations. `define_benchmarks.jl` adds the suite only if
`RUN_REVERSE` is set.

**Configurations.**
- Model: `PrimitiveWetModel` with default physics, L8, Float32 and the default grid.
- Transform: the `SpectralTransform`, built through
  `build_model(...; transform_kind = :default)`.
  - On Reactant, `build_model` passes the `SpectralTransform` explicitly.
  - Its Fourier stage is a truncated DFT.
- Truncations: T32, T64, T128 and T256 (1-based, as in the model code).
- Steps: 50, 100, 200 and 2000.
- That gives 16 configurations.
- On `reactant-cpu`, only truncations ≤ 128 run. This mirrors how benchmark201 caps the matrix
  transform per architecture.
- The suite stores these fields per configuration: truncation, nlayers, nlat, nsteps, checkpointing,
  Δt, forward time, gradient time, forward compile time, gradient compile time and status.

**Steps for one configuration.**
1. **Build and initialize.**
   - Build the model, then run `simulation = initialize!(model)` and
     `initialize!(simulation; steps = nsteps)`. This scales the variables and sets the clock.
   - Keep `vars0 = simulation.variables` as a template. It is never run.
2. **Compile the forward pass.**
   - Compile `loss_n_steps(vars, model, nsteps)` with `@compile`, using `REACTANT_AD_COMPILE_OPTIONS`.
   - These are the gradient's compile options, so Gradient/Forward measures the differentiation and
     not a difference in compile options.
   - Record the compile time.
3. **Compile the gradient.**
   - Compile `Enzyme.autodiff(Reverse, loss_n_steps, Active, Duplicated(vars, dvars), Const(model),
     Const(nsteps), Const((; checkpointing)))` with `@compile`.
   - Record the compile time.
4. **Time each program.**
   - Run one warm-up call, then up to 3 timed samples. Report the minimum.
   - If one sample takes longer than 60 s, stop after it. This keeps 2000-step runs at high
     resolution affordable.
   - Untimed preparation before each call:
     - `vars = deepcopy(vars0)`.
     - For the gradient, also a fresh `dvars = Enzyme.make_zero(vars0)`, because `autodiff`
       accumulates into the shadow.
   - The timed region is the call plus a wait for the device to finish:
     - Forward: read the returned loss as a `Float64`.
     - Gradient: `Reactant.synchronize` on the shadow of the vorticity's fused parent array, an output
       of the program.
5. **Sanity check (untimed).** After the last sample, check that the vorticity gradient is finite and
   not all zero. Record the status as `ok` or `non-finite gradient`.

**Checkpointing.**
- All configurations use one scheme: `Reactant.Binomial(budget)`.
- Phase 0 sets the budget; the provisional value is 32.
- The scheme is stored per configuration, so it appears in the table and can later be varied per
  configuration.
- Why not the alternatives:
  - Storing every step (`false`) needs memory proportional to n. 2000 copies of the state at T256 do
    not fit.
  - `Periodic` does not compile (see Background).

**Loss.**
- The loss is the same as in the tests: `sum(abs2, vars.prognostic.vorticity.data)` after the n steps.
- The gradient is taken with respect to all of `vars` (`Duplicated`). The model is `Const`.
- `loss_n_steps` and `REACTANT_AD_COMPILE_OPTIONS` are short, so they are copied into the benchmark
  file. Including them from the test file is not an option: that file runs tests when included and
  needs `Test`.

**Table.** `write_results` writes one row per configuration:

| T | L | Rings | Δt | Steps | Checkpointing | Forward | Gradient | Gradient/Forward | Compile forward | Compile gradient |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |

- Times are in seconds, with 3 significant digits.
- If part of a configuration fails, its cell shows the status, as benchmark201 does for SYPD.
- The section heading includes the measurement date, because this section can be older than the rest
  of the architecture's section (see 4).

### 4. Storing the results (revised: no carry-over)

A run replaces its architecture's whole record, as before this PR. A run without `--reverse`
therefore drops earlier reverse results. The first implementation copied them over; that carry-over
was removed on request (see the revision log).
- **The reverse suite is a section like any other.** Its table is part of the record's `markdown`,
  after the other suites, so README.md and the docs page show it without any change to
  `docs/generate_benchmarks_page.jl`.
- **Structured numbers.** With `--reverse`, the record also gets a `reverse` field with the suite's
  numbers (`reverse_data`), for later figures:

  ```json
  "reactant-gpu": {
      "meta": {...}, "markdown": "...", "overview": {...},
      "reverse": {"truncation": [...], "nsteps": [...], "forward_time": [...], "gradient_time": [...], ...}
  }
  ```
- **Other architectures' records** are never touched, as before.
- **Unreadable store.**
  - Today, if `benchmark_results.json` exists but can't be parsed, `load_results()` warns and starts
    from an empty store. The run then overwrites the file, which loses every architecture's results,
    the reverse ones included.
  - Instead, the driver checks at startup, before any suite runs, that the store is readable, and
    errors if it is not. This fails within seconds, not after hours of benchmarking.
  - A missing file still starts a fresh store.
  - The merge at the end of the run re-reads the file as before.
- Non-finite values are written as `null`, as in `overview`.

### 5. The `--hi-res` flag (added on request, see the revision log)

- `--hi-res` adds T342 and T512 at 8, 16 and 24 layers to `:benchmark201`, inserted into each layer
  block of the sweep.
- **GPU only.** On `gpu-nvidia`, `gpu-amd` and `reactant-gpu` (the new driver constant `ON_GPU`, which
  also sets the matrix transform's cap) the six rows are added. On CPU architectures the flag is ignored
  with a warning, since these resolutions are too slow there.
- The matrix transform rows are unchanged, as their cap is below these truncations.
- The new rows get 50 time steps (the floor of `n_timesteps`); the multiplier applies to them too.
- A later run without `--hi-res` drops them for that architecture (section 4).

### 6. The `--SPEEDY` flag (added on request, see the revision log)

- **What changes.** `build_model` gives every primitive equation model of every suite (including the
  dynamics-functions and reverse suites) SPEEDY.f90-like components, `speedy_components` in
  `benchmark_suite.jl`. Components given per run (`model_kwargs`) take precedence. Barotropic and
  ShallowWater models are unchanged.
  - Ocean: `SeasonalOceanClimatology` (monthly SST climatology, linear in time) instead of `SlabOcean`.
  - Land surface temperature: `SeasonalLandTemperature` instead of `LandBucketTemperature`.
  - Soil moisture: `SeasonalSoilMoisture` instead of `LandBucketMoisture` (PrimitiveWet only).
- **Stored separately.** The results go under the label `<arch>-speedy`, with
  `meta.configuration = "SPEEDY"` and a note at the top of the section. Default and SPEEDY results
  therefore coexist; the overview table gets its own column and the docs page its own line. Flags are
  case-insensitive.
- **Not with Reactant yet.** `--SPEEDY` with `reactant-*` errors at startup, see "SPEEDY configuration
  with Reactant" below.

**Review: SPEEDY.f90 against SpeedyWeather's options.**

| Aspect | SPEEDY.f90 | SpeedyWeather | In `--SPEEDY` |
| --- | --- | --- | --- |
| Ocean | Monthly SST climatology (optional slab anomaly model) | `SlabOcean` (default), `SeasonalOceanClimatology` | Yes |
| Land surface temperature | Climatology | `LandBucketTemperature` (default), `SeasonalLandTemperature` | Yes |
| Soil moisture | Monthly climatology | `LandBucketMoisture` (default), `SeasonalSoilMoisture` | Yes |
| Sea ice | Concentration climatology | `ThermodynamicSeaIce` (default), driven here by the prescribed SST; no sea-ice climatology (`PrescribedSeaIce` is a placeholder for coupling) | Default kept |
| Snow | Depth climatology | `SnowModel` only | Default kept |
| Vegetation, albedo | Climatologies, snow and sea-ice albedo | `VegetationClimatology`, `OceanLandAlbedo` | Already the default |
| Surface fluxes | Bulk formulae | Implemented after the SPEEDY documentation | Already the default |
| Shortwave | Two bands, diagnostic clouds | `OneBandShortwave` with SPEEDY's `DiagnosticClouds` | Already the default |
| Longwave | Four bands | One band, Frierson or constant transmissivity | No equivalent |
| Convection | Simplified mass flux | `BettsMillerConvection` only | No equivalent |
| Condensation | Relaxation above a humidity threshold | `ImplicitCondensation`, similar | Already the default |
| Vertical diffusion | SPEEDY's diffusion and shallow convection | `BulkRichardsonDiffusion` | No equivalent |
| Grid | Full Gaussian, 96×48 at T31 (SPEEDY's T30) | `OctahedralGaussianGrid` (default), `FullGaussianGrid` | Not included, option |
| Vertical levels | 8 σ levels, half levels 0, 0.05, 0.14, 0.26, 0.42, 0.6, 0.77, 0.9, 1 | `SigmaCoordinates(spectral_grid, σ_half)` | Not included, option |
| Time step | 40 min at T31 | 40 min at T32 (default) | Already the default |

The two options left are numerics rather than physics. Both are cheap to add but change the cost:
the full Gaussian grid has about 45 % more grid points than the octahedral one at the same number of
rings. SPEEDY's σ levels only apply to the L8 configurations, and a grid setting must not override the
Grids suite (`:benchmark400`).

**SPEEDY configuration with Reactant.** Two blockers in the package source, found with
`spike/speedy_config.jl` and `spike/speedy_fix_prototype.jl` (scratchpad):
1. **`Dates.month` of the traced clock branches** on traced values ("non-boolean
   `TracedRNumber{Bool}` used in boolean context"), and so do `firstdayofmonth` and `daysinmonth`. All
   three seasonal components compute month and weight that way. A shared branch-free helper fixes it:
   - `Dates.yearmonthday`, which Reactant provides branch-free, gives year, month and day;
   - days-in-month comes from an arithmetic formula, the leap year from `&`/`|`.

   Prototyped, it gives the same month, next month and weight as the current code for 2507 times
   covering 2000 (leap year) and 2001.
2. **The interpolation kernel indexes the monthly array with a traced month.**
   `interpolate_monthly_climatology_kernel!` takes `this_month` as a kernel argument, which Reactant
   passes as a `CuTracedRNumber`. It has no `to_index` inside the kernel ("invalid LLVM IR"). A fix would
   select the two months outside the kernel (a dynamic slice with Reactant, a view elsewhere) and pass
   only the slices and the weight. Not prototyped yet.

Both are changes to `SpeedyWeather/src/parameterizations` and need tests in `test/reactant`. They are
left for a decision.

### 7. The `reactant-tpu` architecture (added on request, see the revision log)

- **Usage.** `manual_benchmarking.jl reactant-tpu` runs on a TPU host, e.g. a Cloud TPU VM; Reactant
  downloads libtpu itself. It loads CUDA.jl and Reactant and calls
  `Reactant.set_default_backend("tpu")`. CUDA.jl loads without a GPU; Reactant's kernel raising needs
  it, and on TPU Reactant raises every kernel, as SpeedyWeather's compile options already do
  everywhere.
- **Like a GPU.** The new driver constant `ON_ACCELERATOR` replaces `ON_GPU` and covers
  `gpu-nvidia`, `gpu-amd`, `reactant-gpu` and `reactant-tpu`. On a TPU:
  - `--hi-res` adds T342 and T512;
  - the matrix transform goes up to T171;
  - the reverse suite goes up to T256.
- **Like the other Reactant architectures.** `--reverse` works; `--SPEEDY` is refused at startup.
- **Driver cleanup.** A single `ON_REACTANT` constant now replaces the repeated
  `reactant-cpu || reactant-gpu` checks.
- **Machine info.** The README section lists the number and kind of Reactant devices
  (`Reactant.devices()`, `Reactant.XLA.device_kind`), for every Reactant architecture. The benchmarks
  use one device.
- **Ordering.** `reactant-tpu` comes after `reactant-gpu` in README.md and on the docs page
  (`ARCH_ORDER` in both scripts).

What to expect on a TPU (not testable here):
- **Matrix multiply precision.** XLA's default precision for Float32 matrix multiplies on TPU is a
  single bfloat16 pass. On Reactant the Fourier transform (truncated DFT) and the Legendre transform
  are matrix multiplies, so the transforms run at reduced precision unless set otherwise. The
  benchmark uses the backend's default, as on the GPUs (TF32 on the H100), and doesn't set a precision.
  Reactant can raise it per compilation with `Reactant.with_config(; dot_general_precision =
  PrecisionConfig.HIGHEST)`.
- **Float64.** TPUs have no native Float64, so the Float64 configuration of `:benchmark300` will
  probably fail. Failures are recorded per configuration and the other configurations continue.
- **Memory.** The matrix transform at T171 needs about 17 GB (H100 numbers), more than the HBM of
  some TPU chips (e.g. 16 GB on v5e). It would then fail and be recorded the same way.

### Files

| File | Change |
| --- | --- |
| `benchmark/manual_benchmarking.jl` | Architecture `reactant-tpu`, `ON_REACTANT`, `ON_ACCELERATOR`, Reactant devices in the machine info; flag parsing (`--reverse`, `--hi-res`, `--SPEEDY`, case-insensitive); `RESULTS_LABEL` and `meta.configuration` for `--SPEEDY`; conditional `using Enzyme` and include; multiplier for benchmark201 only; the `reverse` numbers in the record; check at startup that the results store is readable; usage header; the README "Running the benchmarks" section lists the flags and the multiplier |
| `benchmark/define_benchmarks.jl` | `:benchmark800`, behind `RUN_REVERSE`; T342 and T512 in `:benchmark201`, behind `RUN_HI_RES` on GPU and TPU (`ON_ACCELERATOR`); 1-based truncations in comments |
| `docs/generate_benchmarks_page.jl` | `reactant-tpu` in `ARCH_ORDER` |
| `benchmark/benchmark_suite.jl` | `speedy_components`, `SPEEDY_CONFIGURATION` and their use in `build_model`; the comment on `timestep_multiplier`; 1-based truncation in the heading of the individual-dynamics-functions tables |
| `benchmark/reverse_benchmark_suite.jl` | New file |
| `benchmark/Project.toml` | Adds `Enzyme` |
| `CHANGELOG.md` | One line under `## Unreleased` |

No package source changes, so no version bumps.

### Phase 0: spike before implementing

A throwaway script in the session scratchpad runs on reactant-cpu locally, and on the H100 if it is
available. It uses the PrimitiveWet model at T32L8 with the `SpectralTransform` and default physics. It
answers three questions:

1. **Does the gradient compile with default physics?**
   - Test whether the gradient of `time_steps!` compiles for the default PrimitiveWet model, including
     convection.
   - If it does not, the fallback is the CI configuration: `convection = nothing` and a constant
     transmissivity, stated in the suite title. This needs the user's decision before going ahead.
2. **Does compile time grow with n or the budget?**
   - Measure the gradient's compile time and run time for n = 50 and n = 2000, with `Binomial(b)` for
     b ∈ {8, 32}.
   - If the revolve schedule is unrolled, compile time would grow with n, and n = 2000 would be
     infeasible.
3. **Does the budget fit in memory?** Check that `Binomial(32)` at 2000 steps fits: T128 on the CPU,
   T256 on the H100.

The spike's outcome is the budget, recorded in the revision log. If compile time grows strongly with n,
come back to the user before implementing.

## Testing and verification

The benchmark scripts have no unit tests, like the existing suites. To verify:

1. **Local run.** Run `manual_benchmarking.jl reactant-cpu --reverse` locally.
   - Temporarily restrict it to the reverse suite at T32 with n = 50 and 100.
   - Check that the table renders, the JSON has a `reverse` field and README.md has the section.
2. **Replacing.** After step 1, run again without `--reverse`, restricted to one quick suite. The
   reverse section and the `reverse` field are gone for that architecture; other architectures are
   unchanged. Then make the JSON unreadable (e.g. truncate a copy and point the driver at it). The
   driver errors at startup and leaves the file unchanged.
3. **Multiplier.**
   - `manual_benchmarking.jl cpu 2` makes the timed runs of benchmark201 twice as long. Check the
     logged step counts; the other suites are unchanged.
   - `cpu --reverse` errors immediately.
4. **Cross-check against benchmark201.** At n = 2000, the forward time per step should match
   benchmark201's SYPD for the same truncation within ±20 %.
5. **Full run on the H100.** Run `reactant-gpu --reverse` and commit the results together with
   README.md.
6. **Docs.** `docs/generate_benchmarks_page.jl` renders the new section.

## Documentation changes

- The driver's usage header and the README section "Running the benchmarks" (written by
  `write_preamble`):
  - document `--reverse`;
  - say that the multiplier applies only to benchmark201.
- The explanation list in the README preamble gets the new columns: Steps, Checkpointing, Forward,
  Gradient and Gradient/Forward.
- `CLAUDE.md`, "Benchmarks" section: one sentence on `--reverse`.

## Known limitations

- **Many compilations.** n is baked into each compiled program, so the suite needs 16 forward and 16
  gradient compilations. Phase 0 measured about 35–40 s per gradient compilation once Julia is warm,
  independent of n, the budget and the resolution. Including the runs, the estimate is 1–1.5 h on the
  H100 and about 1.5 h on reactant-cpu (12 configurations); on the CPU, T128 with 2000 steps alone takes
  about 40 minutes.
- **Different simulated periods.** The same step count covers different periods: 2000 steps are 55 days
  at T32 (Δt = 2400 s) but 7 days at T256 (Δt = 300 s).
- **Gradient/Forward includes recomputation.** The ratio includes revolve's recomputation, so with
  budget 32 it jumps at n = 2000 (r = 3 instead of 2). It reports the cost of this checkpointing scheme,
  not of the adjoint alone.
- **No memory numbers.** Reactant only exposes the process-wide allocator peak
  (`Reactant.XLA.allocatorstats`), not the peak of a single program.
- **The gradient can overflow.** Over 2000 steps the adjoint of the chaotic dynamics can grow large.
  A Float32 overflow is recorded as `non-finite gradient`; the timing is still valid.

## Future work

- Report the peak memory of each gradient program. This needs a wrapper for PJRT's
  `PJRT_Executable_GetCompiledMemoryStats` (`temp_size_in_bytes`), which Reactant only has as a raw C
  API struct. Memory is the other half of the checkpointing trade-off.
- Sweep checkpointing schemes: several budgets, and storing every step for small n.
- Make n a traced runtime value, so each resolution compiles once for all step counts. This works only if
  Enzyme-MLIR's binomial checkpointing supports a dynamic trip count. Reactant's error hint suggests
  that it does.
- Gradients with respect to parameters (`Duplicated(model)`).
- A docs figure from `reverse.data`.
- `--reverse-only`, to rerun only this suite. It would need to carry the other suites' stored results
  over.
- Record the multiplier in benchmark201's README section.

## Decisions

Settled on 2026-10-21:

1. **Resolutions:** T32, T64, T128 and T256 at L8, for now.
2. **Checkpointing:** one `Binomial(budget)` for every configuration. Phase 0 sets the budget.
3. **Opt-in mechanism:** the `--reverse` flag. Originally a run without it kept the stored reverse
   results; reconsidered later the same day: a run replaces everything stored for its label (section 4).
4. **Physics:** the default PrimitiveWet model, the same as benchmark201. If Phase 0 shows that its
   gradient doesn't compile, ask the user before switching to the CI configuration without convection.

## Protocol

### Phase 0 (2026-10-21, Apple M3, 16 GB, Reactant CPU backend)

Setup: Julia 1.12.7, Reactant 0.2.296, default `PrimitiveWetModel` (`BettsMillerConvection`), L8,
`SpectralTransform`, one Julia process per truncation, jobs in the order listed. The script is
`spike/spike.jl` in the session scratchpad (not committed). "Time" is the faster of the timed runs, the
first call excluded.

| T | Run | Compile | Time | per step |
| --- | --- | --- | --- | --- |
| T32 | forward, n = 50 | 161 s (cold) | 0.18 s | 3.6 ms |
| T32 | gradient, n = 50, Binomial(32) | 112 s (cold) | 1.39 s | 28 ms |
| T32 | gradient, n = 2000, Binomial(32) | 38 s | 65 s | 33 ms |
| T32 | gradient, n = 50, Binomial(8) | 35 s | 1.46 s | 29 ms |
| T32 | gradient, n = 2000, Binomial(8) | 38 s | 80 s | 40 ms |
| T32 | forward, n = 2000 | 7 s | 7.3 s | 3.6 ms |
| T128 | forward, n = 50 | 157 s (cold) | 4.1 s | 83 ms |
| T128 | gradient, n = 50, Binomial(32) | 106 s | 21 s | 420 ms |

The answers to the three questions:

1. **The default physics compiles.** The gradient of `time_steps!` for the default PrimitiveWet model,
   with convection, compiles and gives finite, non-zero gradients at every n and budget. No fallback
   needed.
2. **Compile time doesn't grow.** Gradient compile time is independent of n (50 or 2000), of the
   budget (8 or 32) and of the resolution: about 35–40 s once Julia is warm. "Cold" (the first program
   of its kind in a process) adds about 2 minutes of Julia compilation. The revolve schedule is not
   unrolled.
3. **Memory fits.** With 32 checkpoints the process peak was 2.6 GiB at T32 and 3.1 GiB at T128,
   against about 2 GiB without them. T256 on the H100 should need a few GiB for the checkpoints.

Further findings:
- **Budget 32 beats 8.** At 2000 steps the gradient takes 65 s with budget 32 against 80 s with
  budget 8 (19 % faster), as revolve predicts: 3 recomputation sweeps against 6. At 50 steps both take
  the same time. Budget 32 stays.
- **Gradient/Forward.** At T32 the ratio is 7.6 for 50 steps and 9.0 for 2000; at T128 it is 5.1 for
  50 steps.
- **The first call of a compiled program is slow.** In the first program of its kind per process it
  took 14–25 s longer than later calls. Later programs had less of this extra cost or none (up to
  10 s, e.g. the second gradient program). The suite's warm-up call absorbs it.
- **Synchronization works.** The compiled call returns after 0.002–0.25 s while the work takes up to
  85 s. After `Reactant.synchronize` on the fused parent of the vorticity shadow, copying the gradient
  to the host takes under 1 ms. So the timed region covers the whole computation.
- **Chaotic growth of the adjoint.** The maximum of the vorticity gradient grows from about 5e2 at
  50 steps to about 8e7 at 2000 steps (55 days at T32), still finite in Float32. The two budgets give
  the same gradient at 50 steps (to 4 digits) but differ by 9 % at 2000 steps. That fits rounding
  differences between the two programs' recomputations being amplified, the same chaotic amplification
  found before for the forward model. It doesn't affect the timings.

### Local verification of the implementation (2026-10-21, same machine)

The checks ran against a scratchpad copy of `SpeedyWeather/benchmark` with a reduced
`define_benchmarks.jl`: benchmark201 at T32L8 only, the reverse suite at T32 with n = 50 and 100. The
repository's suite definitions and results were not touched.

1. **Failure paths.** Each fails before any suite runs:
   - `cpu --foo` errors with "Unknown flag".
   - `cpu --reverse` errors with "need Reactant".
   - A truncated `benchmark_results.json` errors with "Could not parse" after 7 s and leaves the file
     byte-identical.
2. **`reactant-cpu --reverse`.**
   - The table renders after the regular suites, headed with its measurement date.
   - The JSON record of `reactant-cpu` has a `reverse` field (`timestamp`, `markdown`, `data`).
   - The record's `markdown` contains the section.
   - Both gradients are `ok`.
3. **Carry-over** (removed later on request, see the revision log). `reactant-cpu` without `--reverse`
   updated the record's timestamp, while the `reverse` field stayed byte-identical. README.md still
   showed the section with its original date.
4. **Multiplier.** `cpu 2` with a second `BenchmarkSuite` present: benchmark201 gets multiplier 2
   (2400 instead of 1200 steps at T32L8), the other suite keeps 1.
5. **Docs page.** A copy of `docs/generate_benchmarks_page.jl`, pointed at the test JSON, renders the
   reverse section; its heading carries the generator's ` — reactant-cpu` suffix.
6. **Formatting.** `runic --check` is clean for all four benchmark files.
7. **A transient slowdown.** The forward run of 100 steps in check 2 took 0.90 s (9 ms/step), against
   4.5 ms/step for 50 steps.
   - A follow-up script (`spike/forward_anomaly.jl`) timed n = 50, 100 and 200 with and without
     `GC.gc()` before each sample, and with and without first waiting for the copied input variables.
   - It measured 4.5–4.8 ms/step in every combination, so the harness is sound.
   - The 0.90 s was a slowdown of the laptop lasting a few seconds; it covered all three samples. The
     same machine also varies by about 30 % between processes (3.6 ms/step in Phase 0, 4.8 here).
     Expect steadier numbers on the HPC nodes.
   - The forward time per step agrees with benchmark201's SYPD in the same run (1364 SYPD = 4.8
     ms/step).

Not done locally: the full suite. A run on reactant-cpu would take about 1.5 h on this laptop and is
less informative than the H100 run, which is the next step.

### Verification of `--hi-res`, `--SPEEDY` and the removed carry-over (2026-10-21, same machine)

1. **`--hi-res`.** The suite definition, evaluated for `reactant-gpu`, `gpu-nvidia`, `reactant-cpu` and
   `cpu-arm`, gets the six rows T342/T512 × L8/16/24 on the GPU labels only; without the flag it is
   unchanged (27 configurations on GPU, 21 on CPU). One-off check that the model runs at these sizes
   (plain CPU, 3 steps each, before the GPU-only request): T342 1.3 s/step, T512 3.8 s/step, finite.
2. **`--SPEEDY` components.** PrimitiveDry and PrimitiveWet run 5 days on the plain CPU with all three
   climatologies: SST 240–303 K, soil temperature 228–310 K, soil moisture 0–0.47, finite.
3. **`--SPEEDY` in the driver** (scratchpad copy, reduced suites).
   - `cpu --SPEEDY` stores `cpu-arm-speedy` with `meta.configuration = "SPEEDY"` and the configuration
     note, and leaves `cpu-arm` byte-identical.
   - The overview table gets a `cpu-arm-speedy` column.
   - PrimitiveWet T32L8 needs 6.4 MB instead of 5.9 MB (the monthly fields), at a speed within this
     laptop's noise. ShallowWater rows are unchanged.
4. **Flags.** `reactant-cpu --speedy` and `reactant-gpu --SPEEDY --reverse` error at startup;
   `cpu --Speedy --foo` errors with "Unknown flag `--foo`" (so `--Speedy` was recognised).
5. **No carry-over.** `reactant-cpu --reverse` writes the `reverse` numbers and the section; the next
   `reactant-cpu` run without the flag drops both.
6. **Formatting.** `runic --check` is clean.

Not tested: `--SPEEDY` and `--hi-res` on a GPU (none here). The seasonal components use
KernelAbstractions kernels with plain integer months there, so they are expected to work, but nothing
tests them on GPU yet.
