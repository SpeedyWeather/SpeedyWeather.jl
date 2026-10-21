# Reactant reverse-mode (AD) benchmarks, and a timestep multiplier for benchmark201 only

> Status: **in progress**. Implemented and verified locally on Reactant CPU (Phase 0 confirmed budget
> 32 and the default physics). Outstanding: the full run on the H100 (`reactant-gpu --reverse`) and the
> PR number in the CHANGELOG line.

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

### 4. A run without `--reverse` never overwrites the reverse results

A run replaces its architecture's whole record. Without special handling, a run without `--reverse`
would therefore delete earlier reverse results. To prevent that, the reverse suite gets its own field
in the record:

```json
"reactant-gpu": {
    "meta": {...}, "markdown": "...", "overview": {...},
    "reverse": {"timestamp": "...", "markdown": "...", "data": {"truncation": [...], "nsteps": [...], "forward_time": [...], ...}}
}
```

The rules:
- **With `--reverse`**, the driver builds a new `reverse` field from this run, which replaces the old
  one.
- **Without `--reverse`**, the driver copies the `reverse` field unchanged from the stored record of
  the same architecture, if there is one. For this, `load_results()` moves before the point where the
  record is built.
- **Other architectures' records** are never touched, as today.
- **README.md and the docs page keep the section.**
  - `arch_markdown()` skips `:benchmark800`.
  - `markdown` is the regular suites followed by `reverse.markdown`, whether that is new or copied.
  - README.md and the docs page therefore keep showing the reverse section, without any change to
    `docs/generate_benchmarks_page.jl`.
- **A copied section shows its age.** Its heading carries its own measurement date
  (`reverse.timestamp`), so it is recognisably older than the rest of the architecture's section.
- **Unreadable store.**
  - Today, if `benchmark_results.json` exists but can't be parsed, `load_results()` warns and starts
    from an empty store. The run then overwrites the file, which loses every architecture's results,
    the reverse ones included.
  - Instead, the driver checks at startup, before any suite runs, that the store is readable, and
    errors if it is not. This fails within seconds, not after hours of benchmarking.
  - A missing file still starts a fresh store.
  - The merge at the end of the run re-reads the file as today.
- **Structured data.** `data` holds the numbers in structured form for later figures. Non-finite
  values are written as `null`, as in `overview`.

### 5. The `--hi-res` flag (added on request, see the revision log)

- `--hi-res` adds T342 and T512 at 8, 16 and 24 layers to `:benchmark201`, inserted into each layer
  block of the sweep.
- **GPU only.** On `gpu-nvidia`, `gpu-amd` and `reactant-gpu` (the new driver constant `ON_GPU`, which
  also sets the matrix transform's cap) the six rows are added. On CPU architectures the flag is ignored
  with a warning, since these resolutions are too slow there.
- The matrix transform rows are unchanged, as their cap is below these truncations.
- The new rows get 50 time steps (the floor of `n_timesteps`); the multiplier applies to them too.
- **Not kept.** Unlike the reverse-mode results, these rows are part of benchmark201's table and
  overview. A later run without `--hi-res` therefore drops them for that architecture.

### Files

| File | Change |
| --- | --- |
| `benchmark/manual_benchmarking.jl` | Flag parsing (`--reverse`, `--hi-res`); `ON_GPU`; conditional `using Enzyme` and include; multiplier for benchmark201 only; carrying the reverse record over; check at startup that the results store is readable; usage header; the README "Running the benchmarks" section lists the flags and the multiplier |
| `benchmark/define_benchmarks.jl` | `:benchmark800`, behind `RUN_REVERSE`; T342 and T512 in `:benchmark201`, behind `RUN_HI_RES` on GPU; 1-based truncations in comments |
| `benchmark/benchmark_suite.jl` | The comment on `timestep_multiplier`; 1-based truncation in the heading of the individual-dynamics-functions tables |
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
2. **Carry-over.** After step 1, run again without `--reverse`, restricted to one quick suite.
   - The `reverse` field in the JSON is identical to before.
   - README.md and `docs/generate_benchmarks_page.jl` still show the reverse section, with its original
     date.
   - Then make the JSON unreadable (e.g. truncate a copy and point the driver at it). The driver errors
     at startup and leaves the file unchanged.
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
- `--reverse-only`, to rerun only this suite. With the carry-over from section 4 in place it is about 15
  lines.
- Record the multiplier in benchmark201's README section.

## Decisions

Settled on 2026-10-21:

1. **Resolutions:** T32, T64, T128 and T256 at L8, for now.
2. **Checkpointing:** one `Binomial(budget)` for every configuration. Phase 0 sets the budget.
3. **Opt-in mechanism:** the `--reverse` flag. A run without it never overwrites the reverse results
   (section 4).
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
3. **Carry-over.** `reactant-cpu` without `--reverse` updates the record's timestamp, while the
   `reverse` field is byte-identical to before. README.md still shows the section with its original date.
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
