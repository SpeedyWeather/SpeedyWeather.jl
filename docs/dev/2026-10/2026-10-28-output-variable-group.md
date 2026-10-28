# Output variable group and model diagnostics

> Status: **in progress**. `OutputVariable`, `variables.output` and `model.diagnostics` (`Accumulate`, `Maximum`, `Mean`, reset by the output writer) are implemented, addressing review comments.

Date of initial draft: 2026-10-28

Base revision: e292baf012d70813f9eab23081661bf1e168c3fb

## Originating prompt

> Some precipitation variables are accumulated and intended for output only. Can you open a PR
> with a new variable group "output" that's intended for output purposes only and move those
> variables in there by defining a new `OutputVariable` like `ParameterizationVariable` this
> shouldn't clash with the OutputVariables module. Then I am wondering whether it's possible to
> define a new model.diagnostics field (maybe a named tuple) that decides on which variables to
> add to variables.output. I want to be able to declare easily accumulated variables like
> precipitation vs rate since last output or monthly mean temperature that can be accumulated on
> the fly and which the output writer can be responsible for. Because output.variables is write
> only for the dynamics/parameterizations it's safe to hand this responsibility to the output writer

## Revision log

- 2026-10-28: Split into part 1 (implemented in this PR) and part 2 (design only, needs sign-off).
- 2026-10-28: Option A chosen (the diagnostic owns its reset schedule). Option B dropped: a model has only one output writer.
- 2026-10-28: Signed off: (1) `variables.output` is written to and restored from restart files,
  (2) parameterizations only write rates, accumulation moves into `Accumulate` diagnostics.
  Requested in addition: a `PeriodicMaximum` diagnostic to illustrate a non-accumulating diagnostic.
  Part 2 implemented in the same PR.

> Yes to 1 and 2. Can you also implement something like a DailyMaximum diagnostic just to
> illustrate how a non-accumulating diagnostic would be set up?
- 2026-10-28: `DailyMaximum` renamed to `PeriodicMaximum` as its period is set by its `schedule`;
  fallbacks of `diagnostic!`/`initialize!` untyped in all but the diagnostic to avoid method
  ambiguities with custom diagnostics (docs build failure).

> The CI for docs are failing, have a look? Also given that DailyMaximum takes a period as option,
> what about naming it PeriodicMaximum?
- 2026-10-28: Review: reverses option A. Diagnostics only accumulate, the output writer resets
  them after every output so that the output interval sets the period of the statistic.
  `PeriodicMaximum` (two variables + schedule) becomes `Maximum` (one variable), new `Mean`;
  the nested rate output variables (`ConvectiveRainRate` etc.) and their post-processing are
  removed, precipitation output is the accumulation since the last output. Diagnostics no
  longer declare the variables they read (warn on initialization, skip silently). `Δt` of the
  accumulation is the model time of the step (`elapsed_time`, Δt/2 on leapfrog spin-up steps),
  also for the river runoff rate. Compile-time unrolled loop over the diagnostics.

> I believe the diagnostic should be responsible for the accumulation (with +, or max etc) but the
> output writer for the reset, so that the output interval determines the statistic. E.g. How would
> you output weekly averages if not the output interval is weekly? [PR review, abridged]
- 2026-10-28: Review: `Maximum` resets to -Inf instead of NaN, plain `max` so NaNs propagate.
  Clarified the two `initialize!` methods of diagnostics (before/after variable allocation).

> isn't the `initialize!` called before the variables are allocated?
> Shouldn't the reset be to -Inf? Otherwise you'll have `max(NaN, x) = NaN`?

## Problem description

Accumulated precipitation (`rain_convection`, `snow_convection`, `rain_large_scale`,
`snow_large_scale`) and river runoff (`land.river_runoff`) are declared as
`ParameterizationVariable`s although no parameterization ever reads them: they are only
incremented (`+=`) and read by the output. Mixing them with the parameterization working arrays
hides that they carry memory across time steps and are safe to process or reset by the output.

Rates are currently derived from accumulations only post hoc in `finalize!` of the NetCDF output
(`diff` along time of the accumulated variable), which only works for NetCDF, needs the whole
record in memory, and there is no general way to request other on-the-fly statistics (time means,
maxima, accumulations since last output).

## Background

`Variables` has hardcoded groups (`prognostic`, `grid`, `tendencies`, `dynamics`,
`parameterizations`, `particles`, `scratch`, `fused`). Each group has a variable declaration type
generated in a loop in `variables/variables.jl`, e.g. `ParameterizationVariable`; `Variables(model)`
filters all declared variables by that type and allocates them into the group.
`SpeedyWeather.OutputVariables` (in `output/variables/OutputVariables.jl`) is a module of output
*definitions* (`AbstractOutputVariable`s such as `ConvectiveRain`), the name `OutputVariable`
(singular, a variable declaration) does not clash with it.

## Summary of changes

### Part 1 (this PR)

- New variable declaration type `OutputVariable` (generated alongside the others, exported).
- New group `Variables.output` (between `particles` and `scratch`), added to
  `ALL_VARIABLE_GROUPS`, so `copy!`, `on_architecture`, JLD2/Array output group filtering,
  Enzyme/FiniteDifferences extensions pick it up automatically.
- Moved to `OutputVariable`: `rain_convection`, `snow_convection`, `rain_large_scale`,
  `snow_large_scale`, `land.river_runoff`. Kernels, output paths, tests and docs updated
  (`vars.parameterizations.rain_convection` → `vars.output.rain_convection`).
- Contract (documented): the output group is write-only for dynamics and parameterizations;
  nothing in the model reads it, so output writers/diagnostics may process or reset it without
  changing the simulation.

The instantaneous rates (`rain_rate`, `snow_rate`, `*_rate_convection`, …) stay in
`parameterizations` because the land/snow/cloud schemes read them.

### Part 2 (implemented, as proposed below)

A `diagnostics` field in the primitive equation models (and later all models), a `NamedTuple` of
`AbstractDiagnostic <: AbstractModelComponent`s, e.g.

```julia
model = PrimitiveWetModel(spectral_grid;
    diagnostics = (
        precipitation = AccumulatedSinceOutput(:rain_large_scale, :rain_convection),
        temperature_monthly = TimeMean(:temperature, schedule = MonthlySchedule()),
    ))
```

Each diagnostic

1. declares its storage via `variables(d)` returning `OutputVariable`s (namespace
   `:diagnostics` or the diagnostic's key), so `Variables(model)` allocates them with no extra
   machinery — this is how `model.diagnostics` decides what goes into `variables.output`;
2. implements `accumulate!(vars, d, model)`, called once per time step after the time step and
   filters (so leapfrog/multi-step schemes accumulate the new state once), as a broadcast/kernel
   so it is GPU-compatible;
3. implements `reset!(vars, d)` and `finalize_period!(vars, d)` (e.g. sum → mean = sum/count).

Who resets (revised in review): the diagnostic accumulates, the output writer resets all
diagnostics after writing (`reset_diagnostics!`, called by NetCDF/Zarr/HEALPix, JLD2 and Array
output). The output interval therefore sets the period: precipitation since the last output,
daily maximum with daily output, weekly mean with weekly output. Without output no reset.

Implemented (`output/diagnostics.jl`):

- `AbstractDiagnostic` with `variables(diagnostic)`, `initialize!(vars, diagnostic, model)`,
  `diagnostic!(vars, diagnostic, model)` and `reset_diagnostic!(vars, diagnostic, model)`;
  fallbacks untyped in all but the diagnostic to avoid method ambiguities with custom diagnostics.
  `diagnostics!`/`reset_diagnostics!` unrolled at compile time over the `NamedTuple`.
- `diagnostics!(vars, model)` is called in `time_step!(simulation)` after the time step and before
  the clock is stepped (spin-up steps included), so `elapsed_time(time_stepping, clock)` is the
  model time of this step: Δt/2, Δt/2, then Δt for `Leapfrog` (its clock rewinds by Δt/2 on the
  second step), Δt otherwise. `time_step` would be the leapfrog integration step 2Δt.
- `diagnostics::DG` field in `PrimitiveWetModel` (default `PrecipitationDiagnostics()`) and
  `PrimitiveDryModel` (default `(;)`). `initialize!(vars, model.diagnostics, model)` runs before
  the initial conditions so `StartFromFile` can overwrite them.
- `Accumulate(; rate, name, namespace, units, desc)`: adds `elapsed_time * rate`, reset to 0.
- `Maximum(; variable, group, namespace, name, units)`: running maximum, reset to -Inf (NaNs propagate).
- `Mean(; variable, group, namespace, name, units)`: incremental running mean with a step count
  (`<name>_count`), reset to NaN/0, readable at any time.
- Diagnostics do not declare the variables they read: warn on initialization, skip silently.
- `PrecipitationDiagnostics()`: `Accumulate`s for large-scale and convective rain and snow and
  `land.river_runoff`. Convection and large-scale condensation no longer accumulate; the soil
  moisture filter writes `land.river_runoff_rate` [m/s] = runoff / `elapsed_time`.
- Output: `ConvectiveRain`, `LargeScaleRain`, … are the accumulation since the last output [mm];
  `ConvectiveRainRate` etc., `AbstractRateOutputVariable` and the `finalize!` post-processing removed.
- Restart files store `variables.output`, `StartFromFile` restores the entries present in both.

## Testing and verification

- `test/variables/variables.jl`: new testset checking the moved variables live in `vars.output`
  and not in `vars.parameterizations`, `OutputVariable` type.
- Existing convection, large-scale condensation, land time stepping and output tests updated to
  the new paths and run; convection and soil moisture tests now check rates × Δt.
- `test/output/diagnostics.jl`: leapfrog `elapsed_time`; `Accumulate` adds `Δt * rate`, resets,
  warns and skips without convection; `ArrayOutput` every 6h resets all diagnostics and the sum
  over all outputs equals the accumulation of a run without output; `Maximum`/`Mean` unit tests;
  restart round trip of `variables.output`.

## Documentation changes

`docs/src/variables.md` and `docs/src/variable_system.md` list the new `output` group and
`OutputVariable`; `docs/src/examples_3D.md` uses `simulation.variables.output`. New page
`docs/src/diagnostics.md` on default, added and custom diagnostics.

## Known limitations

Breaking for user code accessing `simulation.variables.parameterizations.rain_convection` etc.
Passing `diagnostics` replaces the defaults, so precipitation needs `PrecipitationDiagnostics()...`
merged in, otherwise e.g. `ConvectiveRain` output fails. `Maximum` and `Mean` are 2D only.
The precipitation rate output variables are removed (breaking), precipitation output is per output
interval instead. Calendar periods (monthly means) require output at calendar intervals.
Barotropic and shallow water models have no `diagnostics` yet.

## Future work

- Monthly means need monthly (calendar) output intervals.
- 3D diagnostics.
- `diagnostics` for barotropic and shallow water models.
