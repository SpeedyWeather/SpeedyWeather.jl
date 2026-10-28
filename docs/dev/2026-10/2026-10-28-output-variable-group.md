# Output variable group and model diagnostics

> Status: **in progress**. `OutputVariable`, `variables.output` and `model.diagnostics` with `Accumulate` and `PeriodicMaximum` are implemented; "rate since last output" and monthly means are future work.

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

Who resets: the diagnostic owns its schedule (e.g. a `PeriodicSchedule`, or "since last output"
bound to `model.output.interval` at `initialize!`), and the output writer only reads. A model has
a single `output` writer, so accumulators are never shared between writers. Callbacks may also
read `variables.output`, but they only read and never reset.

With this, the post-hoc `AbstractRateOutputVariable` diff in NetCDF `finalize!` can be replaced by
an on-the-fly "accumulated since last output / elapsed time" diagnostic that works for every writer.

Implemented (`output/diagnostics.jl`):

- `AbstractDiagnostic`, `diagnostic!(vars, diagnostic, model)` and `initialize!(vars, diagnostic, model)`
  fallbacks; `diagnostics!(vars, model)` called in `time_step!(simulation)` after the clock is
  stepped (so spin-up steps are included, as before) and before `output!`.
- `diagnostics::DG` field in `PrimitiveWetModel` (default `PrecipitationDiagnostics()`) and
  `PrimitiveDryModel` (default `(;)`). `initialize!(model.diagnostics, model)` in the model's
  `initialize!`, `initialize!(vars, model.diagnostics, model)` before the initial conditions so
  `StartFromFile` can overwrite them.
- `Accumulate(rate, name; namespace, units, desc)`: adds `Δt * rate` every time step, `Δt` being
  `model.time_stepping.Δt` that the parameterizations used to compute their rates, so results are
  unchanged. Also declares the rate so that it works (as zero) when the component is `nothing`.
- `PrecipitationDiagnostics()`: `Accumulate`s for large-scale and convective rain and snow and
  `land.river_runoff`. Convection and large-scale condensation no longer accumulate; the soil
  moisture filter writes a new `land.river_runoff_rate` [m/s] instead of accumulating.
- `PeriodicMaximum(; variable, group, namespace, name, units, schedule)`: running maximum
  `<name>_running` (initialized to -Inf), copied to `<name>` (initialized to NaN) and reset when
  `schedule` (`PeriodicSchedule(every = Day(1))`) is scheduled, branchless for Reactant.
- Restart files store `variables.output`, `StartFromFile` restores the entries present in both.

Not yet: "rate since last output" (a diagnostic resetting on the output schedule) and monthly means.

## Testing and verification

- `test/variables/variables.jl`: new testset checking the moved variables live in `vars.output`
  and not in `vars.parameterizations`, `OutputVariable` type.
- Existing convection, large-scale condensation, land time stepping and output tests updated to
  the new paths and run; convection and soil moisture tests now check rates × Δt.
- `test/output/diagnostics.jl`: `Accumulate` adds `Δt * rate` and works without convection,
  `PeriodicMaximum` running maximum and end-of-day reset (unit and in a 30h run), restart round trip
  of `variables.output`.

## Documentation changes

`docs/src/variables.md` and `docs/src/variable_system.md` list the new `output` group and
`OutputVariable`; `docs/src/examples_3D.md` uses `simulation.variables.output`. New page
`docs/src/diagnostics.md` on default, added and custom diagnostics.

## Known limitations

Breaking for user code accessing `simulation.variables.parameterizations.rain_convection` etc.
Passing `diagnostics` replaces the defaults, so precipitation needs `PrecipitationDiagnostics()...`
merged in, otherwise e.g. `ConvectiveRain` output fails. `PeriodicMaximum` is 2D only.
Barotropic and shallow water models have no `diagnostics` yet.

## Future work

- "Rate since last output" diagnostic bound to `model.output.interval`, replacing the post-hoc
  `AbstractRateOutputVariable` diff in the NetCDF `finalize!`.
- Time means over calendar periods (monthly mean temperature), needs a calendar-aware schedule.
- `diagnostics` for barotropic and shallow water models.
