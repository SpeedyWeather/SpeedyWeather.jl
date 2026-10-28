# Output variable group and model diagnostics

> Status: **in progress**. Part 1 (`OutputVariable` and `variables.output`) is implemented; part 2 (`model.diagnostics`) is a proposal awaiting sign-off.

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

### Part 2 (proposal, not implemented)

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

Open questions for sign-off: whether `output` should be written to restart files
(for continuous accumulations across restarts); whether existing parameterization accumulations
(`rain_convection += …`) should move into a diagnostic that accumulates the instantaneous rate
instead, so the parameterizations only write rates.

## Testing and verification

- `test/variables/variables.jl`: new testset checking the moved variables live in `vars.output`
  and not in `vars.parameterizations`, `OutputVariable` type.
- Existing convection, large-scale condensation, land time stepping and output tests updated to
  the new paths and run.

## Documentation changes

`docs/src/variables.md` and `docs/src/variable_system.md` list the new `output` group and
`OutputVariable`; `docs/src/examples_3D.md` uses `simulation.variables.output`.

## Known limitations

Breaking for user code accessing `simulation.variables.parameterizations.rain_convection` etc.

## Future work

Part 2 above after sign-off.
