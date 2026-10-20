# Move output variables into an exported `OutputVariables` module

> Status: **planned**. Plan drafted, awaiting sign-off before implementation.

Date of initial draft: 2026-10-19

Base revision: e410eed05c967a4e117b91c91591d17da155a774

## Originating prompt

> Can you restructure the output variables? I would like to have a PR with the module
> OutputVariables that's always exported with all the TemperatureOutput etc in it? So that
> the access would be OutputVariables.Temperature() not SpeedyWeather.TemperatureOutput()?
> This would be breaking!

## Revision log

## Problem description

The ~40 output variable types (`VorticityOutput`, `TemperatureOutput`, …) and the
collection functions that return tuples of them (`DynamicsOutput()`, `LandOutput()`, …)
are defined straight in `SpeedyWeather`. None of them is exported, so users write
`SpeedyWeather.TemperatureOutput()`. That is long, the `Output` suffix repeats itself,
and output *variables* look the same as output *writers* (`NetCDFOutput`, `ZarrOutput`,
`HEALPixOutput`), which also end in `Output`.

## Background

- Definitions live in `SpeedyWeather/src/output/variables/*.jl`. `output_variables.jl`
  includes them and defines `AllOutputVariables()` plus the shared helpers (`get_indices`,
  `is2D`, `path_or_nothing`, …).
- Each variable is an `AbstractOutputVariable` with a `path(::XOutput, simulation)` method,
  and some also have `output!` / `transform` methods.
- The writers' defaults in `output/writers/general.jl` call these constructors, e.g.
  `default_output_variables`.
- The Terrarium extension adds `TerrariumOutput` (a function stub that is exported) and
  `TerrariumOutputVariable`.

## Summary of changes

1. New submodule `SpeedyWeather.OutputVariables`, defined in
   `src/output/variables/OutputVariables.jl`. It wraps the existing variable files and
   imports the parent's abstract types and generic functions to extend
   (`AbstractOutputVariable`, `path`, `output!`, `DocStringExtensions`, …).
   The helpers (`get_indices`, `is2D`, `path_or_nothing`, …) stay in `SpeedyWeather`,
   because the writers use them.
2. Drop the `Output` suffix inside the module:

   | Old | New |
   |---|---|
   | `SpeedyWeather.TemperatureOutput()` | `OutputVariables.Temperature()` |
   | `SpeedyWeather.VorticityOutput()` | `OutputVariables.Vorticity()` |
   | `SpeedyWeather.TracerOutput(:co2)` | `OutputVariables.Tracer(:co2)` |
   | `SpeedyWeather.DynamicsOutput()` (tuple) | `OutputVariables.Dynamics()` |
   | `SpeedyWeather.LandOutput()` (tuple) | `OutputVariables.Land()` |
   | `SpeedyWeather.AllOutputVariables()` | `OutputVariables.All()` |
   | … the same rule for all ~40 types and 9 collections | |

   Writers keep their names (`NetCDFOutput`, `ZarrOutput`, `HEALPixOutput`, `ArrayOutput`).
3. `SpeedyWeather` exports `OutputVariables`, the module name only. The module's own names
   are not exported, so they cannot clash with `Temperature`-like names elsewhere.
4. Terrarium: `TerrariumOutput(...)` becomes `OutputVariables.Terrarium(...)`. The stub
   moves into the module and the extension adds methods to it.
5. Breaking: the old names are removed, with no deprecation aliases (per the prompt).
   Version bump: `SpeedyWeather` minor, e.g. `0.x.y` → `0.(x+1).0-DEV`.
6. Every use in `src/`, `ext/`, `test/`, `docs/src/` and the benchmarks is updated.
   The CHANGELOG entry gets a migration note.

## Testing and verification

- The existing output tests (`test/output/*.jl`, plus the test files that use output
  variables) pass after the rename.
- One short new test checks that `OutputVariables` is exported, that
  `OutputVariables.All()` returns `AbstractOutputVariable`s, and that the old names are
  gone.
- Docs build: all `@ref`s and `@example` blocks resolve.

## Documentation changes

`docs/src/output.md`, `custom_netcdf_output.md`, `land.md`, `ocean.md` and `sea_ice.md` are
updated to the new names, plus a short paragraph on the `OutputVariables` namespace and on
how to define a custom variable: `struct MyVar <: SpeedyWeather.AbstractOutputVariable`
with `SpeedyWeather.path(::MyVar, simulation) = …`.

## Known limitations

- Breaking for every user script that names output variables.
- A short name like `OutputVariables.Temperature` no longer says "output" by itself; the
  module qualifier does that job.

## Future work

- Possibly move the writers into a similar `OutputWriters` namespace.
