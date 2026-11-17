# Replace ProgressMeter.jl by a lightweight progress meter in SpeedyWeatherInternals

> Status: **planned**. Plan drafted, awaiting sign-off. Drop the ProgressMeter.jl dependency (and its
> `[sources]` entries pointing at the unreleased ProgressMeter.jl#369 branch) by implementing the subset
> SpeedyWeather uses in a new submodule `SpeedyWeatherInternals.ProgressMeters`. Continues
> [progress-elements.md](../2026-10/progress-elements.md) in the same PR (#1298).

Date of initial draft: 2026-11-17

Base revision: eb3a7080472dc558e82af58a39eb83f3523373f4 (head of `mc/progress-elements`, 70 commits behind `main` at 19a6ba18)

## Originating prompt

> Can you have a look at PR1298, there's no progress upstream in ProgressMeter.jl and so I am
> considering to break with that dependency and implement a lightweight version with only the features
> we need in a small module in SpeedyWeatherInternals. Have a look at the subset of functionality we
> use in SpeedyWeather and create a plan. All of this should go into PR1298.

## Revision log

## Problem description

PR #1298 builds the progress line from `ProgressMeter.Elements`, which only exists in the open
ProgressMeter.jl#369. That PR has not moved for weeks, so #1298 needs `[sources]` entries pointing at a
fork branch in 9 environments, fails on Julia 1.10 (no `[sources]`), and cannot be merged or released.
At the same time SpeedyWeather uses only a small part of ProgressMeter.jl (one `Progress` type, `next!`,
`finish!`, a handful of elements), and the part it does use already had to be worked around: the
`speedstring` piracy (removed in #1298), an Enzyme `make_zero` rule for the `ReentrantLock` inside
`ProgressCore`, and `remaining_time`, a copy of ProgressMeter's ETA logic for `ProgressTxt`.

## Background

What SpeedyWeather uses of ProgressMeter.jl today (on this branch):

| Where | What |
|---|---|
| `SpeedyWeather/src/output/feedback.jl` | `Progress(n; enabled, elements, desc, color, barlen, barglyphs, dt)`, `BarGlyphs(" ━━  ")`, `next!`, `finish!`, `durationstring`; fields `counter`, `n`, `start`, `tinit`, `tlast` |
| `SpeedyWeather/src/output/progress_elements.jl` | `Elements.AbstractProgressElement`, `print_element`, `Colored`, `Description`, `Percentage`, `Bar`, `ETA`, `Speed`, `ElapsedTime`, `Counter`; fields `tcurrent`, `tinit`, `counter`, `start` |
| `SpeedyWeather/ext/SpeedyWeatherEnzymeExt.jl` | `make_zero(::Type{<:ProgressCore}, ...)` returning the primal (the lock can't be zeroed) |
| `SpeedyTransforms/src/matrix_transform.jl` | `Progress(n; dt = 2, desc)` and `next!` while precomputing the transform matrices |
| tests | `test/output/feedback.jl` refers to `ProgressMeter.speedstring` and `ProgressMeter.Elements.Bar` |

Not used: `ProgressThresh`, `ProgressUnknown`, `@showprogress`, `progress_map`/`progress_pmap`,
`update!` with options, `showvalues`, `offset` (multiple stacked meters), `cancel`, `ijulia_behavior`/
`IJulia.clear_output`, the threading lock, the adaptive `check_iterations` predictor, the
`Threshold`/`Spinner` elements, the `legacy_width` bar.

How ProgressMeter draws (to be reproduced): every `dt` seconds the line is rendered by calling
`print_element(element, p)` for every element, bars last so they fill the width the others leave
(`displaysize(output)`), neighbouring elements of the same color merged, then written in one go as
`"\r" * line * "\e[K"` to `stderr`, colored only if the stream has `:color`. When the counter reaches
`n` (or `finish!`), the line is drawn a last time with `finished = true` (100%, full bar, ETA becomes
`Time: <elapsed>`) followed by a newline. Nothing is printed when `enabled = false`, and the final
line only if a line has been printed before.

## Summary of changes

### 0. Merge `main` into `mc/progress-elements`

The branch is 70 commits behind. Merge first (resolve CHANGELOG, Project.toml conflicts), so the new
work starts from current `main`.

### 1. New submodule `SpeedyWeatherInternals.ProgressMeters` (~250 lines)

`SpeedyWeatherInternals/src/ProgressMeters/ProgressMeters.jl`, included in `SpeedyWeatherInternals.jl`
like the other submodules. New stdlib dependency `Printf` for `SpeedyWeatherInternals`.

```julia
module ProgressMeters
export Progress, BarGlyphs, next!, finish!,
    AbstractProgressElement, print_element, Colored, Description, Percentage, Bar, ETA, Speed,
    ElapsedTime, Counter, default_elements
```

- `BarGlyphs(leftend, fill, front, empty, rightend)`, all `Char`, and `BarGlyphs(s::AbstractString)`
  from 5 characters. Default `BarGlyphs("|██ |")`, a single front glyph; the sub-character
  front glyphs (`['▏','▎',...]`) of ProgressMeter are dropped (SpeedyWeather uses `" ━━  "`).
- `mutable struct Progress{O <: IO}` with fields
  `n, counter, desc, color, enabled, dt, output, barlen::Union{Int, Nothing}, barglyphs, elements::Tuple,
  tinit, tlast, tcurrent, printed, finished`. No `start` (always 0 in SpeedyWeather; `SimulationSpeed`
  and the ETA use `counter` directly), no lock, no `ProgressCore` split.
  Constructor `Progress(n; desc = "Progress: ", color = :green, enabled = true, dt = 0.1, output = stderr,
  barlen = nothing, barglyphs = BarGlyphs(...), elements = default_elements())`. `output::IO` is stored
  as given; with the default `stderr` the field is `IO`-typed (dynamic dispatch only when printing).
- `next!(p; step = 1)`: increments `counter`; returns immediately if `!enabled`; draws if
  `time() > tlast + dt`, and the final line once `counter >= n`. `time()` is called on every `next!`
  instead of ProgressMeter's adaptive `check_iterations`: a time step takes ms, `time()` ~20 ns.
  Not drawing must not allocate (tested), as `next!` is called every time step.
- `finish!(p)`: sets `counter = n` and draws the final line if anything was printed before.
- Elements, same interface as ProgressMeter.jl#369 so `ProgressElements` and user elements carry over:
  `abstract type AbstractProgressElement`, `print_element(element, p)::String` (strings print
  themselves), `Colored(element, color)`, `Description()`, `Percentage()`, `Bar()`, `ETA()`, `Speed()`,
  `ElapsedTime()`, `Counter()`, `default_elements() = (Description(), Percentage(), Bar(), ETA())`.
- Internal helpers: `render_line(p)` (segments, bars last), `draw!(p)` (single write of
  `"\r" * colored line * "\e[K"`, plus `"\n"` when finished, then `flush`), `barstring`,
  `durationstring`, `speedstring` (per-iteration, for `Speed`), `eta_string(p)` (shared by `ETA` and
  `ProgressTxt`, replacing `remaining_time`). Coloring by hand (ANSI codes from `Base.text_colors`
  if `get(io, :color, false)`), as ProgressMeter does, to write the line at once.
- Explicitly not implemented: threading lock, `showvalues`, `offset`, IJulia `clear_output`/stdio
  throttling reset, `ProgressThresh`/`ProgressUnknown`, macros.

### 2. SpeedyWeather

- `Project.toml`: remove `ProgressMeter` from `[deps]`, `[sources]`, `[compat]`.
- `SpeedyWeather.jl`: `import SpeedyWeatherInternals.ProgressMeters` instead of `import ProgressMeter`.
- `output/progress_elements.jl`: import the generic elements from `SpeedyWeatherInternals.ProgressMeters`
  instead of `ProgressMeter.Elements`; `SimulationSpeed` uses `p.counter` (no `p.start`). Everything
  else (bound elements, separators, `VerticalCourantNumber`) unchanged.
- `output/feedback.jl`: `ProgressMeters.Progress`, `ProgressMeters.BarGlyphs`; `ProgressTxt` uses
  `ProgressMeters.eta_string`, `remaining_time` removed. The progress line looks exactly as before.
- `ext/SpeedyWeatherEnzymeExt.jl`: the `make_zero` rule moves from `ProgressMeter.ProgressCore` to
  `ProgressMeters.Progress` (no lock any more, but the `IO` field still should not be zeroed/copied);
  drop it if the differentiability tests pass without it.

### 3. SpeedyTransforms

`matrix_transform.jl` uses `SpeedyWeatherInternals.ProgressMeters.Progress(n; dt = 2, desc = ...)`;
`ProgressMeter` removed from `SpeedyTransforms/Project.toml` (`[deps]`, `[compat]`).

### 4. Environments

Remove the `ProgressMeter` entries (`[deps]` and `[sources]`) from `SpeedyWeather/test`,
`test/GPU/{CUDA,AMDGPU,MetalGPU}`, `test/differentiability`, `test/reactant`, `benchmark`,
`benchmark/CUDA` and `docs`. With that the PR has no unreleased dependency and Julia 1.10 CI works again.

### 5. Versions, CHANGELOG, PR

- `SpeedyWeatherInternals`: `0.2.3` → `0.2.3+DEV` (additive, non-breaking). On release it becomes
  `0.2.4` and the `SpeedyWeatherInternals` compat of `SpeedyWeather` and `SpeedyTransforms` must be
  raised to `0.2.4`, otherwise `0.2.3` (without `ProgressMeters`) is allowed. (Alternative: `0.3.0-DEV`
  now with compat `"0.3"`, which also requires compat updates in RingGrids and LowerTriangularArrays.)
- `SpeedyTransforms` already `0.2.1+DEV`, `SpeedyWeather` already `0.23.0-DEV`: no bump.
- CHANGELOG: rewrite the existing #1298 line ("... removes the ProgressMeter.jl dependency, progress
  meter now in `SpeedyWeatherInternals.ProgressMeters`" instead of "Depends on ProgressMeter.jl#369").
- PR description: drop the draft/dependency note, describe the new module; mark ready for review.
- Old plan `progress-elements.md`: revision-log entry pointing here, status updated.

## Testing and verification

- `SpeedyWeatherInternals/test/progress_meters.jl` (new, in `runtests.jl`), all into an `IOBuffer`:
  rendering of each element and of a custom element; `Colored` only emits ANSI codes with
  `IOContext(io, :color => true)`; bar width from `barlen` and from `displaysize` (via
  `IOContext(io, :displaysize => (24, 80))`); `dt` throttling (`dt = Inf` prints nothing until finished,
  `dt = 0` prints every step); final line has `100%`, a full bar, `Time:` and a trailing newline;
  `enabled = false` prints nothing; `finish!` is a no-op output-wise if nothing was printed;
  `durationstring`/`barstring` edge cases; `@allocated next!(p) == 0` while not drawing.
- `SpeedyWeather/test/output/feedback.jl`: replace the `ProgressMeter.speedstring`/`Elements.Bar`
  references (`ProgressElements.Bar === SpeedyWeatherInternals.ProgressMeters.Bar`); add a run with
  `Feedback(verbose = true)` writing to an `IOBuffer` and check the final line matches the documented
  layout. Existing element tests stay.
- SpeedyTransforms matrix-transform tests unchanged, must still pass.
- Run with `--check-bounds=yes`: SpeedyWeatherInternals tests, SpeedyTransforms tests, SpeedyWeather
  `output/feedback.jl`, `output/callbacks.jl`, `dynamics/dispatch.jl` (JET), and the differentiability
  tests touching `Feedback` (Enzyme rule). Julia 1.12 locally; CI covers 1.10.
- Visual check: run a T31 `PrimitiveWetModel` in a terminal before/after and compare the line
  (screenshot/recording in the PR), incl. the ETA → `Time:` switch on finish.
- Runic formatting.

## Documentation changes

- `docs/src/feedback.md`: say the elements come from `SpeedyWeatherInternals.ProgressMeters`
  (re-exported by `ProgressElements`) instead of ProgressMeter.jl, and document the custom-element
  interface there instead of pointing to ProgressMeter.jl's docs.
- Docstrings for `Progress`, `next!`, `finish!` and every element.

## Known limitations

- No thread safety: `next!` must be called from one task (true for `Feedback` and the matrix transform).
- No IJulia-specific handling (ProgressMeter resets IJulia's stdio byte counter and can clear the cell);
  in Jupyter the line is still overwritten via `\r`, but very long runs may hit IJulia's output throttling.
- Non-terminal outputs (log files, CI) get the `\r`-separated lines like with ProgressMeter.
- Only a single front glyph in bars (no sub-character resolution).
- Loses upstream fixes in ProgressMeter.jl; the module is small enough to maintain.

## Future work

- If ProgressMeter.jl#369 is merged, the element interface is the same, so moving back is mechanical,
  but there is little reason to.
- A plain (no `\r`) mode for non-TTY outputs, e.g. one line every N percent, which could also replace
  `ProgressTxt`'s own percentage logic.
