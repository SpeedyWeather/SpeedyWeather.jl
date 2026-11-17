# Replace ProgressMeter.jl by a lightweight progress meter in SpeedyWeatherInternals

> Status: **in progress**. Plan signed off, implementing. Drop the ProgressMeter.jl dependency (and its
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

- 2026-11-17, sign-off: "Agree with version, and do feel free to simplify ProgressMeters wrt to the
  ProgressElements. In that ProgressMeter.jl PR some design choices have been made for a non-breaking
  change, here we'd be the only user so we can break for the sake of a simpler interface or a simpler
  default bar." `SpeedyWeatherInternals` gets `0.2.3+DEV`. Elements print only their value and
  `Progress` joins them with its `separator` (moved from `Feedback`/`bind_elements`); the bar carries
  its length and glyphs (`Bar(; length, fill, empty)`, no `BarGlyphs`, no partial glyphs); `Description`,
  `Colored`, `Speed` and `start` dropped; `Feedback.progress_bar_length` removed (no bar by default,
  as before).

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

### 1. New submodule `SpeedyWeatherInternals.ProgressMeters` (~200 lines)

`SpeedyWeatherInternals/src/ProgressMeters/ProgressMeters.jl`, included in `SpeedyWeatherInternals.jl`
like the other submodules. New stdlib dependency `Printf` for `SpeedyWeatherInternals`.

SpeedyWeather is the only user, so the interface is simplified instead of copying ProgressMeter.jl#369,
whose choices (`Description`/`Colored` elements, `barlen`/`barglyphs` on the meter, leading spaces inside
elements, `start`) were made to stay non-breaking for ProgressMeter's existing layouts.

```julia
module ProgressMeters
export Progress, next!, finish!,
    AbstractProgressElement, print_element, Percentage, Bar, ETA, ElapsedTime, Counter, default_elements
```

- `mutable struct Progress` with fields `n, counter, description, elements::Tuple, separator, color,
  enabled, dt, output::IO, tinit, tlast, tcurrent, printed, finished`. Constructor
  `Progress(n; description = "", elements = default_elements(), separator = ", ", color = :blue,
  enabled = true, dt = 0.1, output = stderr)`. No `start`, no lock, no `ProgressCore` split.
- **Elements print only their value** (`"23%"`, `"ETA: 0:00:07"`, `"2000-01-03"`), no leading spaces
  or separators. `Progress` assembles the line: the `description` (if not empty), then the
  non-empty element texts joined by `separator`, except next to a `Bar`, which is separated by a single
  space. This moves the separator logic of #1298 (`bind_elements` putting `Feedback.separator` in front
  of `AbstractSimulationElement`s) from SpeedyWeather into the meter and makes it uniform:
  ```
  23% ━━━━━━━━━━          ETA: 0:00:07, 2000-01-03, 36.99 years/day, 42 m/s, [-64, 21] ˚C
  ```
  Strings in `elements` print themselves (and are separated like any element).
- `print_element(element, p)::String` with `p.counter`, `p.n`, `p.tinit`, `p.tcurrent` (time of this
  redraw) and `p.finished` (final redraw) available, as in #369.
- Generic elements: `Percentage()`, `ETA()` (`Time: <elapsed>` once finished), `ElapsedTime()`,
  `Counter()`, and `Bar(; length = nothing, fill = '━', empty = ' ')`. The bar carries its own length
  and glyphs (no `BarGlyphs`, no `barlen` on the meter, no partial front glyphs, no end glyphs);
  `length = nothing` fills the terminal width the other elements leave, several bars share it.
  `default_elements() = (Percentage(), Bar(), ETA())`.
- Dropped compared to #369: `Description` (now the `description` field), `Colored` (one `color` for the
  whole line), `Speed` (s/it, SpeedyWeather has `SimulationSpeed`).
- `next!(p)`: increments `counter`; returns immediately if `!enabled`; redraws if
  `time() > tlast + dt`, and draws the final line once `counter == n` if something was printed before
  (as ProgressMeter: runs shorter than `dt` print nothing). `time()` is called on every `next!` instead
  of ProgressMeter's adaptive `check_iterations` (a time step takes ms, `time()` ~20 ns). `next!` must
  not allocate when not redrawing (tested).
- `finish!(p)`: sets `counter = n` and draws the final line if something was printed before.
- Drawing: the line is written at once as `"\r" * line * "\e[K"` (plus `"\n"` when finished), in
  `color` only if `get(output, :color, false)`, then `flush`.
- Helpers: `durationstring(seconds)` (`h:mm:ss`, days for long runs), `eta(p)` (remaining seconds,
  used by `ETA` and by `ProgressTxt`, replacing `remaining_time`).
- Not implemented: threading lock, `showvalues`, `offset`, IJulia `clear_output`/stdio throttling,
  `ProgressThresh`/`ProgressUnknown`, macros.

### 2. SpeedyWeather

- `Project.toml`: remove `ProgressMeter` from `[deps]`, `[sources]`, `[compat]`.
- `SpeedyWeather.jl`: `import SpeedyWeatherInternals.ProgressMeters` instead of `import ProgressMeter`.
- `output/progress_elements.jl`: `ProgressElements` re-exports the generic elements of
  `ProgressMeters`; SpeedyWeather's elements unchanged except that `bind_elements` only binds and
  filters (no separators any more). `AbstractSimulationElement` stays for the compact `show` of bound
  elements. `default_elements() = (Percentage(), ETA(), SimulationTime(), SimulationSpeed(),
  MaximumWindSpeed(), TemperatureRange())`, i.e. no bar, as before (`progress_bar_length = 0`).
- `output/feedback.jl`: **breaking**, `Feedback.progress_bar_length` removed, add `ProgressElements.Bar()`
  (or `Bar(length = 20)`) to `elements` instead. `Feedback.separator` is passed to `Progress`.
  `ProgressTxt` uses `ProgressMeters.eta` and `time() - tinit` (the previous `tlast - tinit` is never
  updated with `verbose = false`, giving a wrong speed in `progress.txt`); `remaining_time` removed.
- `ext/SpeedyWeatherEnzymeExt.jl`: the `make_zero` rule moves from `ProgressMeter.ProgressCore` to
  `ProgressMeters.Progress` (the `IO` field must not be zeroed or copied); drop it if the
  differentiability tests pass without it.

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
  raised to `0.2.4`, otherwise `0.2.3` (without `ProgressMeters`) is allowed.
- `SpeedyTransforms` already `0.2.1+DEV`, `SpeedyWeather` already `0.23.0-DEV`: no bump.
- CHANGELOG: rewrite the existing #1298 line ("... removes the ProgressMeter.jl dependency, progress
  meter now in `SpeedyWeatherInternals.ProgressMeters`" instead of "Depends on ProgressMeter.jl#369").
- PR description: drop the draft/dependency note, describe the new module; mark ready for review.
- Old plan `progress-elements.md`: revision-log entry pointing here, status updated.

## Testing and verification

- `SpeedyWeatherInternals/test/progress_meters.jl` (new, in `runtests.jl`), all into an `IOBuffer`:
  rendering of each element and of a custom element; separators (none next to a bar, empty texts
  skipped); ANSI codes only with `IOContext(io, :color => true)`; bar width from `length` and from `displaysize` (via
  `IOContext(io, :displaysize => (24, 80))`); `dt` throttling (`dt = Inf` prints nothing until finished,
  `dt = 0` prints every step); final line has `100%`, a full bar, `Time:` and a trailing newline;
  `enabled = false` prints nothing; `finish!` is a no-op output-wise if nothing was printed;
  `durationstring` edge cases; `@allocated next!(p) == 0` while not drawing.
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
- Bars have no sub-character resolution and no end glyphs.
- One color for the whole line.
- Loses upstream fixes in ProgressMeter.jl; the module is small enough to maintain.

## Future work

- If ProgressMeter.jl#369 is merged, moving back would mean re-adding leading spaces/separators to
  the elements; there is little reason to.
- A plain (no `\r`) mode for non-TTY outputs, e.g. one line every N percent, which could also replace
  `ProgressTxt`'s own percentage logic.
