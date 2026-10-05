# Progress line from ProgressMeter elements, no type piracy, optional vertical Courant number

> Status: **in progress**. Implemented and tested locally, draft PR waiting on ProgressMeter.jl#369. Replace the `ProgressMeter.speedstring(::AbstractFloat)` piracy and the global
> `FEEDBACK_*` refs by `ProgressMeter.Elements.AbstractProgressElement`s (ProgressMeter.jl PR #369) and add an
> optional `VerticalCourantNumber` element. Draft PR, depends on the unreleased ProgressMeter branch.

Date of initial draft: 2026-10-05

Base revision: e82b41763d831f2356ad4e661650f426ab9762ff

## Originating prompt

> Previously you've worked on https://github.com/timholy/ProgressMeter.jl/pull/369, given the current
> state of that branch can you create a PR in SpeedyWeather that would avoid the type piracy we're
> currently committing in order to customize the progress line printing and define an optional
> vertical Courant number element like PR1286

> Draft PR, [sources] to branch. Definition as given [max |σ̇| Δt/Δσ], but don't create a new flag, we
> want to be able to do something like
> `feedback = Feedback(spectral_grid, progress = Progress(..., elements = (..., VerticalCourant(...))))`

## Revision log

- PR #1286 turned out to be the vertical-advection fix (no Courant element); the Courant number is
  defined as max over the grid of |σ̇|Δt/Δσ instead.
- `Feedback` is created before the model and `Variables` exist and `Progress` needs `clock.n_steps`, so
  the user passes `elements` (not a finished `Progress`) to `Feedback`; `Progress` is still built in
  `initialize!`.
- `w` is radius-scaled inside the dynamical core (like divergence), so the element uses
  `Δt / vars.prognostic.scale[]`; a first version without it printed Cᵥ ≈ 3.6e5.

- Review (Opus, after the first push): CI failed in the Enzyme env because only the package and
  main test env had the `[sources]` entry, so ProgressMeter is now sourced from the branch in every env
  that resolves SpeedyWeather by path (docs, benchmark, benchmark/CUDA, test/GPU/*, test/differentiability,
  test/reactant). Bound elements printed the whole `Variables` (1 s to show a `Feedback`, 125k
  characters), so they now share `AbstractSimulationElement` with a compact `show`. The pirated
  `Base.show(::IO, ::ProgressMeter.Progress)` is removed too. `vertical_courant_number` reduces per
  layer (`maximum(abs, w; dims=1)`) instead of allocating two full-size temporaries per redraw. Docs
  section "Progress line" added to `how_to_run_speedy.md`.
- 2026-10-05, ProgressMeter.jl#369 changed after MarcMush's review ("In SpeedyWeather, can you update
  PR1298 to account for these changes?"): the elements now live in the exported submodule
  `ProgressMeter.Elements`, `ProgressStatus` is gone and `print_element(element, p)` takes two
  arguments, with the redraw time and finished state in `p.tcurrent`/`p.finished`. `SimulationSpeed`
  uses `p.tcurrent - p.tinit` instead of `status.elapsed`. Merging ProgressMeter master (#367) made
  `ProgressCore` parametric on the output type, so the Enzyme extension's `make_zero` now dispatches
  on `Type{<:ProgressCore}`.
- 2026-10-05, docs: "the usage of a custom Progress line is advanced and not something a first-time
  user should be worried about, can you move this section to Advanced -> Feedback?" The section
  "Progress line" moved from `how_to_run_speedy.md` to a new page `feedback.md` under Advanced, with
  a one-line pointer left in `how_to_run_speedy.md`.
- 2026-10-05: "can we group the progress elements? Could there be a local module call ProgressElements
  that contains the one from ProgressMeter.jl plus the ones we define? I want the user interface to be
  ProgressElements.MaximumWindSpeed() [...] Also remove the `showspeed`, `show_time`, `show_umax` kwargs
  from Feedback (they aren't considered public API), that should be now just controlled by `elements`".
  New submodule `ProgressElements` (`output/progress_elements.jl`, exported by SpeedyWeather, its names
  not) re-exports ProgressMeter's elements and holds SpeedyWeather's, `bind_element`, `speedstring`,
  `vertical_courant_number` and `default_elements()` (no argument). `Feedback.elements` defaults to
  `ProgressElements.default_elements()`; `showspeed`, `show_time`, `show_umax`,
  `show_temperature_range` are removed.
- 2026-10-05: "Can you rewrite the VerticalCourantNumber to avoid allocations? Can you use a vertical
  scratch vector for the `maximum` call? And given the smoothness of sigma, I'd probably not bother
  checking for both above and below but just estimate by using one [...] mark this points as comments".
  `VerticalCourantNumber` declares a `Vertical1D` scratch variable (`vars.scratch.vertical_velocity_maximum`,
  collected via `variables(::Feedback)` from its elements); `bind_element` stores a 1×nlayers reshape of
  it, a CPU buffer and Δσ on the CPU. `vertical_courant_number!` reduces with
  `maximum!(abs, w_max, w.data; init = false)`, copies to the CPU and uses only the interface below each
  layer. Allocation-free on CPU (tested); on CUDA the reduction leaves a 64 B device temporary from
  GPUArrays plus the host allocations of the kernel launch.

## Problem description

`src/output/feedback.jl` customizes the progress line by redefining
`ProgressMeter.speedstring(::AbstractFloat)` (type piracy: a method on a ProgressMeter function for a
ProgressMeter-owned call signature) and by smuggling the simulation time, time step, maximum wind speed
and temperature range through global `Ref`s (`FEEDBACK_DT_IN_SEC`, `FEEDBACK_TIME`, `FEEDBACK_UMAX`,
`FEEDBACK_TMIN`, `FEEDBACK_TMAX`). The diagnostics are computed every `check_iterations` steps whether or
not the line is drawn, and adding a diagnostic means editing `progress_string`.

## Background

ProgressMeter.jl PR #369 (branch `mc/modular-progress-callbacks` of `milankl-claude/ProgressMeter.jl`)
assembles the status line from `AbstractProgressElement`s with `print_element(element, p, status)`,
which is only called when the line is redrawn. It is open, conflicting and under design discussion, so
the element API may still change, and it is in no release. Hence a draft PR with a `[sources]` entry
pointing at the branch (Julia ≥ 1.11 only).

## Summary of changes

- New submodule `ProgressElements` (`output/progress_elements.jl`, the module is exported, its names
  are not, so the interface is `ProgressElements.MaximumWindSpeed()`). It re-exports the generic
  elements of `ProgressMeter.Elements` and defines SpeedyWeather's: `SimulationTime`,
  `SimulationSpeed`, `MaximumWindSpeed`, `TemperatureRange`, `VerticalCourantNumber`. Created unbound
  by the user, bound to `Variables` and `model` by `bind_element` in `initialize!(::Feedback, ...)`,
  so `Feedback` never holds a reference to `Variables` (Enzyme `make_zero` of the model) and the
  diagnostics are computed on redraw only.
- `Feedback.elements::Tuple = ProgressElements.default_elements()` (the previous default layout)
  replaces the options `showspeed`, `show_time`, `show_umax`, `show_temperature_range`, e.g.
  `Feedback(elements = (ProgressElements.default_elements()..., ProgressElements.VerticalCourantNumber()))`.
- Remove the `speedstring(::AbstractFloat)` method, `progress_string` and the `FEEDBACK_*` globals;
  `ProgressTxt` takes the time step from `model.time_stepping.Δt`.
- `progress!` no longer evaluates `max_speed`/`temperature_range` every few steps.
- `SpeedyWeather/Project.toml` and `SpeedyWeather/test/Project.toml`: `[sources] ProgressMeter` →
  fork branch.

## Testing and verification

New `SpeedyWeather/test/output/feedback.jl`: default elements reproduce the old line pieces, a custom
tuple with `VerticalCourantNumber` renders a finite number for `PrimitiveDryModel` and nothing breaks for
`Barotropic`, no `ProgressMeter.speedstring` method owned by SpeedyWeather remains.

## Documentation changes

Docstrings of the elements, `Feedback.elements`, new page "Feedback" (`feedback.md`, under Advanced) with the section "Progress line", linked from `how_to_run_speedy.md`;
CHANGELOG.

## Known limitations

- Depends on unreleased ProgressMeter API; CI on Julia 1.10 cannot resolve the `[sources]` entry until
  ProgressMeter is released and compat is raised.
- Courant number uses `Δt`, not the leapfrog span `2Δt`.

## Future work

Raise ProgressMeter compat and drop `[sources]` after release; more elements (e.g. CFL of horizontal
advection).
