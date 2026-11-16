# Unroll the tracer loops over compile-time tracer names

> Status: **in progress**. Implemented on branch `mg/unroll-tracer-loops` and verified locally (see "Results" below); no PR yet.

Date of initial draft: 2026-11-16

Base revision: b06af911

## Originating prompt

> investigate the warning I get in the Reactant / Enzyme tests
>
> [CI log of the Reactant Correctness workflow: `Warning: .../time_stepping/time_integration.jl:122:0:
> in function preprocess_julia_update_prognostic__576479 ...: Enzyme: Cannot deduce type of load
> %.sroa.8.16.copyload22 = load i64, ptr addrspace(11) %88, ... assumed {[-1]:Pointer}`, twice each
> for the BarotropicModel differentiability tests with `SpectralTransform` and `MatrixSpectralTransform`]

> yes, draft the plan for unrolling the tracer loops

> go with your recommendations, implement it

## Revision log

- 2026-11-16: initial draft.
- 2026-11-16: signed off with the recommended answers to all three open questions: remove the custom
  `iterate(::TRACER_DICT)` rules, use one generic `foreach_active_tracer` driver, and add a Reactant correctness case
  with a tracer if it passes on main.
- 2026-11-16: the three `transform!` sites share one `transform_tracer!` that always passes `model` to
  `get_prognostic_step`. The primitive-equation site passed no model before, but the only model-specific
  `which_prognostic_step` methods are for `TwoDModels`, so for primitive equation models the 4-argument call falls back
  to the 3-argument one and nothing changes.
- 2026-11-16: Reactant correctness test with tracers. `compare_arrays` only compares top-level arrays, so tracers (in
  their own namespace) would not have been compared. Added `compare_tracers` for the `tracers`/`grid_tracers`
  namespaces, a `tracers` keyword to `test_model`, and a check that the tracers are advected, not merely frozen.
  Cases: BarotropicModel and PrimitiveDryModel with one tracer each; the latter covers vertical advection and the
  primitive-equation transform and diffusion.
- 2026-11-16: the same Reactant tracer test also passes on main. Reactant traces the time step once, so on main the
  `Dict` loop (which no longer exists on this branch) runs in plain Julia at compile time and the compiled program
  already contains one block of operations per active tracer. On this branch the `@generated` `foreach_active_tracer`
  does that unrolling in Julia, and only the `tracer_active` lookup is evaluated during tracing. Both compile to the
  same Reactant program, so under Reactant the change adds test coverage rather than fixing anything.
- 2026-11-16: differentiating `time_step!` in forward and reverse mode on Julia 1.12 with Enzyme 0.13.220 prints many
  "Cannot deduce type of store `store i64 -1`" warnings in `_transform!` and `dynamics_tendencies!` (and, with a
  tracer, in `horizontal_advection_tracer!`, which calls the same transforms). They appear without tracers too, so
  they are unrelated to this change. CI runs these tests on Julia 1.10 and 1.11 and doesn't show them.
- 2026-11-16: only the PrimitiveDryModel tracer case was kept in CI (the BarotropicModel one was dropped in commit
  f5b9047e). Its tracer was a σ-independent blob on the equator, where the default (Jablonowski–Williamson) wind is
  zero, and "advected" was only checked as `!=`, which diffusion alone would satisfy. Changed to a blob at 45˚N in
  the midlatitude jet that increases linearly with σ, so vertical advection acts from the first step, and to require
  a relative change of more than 0.05 over the compared steps. See Results for the measurements behind this.

## Problem description

The time stepping loops over the tracer registry `model.tracers::Dict{Symbol, Tracer}` at runtime in 9 places:

| File | Function | Line |
|------|----------|------|
| `time_stepping/time_integration.jl` | `update_prognostic!` | 136 |
| `time_stepping/transform.jl` | `transform!` (Barotropic, ShallowWater, PrimitiveEquation) | 42, 100, 174 |
| `dynamics/horizontal_diffusion.jl` | `horizontal_diffusion!` (Barotropic, ShallowWater, PrimitiveEquation) | 271, 308, 357 |
| `dynamics/tendencies.jl` | `tracer_advection!` | 804 |
| `dynamics/vertical_advection.jl` | `vertical_advection!` | 36 |

Each loop has the form `for (name, tracer) in model.tracers; tracer.active && do_something(vars.….tracers[name]); end`.
Iterating a `Dict` whose value type `Tracer` is abstract (`Tracer{B}` is a `UnionAll`) is a problem for Enzyme in every
differentiated time step:

1. **Warnings since Enzyme_jll 0.0.310.** EnzymeAD/Enzyme#3712 now always warns when Enzyme cannot deduce the type of
   an active load and falls back to guessing a pointer. With Enzyme.jl 0.13.220 / Enzyme_jll 0.0.310 every CI run
   prints these warnings: 4 in the Reactant Correctness workflow (`update_prognostic!`) and 12 in the Enzyme workflow
   (`update_prognostic!`, `vertical_advection!`, `move_prognostic_grid_variables_back!`). Main on Enzyme 0.13.219 /
   jll 0.0.307 showed none; the fallback itself is unchanged, only the warning is new.
2. **Custom Enzyme rules as a workaround.** `SpeedyWeatherEnzymeExt.jl` defines forward and reverse rules for
   `iterate(::TRACER_DICT, ...)`, marked as temporary. They exist only because forward mode fails with an internal
   error compiling `iterate(::Dict)` inside the time stepping.
3. **Forward mode with tracers is broken.** Forward mode over a populated registry throws `EnzymeNoShadowError`
   (`@test_broken` in `test/differentiability/parameters.jl`). Forward-mode differentiation of a model with
   tracers is therefore not possible.
4. **Inconsistent with the rest of the time stepping.** The atmospheric, ocean and land variables are already unrolled
   over compile-time names from the `Variables` type (`update_prognostic_atmosphere!`,
   `update_prognostic_namespaces!`, `initialize_steps_unrolled!`, `move_grid_tracers_back!`). A TODO at
   `time_integration.jl:134` asks for the same for the tracers.

The warnings themselves are benign (see Background). The aim is to remove the runtime `Dict` iteration from the
differentiated code. That removes the warnings, makes the custom `iterate` rules unnecessary, and should make forward
mode with tracers work.

## Background

**What the warning is.** In the IR of `update_prognostic!` (dumped with `Enzyme.Compiler.DumpPreEnzyme[] = true`,
Julia 1.12, Enzyme 0.13.220, BarotropicModel T32, reproduced locally on aarch64):

- `%51 = call @julia_iterate_…(model.tracers, i64 %state)` is `iterate(::Dict{Symbol, Tracer}, state)`. Because `Tracer`
  is abstract, it returns a *boxed* `Union{Nothing, Tuple{Pair{Symbol, Tracer}, Int}}`.
- The flagged `%.sroa.8.16.copyload22 = load i64 … %74` reads byte offset 16 of that tuple, i.e. the `Int` iteration
  state, which only feeds the next `iterate` call. The pointer guess is therefore harmless: no derivative flows through
  it, and the test's Reactant-vs-Enzyme gradient comparison (rtol 1e-3) passes.
- The state goes around the loop (load → phi → argument to `iterate` → load), and the custom `iterate` rule returns no
  shadow. Enzyme's activity analysis therefore cannot prove the state inactive. This is the most likely reason only this
  load is flagged.
- `tracer.active` compiles to a builtin `jl_f_getfield` call on an abstract `Tracer`.

The `vertical_advection!` warnings in the Enzyme CI have the same shape (an `i64` unpacked from the boxed `iterate`
result) and almost certainly the same cause. The `move_prognostic_grid_variables_back!` warning is a different load
(`size` of an array in Base `array.jl:194` on Julia 1.10), has nothing to do with the tracers, and is out of scope.

**Type stability is not the issue.** JET's `report_opt` on `time_step!` of a BarotropicModel *with* a tracer reports no
runtime dispatch on main. All tracer variables have the same type, so `vars.prognostic.tracers[name]` infers
concretely even with a runtime `Symbol`.

**Relation between `model.tracers` and the variables.** `Variables(model)` allocates variables per tracer in
`model.tracers` at `initialize!(model)`, in the namespaces `prognostic.tracers`, `grid.tracers`, `tendencies.tracers` and
`tendencies.grid_tracers`. `_namespace_names(T, :tracers)` / `tracer_tendency_names(vars)` give these names at compile
time. The docs (`docs/src/tracers.md`) state that tracers must be added before `initialize!(model)`. The `active` flag can
be toggled at any time and must stay a runtime check.

**Reactant.** `Tracer.active` is a host `Bool` (numbers are not tracked by `@compile`), so the current loop is resolved at
trace time. The same will hold for the unrolled version. No Reactant test currently uses tracers.

## Summary of changes

All changes are in `SpeedyWeather/` (`src/` and `ext/`).

### 1. Two helpers in `dynamics/tracers.jl`

```julia
"""$(TYPEDSIGNATURES)
Whether the tracer `name` is in `model.tracers` and active. A tracer deleted after initialization is
treated as inactive (frozen). Host-side lookup in the tracer dictionary, inactive for Enzyme."""
tracer_active(model::AbstractModel, name::Symbol) = tracer_active(model.tracers, name)
tracer_active(tracers::TRACER_DICT, name::Symbol) = (haskey(tracers, name) && tracers[name].active)::Bool

"""$(TYPEDSIGNATURES)
Call `f(Val(name), vars, model, args...)` for every tracer `name` that has variables in `vars` and is
active in `model.tracers`. Unrolled at compile time over the tracer names in the type of `vars`, so
that no `Dict` is iterated at runtime."""
@generated function foreach_active_tracer(
        f::F, vars::Variables{Po, G, T}, model::AbstractModel, args...
    ) where {F, Po, G, T}
    calls = [
        :(tracer_active(model, $(QuoteNode(name))) && f(Val($(QuoteNode(name))), vars, model, args...))
            for name in _namespace_names(T, :tracers)
    ]
    return Expr(:block, calls..., :(return nothing))
end
```

For models without tracers the generated body is just `return nothing`.

### 2. Replace the 9 loops

Each loop body moves verbatim into a named per-tracer function with a `::Val{name}` first argument, so that
`getfield(vars.….tracers, name)` constant-folds. Each site keeps its exact `get_*_step` call, including whether it passes
`model`. Example for `update_prognostic!`:

```julia
# in update_prognostic!, replacing the loop and its TODO
foreach_active_tracer(update_prognostic_tracer!, vars, model, clock, time_stepping, model.implicit, scale)

function update_prognostic_tracer!(::Val{name}, vars, model, clock, time_stepping, implicit, scale) where {name}
    var = getfield(vars.prognostic.tracers, name)
    tendency = getfield(vars.tendencies.tracers, name)
    return update_prognostic!(var, tendency, clock, time_stepping, implicit, model, scale)
end
```

| Site | Per-tracer function |
|------|---------------------|
| `update_prognostic!` | `update_prognostic_tracer!` |
| `transform!` ×3 | `transform_tracer!` (one method per call signature currently used) |
| `horizontal_diffusion!` ×3 | `horizontal_diffusion_tracer!` (the three bodies are identical) |
| `tracer_advection!` | `horizontal_advection_tracer!` |
| `vertical_advection!` | `vertical_advection_tracer!` |

Each call stays at the same position in its function, so the order relative to other operations is unchanged. This
matters e.g. for the tracer transform before `move_prognostic_grid_variables_back!` in the primitive `transform!`.

### 3. Fail early for tracers added after initialization

Add `check_tracers(variables, model)`, called in `initialize!(::AbstractSimulation)` (`models/simulation.jl`). It runs on
the host once per `run!` and throws an `ArgumentError` naming every tracer in `model.tracers` that has no variables in
`vars`. Without the check, such a tracer would be silently skipped by the unrolled loops. Today it fails with a missing
field error inside the first time step.

### 4. Enzyme extension

- Add `EnzymeRules.inactive(::typeof(SpeedyWeather.tracer_active), args...) = nothing`, and the same for
  `check_tracers`, so Enzyme never differentiates into the `Dict` lookups.
- Remove the forward/augmented-primal/reverse rules for `iterate(::TRACER_DICT, ...)` and their comment block, as
  nothing in the time stepping iterates the registry any more (see Open questions).

### 5. Versioning and changelog

`SpeedyWeather` is at `0.25.0-DEV`, so no version bump is needed. Add a line to `CHANGELOG.md` under `## Unreleased`.

## Behaviour changes

| Situation | Before | After |
|-----------|--------|-------|
| Active tracer | advanced | advanced (same arithmetic, bitwise identical) |
| Deactivated tracer | frozen | frozen |
| Tracer deleted after `initialize!(model)` | frozen (not iterated) | frozen (`tracer_active` is false) |
| Tracer added after `initialize!(model)` | missing field error in the first time step | `ArgumentError` from `run!`/`initialize!(simulation)` |
| Order in which tracers are processed | `Dict` order | `Variables` order (= `Dict` order at initialization). Tracers are independent, so results do not depend on it |

## Testing and verification

Unit tests (main test suite):

- `test/dynamics/tracer_advection.jl`: the existing `add!`/`delete!`, `activate!`/`deactivate!` and primitive-model
  testsets must pass unchanged. Add a testset with two tracers on a BarotropicModel:
  - a deactivated tracer stays bitwise frozen while the active one changes;
  - a tracer deleted after `initialize!` stays frozen;
  - a tracer added after `initialize!` makes `run!` throw an `ArgumentError`.
- `test/dynamics/dispatch.jl`: add a BarotropicModel-with-a-tracer case to the JET runtime-dispatch test. It already
  passes on main; it guards the new generated code.

Differentiability tests (extended suite, `test/differentiability/`):

- `parameters.jl`: replace the "Tracer registry iterate rules" testset (including its `@test_broken`) with a test that
  `tracer_active` contributes no derivative in forward and reverse mode with a populated registry.
- `forward_mode.jl`: add a BarotropicModel-with-one-tracer case to the forward-vs-reverse dot-product test, seeding the
  tracer too. Expected to work after the change; this is the main check for item 3 of the problem description. If it
  does not, record why in the revision log.

Checks done once during development (results go in the PR):

- **Bitwise regression.** A script comparing main and the branch for BarotropicModel, ShallowWaterModel and
  PrimitiveWetModel with two tracers (one deactivated halfway) over 20 steps. All prognostic and grid tracer arrays must
  be identical.
- **Warnings gone.** Re-run the CPU reference gradient of `test/reactant/differentation.jl` with Enzyme 0.13.220. There
  should be no "Cannot deduce type of load" warning from `update_prognostic!`. In the PR's CI logs, the Reactant
  workflow should have none and the Enzyme workflow only the `move_prognostic_grid_variables_back!` ones.
- **Reactant.** Run one time step of a BarotropicModel with a tracer on CPU and with Reactant and compare them, as in
  `test/reactant/test_correctness.jl`. Run it on main first: if it already fails there, record that as a pre-existing
  limitation; otherwise add it to the Reactant correctness tests (see Open questions).

No benchmark run is planned. The removed per-step work is one iteration over a (usually empty) `Dict` per site,
negligible compared with the transforms.

### Results (2026-11-16, Julia 1.12.7, aarch64 macOS)

- **Bitwise regression:** identical. All 42 arrays (prognostic, grid and tendency tracers, vorticity, `u`) of
  BarotropicModel, ShallowWaterModel and PrimitiveWetModel with three tracers are bitwise equal to main after 20 steps.
  One tracer was deactivated after 10 steps and one deleted after initialization. Two runs on main were identical too,
  so the comparison is deterministic.
- **Unit tests** (`--check-bounds=yes`):
  - `test/dynamics/tracer_advection.jl`: 6377 passed, including the 8 new ones.
  - `test/dynamics/dispatch.jl`: all 3 JET tests passed, including the new tracer case.
- **Differentiability** (Enzyme 0.13.220):
  - `parameters.jl` (including "tracer_active is inactive") and `nested_ad.jl` pass.
  - The new tracer case of `forward_mode.jl` passes. Forward mode with a populated tracer registry works now, and the
    tracer block satisfies the forward/reverse dot-product identity. 49 of 49 passed.
- **Warnings:** the CPU reference gradient of `test/reactant/differentation.jl` compiles without any
  "Cannot deduce type of load" warning; main with the same Enzyme gave 2.
- **Reactant** (0.2.308, versions as in CI): the new correctness cases pass. BarotropicModel with a tracer: 13 of 13.
  PrimitiveDryModel with a tracer: 22 of 22. In both, the tracer differences between CPU and Reactant are comparable to
  those of vorticity, within the default tolerances.
- **Choice of the Reactant test tracer** (PrimitiveDryModel, CPU, test setup: 10-day spin-up, then the 100 compared
  steps; variants with advection switched off by tracer name). The relative change of the spectral tracer over the
  100 compared steps:

  | Tracer | Advected | Diffusion only | Vertical advection's share |
  |---|---|---|---|
  | Equator, σ-independent (before) | 0.38 | < 0.0001 | 0.10 |
  | 45˚N, ∝ σ (now) | 0.74 | 0.0016 | 0.59 |

  The threshold of 0.05 therefore separates advection from diffusion by more than an order of magnitude either way.
  With the new tracer the Reactant case still passes 22 of 22; the tracer differs between CPU and Reactant by at
  most about 1e-5 (mean about 2e-7), far within the norm-relative tolerance of 1e-3.
- **Over- and undershoots.** The tracer has no positivity limiter, so spectral advection over- and undershoots at
  T32. With the new tracer, the minimum reaches −0.61 at day 10 in the lowest layer at 69˚N, where the baroclinic wave
  breaks; with horizontal advection only it is −0.14. This comes from the dispersive default
  `CenteredVerticalAdvection`. With `UpwindVerticalAdvection` the minimum is −0.18. A constant tracer stays exactly 1
  and c = σ stays within 0.00–1.16, so the vertical advection of tracers is consistent. Temperature uses the same
  scheme. For a CPU-vs-Reactant comparison the undershoots don't matter.

## Documentation changes

- `docs/src/tracers.md`, "Add/delete tracers": say that `run!` now raises an error for tracers added after
  `initialize!(model)`, and that a tracer deleted after initialization is frozen.
- Docstrings for `tracer_active` and `foreach_active_tracer`, and a short comment at each call site in the existing
  "unrolled per name" style.

## Known limitations

- The `move_prognostic_grid_variables_back!` warning in the Enzyme CI remains. It has a different cause, not inspected.
- Tracers still cannot be added after `initialize!(model)`; the change only turns a confusing late error into a clear
  early one.
- With Reactant, `active` is read when the program is compiled, as on main. Every `run!` compiles the time stepping
  again (`@compile` in the extension's `time_stepping!`), so (de)activating a tracer between `run!` calls works. A
  switch during a run, e.g. from a callback, does not reach the compiled program.
- Forward/reverse differentiation of `time_step!` on Julia 1.12 prints "Cannot deduce type of store" warnings from the
  spectral transforms (see revision log), unrelated to the tracers.
- `tracer_active` asserts a `Bool`. A user-constructed `Tracer{B}` with a non-`Bool` `B` (the parameter introduced in
  #970) would fail there; nothing in the package creates one.
- Compile time grows by one specialization per tracer name and site. This is negligible for a handful of tracers.

## Future work

- Give the registry a concrete value type (`Dict{Symbol, Tracer{Bool}}`) or replace it with an immutable collection
  fixed at initialization. That would also make `tracer.active` a plain field load.
- Investigate the remaining `move_prognostic_grid_variables_back!` warning and the Julia 1.12 store warnings in the
  spectral transforms.

## Open questions

Resolved at sign-off (see revision log): all three recommendations taken.

1. **Remove the custom `iterate(::TRACER_DICT)` rules?** Recommended. After this change nothing in the package iterates
   the registry in differentiated code, and the rules are marked temporary. The cost: user code that iterates
   `model.tracers` inside `autodiff` in forward mode would hit Enzyme's internal error again (the populated-registry
   case is broken even with the rules).
2. **One generic `foreach_active_tracer` driver or one `@generated` function per site?** The generic driver is
   recommended (9 sites, one generated function). Per-site generated functions, as in `move_grid_tracers_back!`, are
   more explicit but repeat the name handling nine times.
3. **Add a Reactant correctness case with a tracer to CI?** Recommended if it passes on main; it costs roughly a minute
   of compile time in a ~1 h workflow.
