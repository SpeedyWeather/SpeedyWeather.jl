# Unroll the tracer loops over compile-time tracer names, tracers in a NamedTuple

> Status: **in progress**. Part 1 (unrolled loops) and part 2 (tracers in a NamedTuple, breaking) are implemented and
> verified locally on the branch of draft PR #1409; part 2 is not committed yet.

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

> we are now back in the tracer work. I want to change the tracer dictionary to a tracer NamedTuple

> [answers to three design questions: tracers are set with a constructor keyword, like `greenhouse_gases`; the
> NamedTuple key is the name, so `Tracer` drops its `name` field; the change goes into the same PR #1409]

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
- 2026-11-17: part 2 drafted: replace the `Dict` registry with a NamedTuple fixed at model construction (the first item
  of Future work), in the same PR. See "Part 2" at the end.
- 2026-11-17: part 2 open questions resolved. Remove the type parameter `B` of `Tracer{B}` (`active::Bool`). It came in
  #970 without a stated reason, presumably for Reactant tracing the model. But a minimal `@compile` of a mutable model
  holding a NamedTuple of tracers gives the right result with a plain `Bool` field as well (Reactant 0.2.294, active
  and deactivated tracer). `tracers = nothing` is accepted and treated as `(;)`. Keep the Enzyme inactive rules.
  `tracer_active` keeps its `haskey` (it constant-folds for a NamedTuple) so that it is `false` for names that are not
  tracers, as the "tracer_active is inactive" test relies on.
- 2026-11-17: part 2 implemented as planned. Additionally: `tracer_names` (keys of the tracers, `()` for `nothing`)
  used by `check_tracers`; the frozen-tracer test also checks that a simulation assembled from the variables of a model
  with other tracers errors. See Results of part 2.

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
  fixed at initialization. That would also make `tracer.active` a plain field load. Taken up as part 2 below.
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

---

# Part 2: tracers in a NamedTuple

Date of initial draft: 2026-11-17

Base revision: 174afe29 (`mg/unroll-tracer-loops`, part 1 complete)

## Problem description

After part 1, the time stepping no longer iterates `model.tracers::Dict{Symbol, Tracer}`, but it still looks up every
tracer in it, once per tracer and site per time step: `tracer_active` does `haskey` and `getindex` on a `Dict` whose
value type `Tracer` is abstract. That leaves three loose ends:

- `tracer_active` must be marked inactive for Enzyme to keep the `Dict` lookups out of the differentiated code, and
  `tracer.active` is still a builtin `getfield` on an abstract type.
- The set of tracers is mutable after initialization although the variables are not. Part 1 covers this with
  `check_tracers` (tracers added later raise an error) and the rule that deleted tracers are frozen.
- It is inconsistent with the greenhouse gases, the other named collection of components in the model, which are a
  NamedTuple set when the model is constructed: `PrimitiveWetModel(spectral_grid; greenhouse_gases = (; co2 = CO2(...)))`.

The aim is to store the tracers as a NamedTuple in a type parameter of the model, fixed when the model is
constructed. Then the tracer names are in the types of both the model and the variables, `tracer_active` reduces to
a field load that constant-folds, and tracers can no longer be added or deleted once the variables exist.

## Background

- All four models are `@parameterized @kwdef mutable struct`s with a concretely typed `tracers::TRACER_DICT` field.
  With a type parameter `tracers::TR`, assigning a NamedTuple with other names to `model.tracers` fails, so
  `add!`/`delete!` cannot work in place any more. Hence the constructor keyword (design decision 1).
- Greenhouse gases use the NamedTuple key as the name: `variables(::Pair{Symbol, <:AbstractGreenhouseGas}, model)`
  passes the key on, and `variables(nt::NamedTuple, model)` in `variables/variables.jl` already iterates the pairs of
  any NamedTuple component. Tracers can use the same mechanism (design decision 2).
- `reconstruct(model, parameters)` (used by `initialize!(model, ps)`) rebuilds the model with
  `ConstructionBase.setproperties`, which keeps unpatched fields and their types, so the new type parameter needs
  no special handling there. `tracers` is not a `@component`, as `Tracer` has no parameters.
- A NamedTuple supports indexing by `Symbol`, so `model.tracers[:abc].active` keeps working, as does
  `model.tracers.abc.active`.
- `SpeedyWeather` is at `0.25.0-DEV`, an unreleased breaking version, so no version bump is needed.

## Summary of changes

All in `SpeedyWeather/` plus docs, tests and changelog.

### 1. `Tracer` without a name (`dynamics/tracers.jl`)

```julia
@kwdef mutable struct Tracer <: AbstractTracer
    active::Bool = true
end

# `Tracer(:abc)` was the old way to define a tracer, give a migration hint
Tracer(name::Symbol; kwargs...) = throw(ArgumentError(
    "Tracers are not named any more, define them through the `tracers` keyword when constructing the model, " *
    "e.g. `ShallowWaterModel(spectral_grid; tracers = (; $name = Tracer()))`. (De)activate them with " *
    "`activate!(model, :$name)` or `deactivate!(model, :$name)`."))
```

All old code calls `Tracer(:name)` (in `add!`, `delete!`, `activate!` and `deactivate!` alike), so this one error
covers every migration case. `Tracer(active = false)` defines a tracer that starts deactivated.

The type parameter `B` (the type of `active`, from #970) is removed: `active` is always a `Bool`. A minimal Reactant
check traced a model with a NamedTuple of such tracers correctly (see revision log); the full check is the Reactant
tracer test below.

### 2. Models (`models/barotropic.jl`, `shallow_water.jl`, `primitive_dry.jl`, `primitive_wet.jl`)

A new type parameter `TR, # NamedTuple of Tracers` and `tracers::TR = (;)` in place of
`tracers::TRACER_DICT = TRACER_DICT()`, at the same position, under `# VARIABLES`.

`tracers = nothing` is accepted and means the same as `(;)`. It cannot be converted to `(;)` on construction:
`@kwdef` passes the keyword straight to the default constructor, which infers `TR` from it, and the keyword form
`Model(; spectral_grid, tracers)` bypasses `Model(spectral_grid; kwargs...)`. So `nothing` is stored, and the
functions that read `model.tracers` treat it as empty: `variables(nothing, model)` already returns `()`, and
`tracer_active`, `set_tracers_active!` and `check_tracers` get a `Nothing` method each (below). This is how
`greenhouse_gases = nothing` works too.

### 3. Registry functions (`dynamics/tracers.jl`)

Removed: `TRACER_DICT`, both `add!` methods for tracers, `delete!(model, ::Tracer)`,
`variables(::TRACER_DICT, model)`, `_activate!`.

New or changed:

```julia
# (de)activate by name, on the model or the simulation
activate!(simulation::AbstractSimulation, names::Symbol...) = activate!(simulation.model, names...)
activate!(model::AbstractModel, names::Symbol...) = set_tracers_active!(model.tracers, true, names...)
deactivate!(simulation::AbstractSimulation, names::Symbol...) = deactivate!(simulation.model, names...)
deactivate!(model::AbstractModel, names::Symbol...) = set_tracers_active!(model.tracers, false, names...)

function set_tracers_active!(tracers::NamedTuple, value::Bool, names::Symbol...)
    for name in names
        haskey(tracers, name) || throw(ArgumentError("No tracer :$name in model.tracers, tracers are $(keys(tracers))"))
        tracers[name].active = value
    end
    return tracers
end
set_tracers_active!(::Nothing, value::Bool, names::Symbol...) = set_tracers_active!((;), value, names...)

# the NamedTuple key is the name of the tracer, as for greenhouse gases
variables(name_tracer::Pair{Symbol, <:Tracer}, model::AbstractModel) =
    variables(name_tracer.second, name_tracer.first, model)
variables(::Tracer, name::Symbol, model::AbstractModel) = (...)  # the current body, `T.name` replaced by `name`

# a field load: `name` is a compile-time constant in `foreach_active_tracer`, so `haskey` and
# `getfield` constant-fold; `false` for names that are not tracers
tracer_active(model::AbstractModel, name::Symbol) = tracer_active(model.tracers, name)
tracer_active(tracers::NamedTuple, name::Symbol) = haskey(tracers, name) && getfield(tracers, name).active
tracer_active(::Nothing, name::Symbol) = false
```

`foreach_active_tracer` stays as it is (unrolled over the tracer names in the type of the variables).

`check_tracers(vars, model)` stays in `initialize!(simulation)` but checks something else: that `model.tracers` is a
NamedTuple of `Tracer`s (or `nothing`, as `(;)`) whose names equal the tracer names in `vars`. With the constructor keyword, these can only
differ if a simulation is assembled by hand from the variables of one model and another model. The more likely
mistake it catches is `tracers = (Tracer(), Tracer())`, a `Tuple`: its elements have no names, so `Variables`
silently allocates nothing for them (`variables(::Any) = ()`), and without the check the tracers would be ignored
without any error.

### 4. Enzyme extension

Keep `EnzymeRules.inactive` for `tracer_active` and `check_tracers`. A `Bool` field load is inactive anyway, so the
rule is not needed for correctness, but it costs nothing and states the intent. Update the comment above them: the
registry is a NamedTuple of `active` flags that carries no derivative.

### 5. Changelog

Rewrite the line of this PR (and replace its `#NNN` placeholder with `#1409`) to mention the breaking change: tracers
are defined with the `tracers` keyword of the model constructor as a NamedTuple, `Tracer` has no name, and
`add!`/`delete!` for tracers are removed.

## Behaviour changes

| Situation | Before (part 1) | After |
|-----------|-----------------|-------|
| Define tracers | `add!(model, Tracer(:abc), Tracer(:co2))` | `ShallowWaterModel(spectral_grid; tracers = (; abc = Tracer(), co2 = Tracer()))` |
| Start deactivated | `Tracer(:abc, active = false)` | `(; abc = Tracer(active = false))` |
| Delete a tracer | `delete!(model, Tracer(:abc))` | construct a new model without it, or deactivate it |
| (De)activate | `activate!(model, Tracer(:abc))` | `activate!(model, :abc)` (also on `simulation`) |
| Look up | `model.tracers[:abc]` | `model.tracers[:abc]` or `model.tracers.abc` |
| `Tracer(:abc)` | defines a tracer | `ArgumentError` with a migration hint |
| Tracer added after `initialize!(model)` | `ArgumentError` from `run!` | not possible |
| Tracer deleted after `initialize!(model)` | frozen | not possible |
| `tracers = (Tracer(), Tracer())` (a `Tuple`) | — | `ArgumentError` from `run!`/`initialize!(simulation)` |
| `tracers = nothing` | — | no tracers, as `(;)` |
| `Tracer{B}` with a non-`Bool` `active` | possible, fails in `tracer_active` | not possible, `active::Bool` |
| Results of active and deactivated tracers | | unchanged, bitwise |

## Testing and verification

Unit tests:

- `test/dynamics/tracer_advection.jl`:
  - "add/delete" testset becomes "tracers keyword": models with no, one and three tracers have the right
    `keys(model.tracers)` and tracer variables; `Tracer(:abc)` throws an `ArgumentError`; a `Tuple` of tracers makes
    `run!` throw an `ArgumentError`; a model with `tracers = nothing` runs, has no tracer variables, and
    `activate!(model, :abc)` throws an `ArgumentError`; `reconstruct(model, parameters(model))` keeps `model.tracers`.
  - (de)activate testset: by `Symbol` on the model and on the simulation; `model.tracers[:abc].active` reflects it;
    an unknown name throws an `ArgumentError`.
  - "deactivated or deleted ... added after initialization" testset: keep the active and the deactivated tracer,
    drop the deleted and the late tracer (both no longer possible).
  - primitive equation testset: construct with the keyword.
- `test/dynamics/dispatch.jl`: construct the tracer case with the keyword. Its JET test now also covers
  `tracer_active` as a field load.

Other tests that only need the new construction: `test/differentiability/parameters.jl` ("tracer_active is
inactive", the populated case becomes a model constructed with two tracers), `test/differentiability/forward_mode.jl`
(tracer case), `test/reactant/test_correctness.jl` (`test_model` passes the tracers to the model constructors via
`create_cpu_model`/`create_reactant_model`).

Checks during development (results go in the PR):

- **Bitwise regression against main:** BarotropicModel, ShallowWaterModel and PrimitiveWetModel with two tracers,
  one deactivated after 10 of 20 steps; main defines them with `add!`, the branch with the keyword. All prognostic,
  grid and tendency arrays must be identical.
- **Differentiability:** `parameters.jl` and `forward_mode.jl` pass (Enzyme 0.13.220, Julia 1.12).
- **Reactant:** the PrimitiveDryModel tracer case of `test_correctness.jl` passes. This is also the check on the full
  model that `Tracer` needs no type parameter.
- **Docs:** the code blocks of `docs/src/tracers.md` run as a script (the full docs build is not needed for this).

### Results (2026-11-17, Julia 1.12, aarch64 macOS)

- **Unit tests** (`--check-bounds=yes`): `test/dynamics/tracer_advection.jl` 6383 passed; `test/dynamics/dispatch.jl`
  all 3 JET tests passed; `test/parameters.jl` 21 and `test/dynamics/simulation_constructor.jl` 10 passed.
- **Bitwise regression:** all 30 arrays (prognostic, grid and tendency tracers, vorticity, `u`) of BarotropicModel,
  ShallowWaterModel and PrimitiveWetModel with two tracers, one deactivated after 10 of 20 steps, are bitwise equal
  to b06af911 (main before part 1, tracers defined with `add!`). Two runs of the branch are bitwise equal too. The
  default `RandomVelocity` initial conditions of the BarotropicModel are randomly seeded (`seed = 0`), so the
  comparison fixes `seed = 1`; without that, BarotropicModel runs differ between any two runs.
- **Reactant** (0.2.308, Enzyme 0.13.220, the CI setup): BarotropicModel without tracers 8 of 8, PrimitiveDryModel with
  a tracer 22 of 22 (as in part 1); the tracer differs between CPU and Reactant by at most about 1e-5 (relative). So
  the tracer needs no type parameter `B`. The Reactant Correctness CI of part 1 (174afe29) failed at its first
  case: `test_model` called `add!(model)` without tracers for models without tracers, which has no method. Passing
  the tracers to the model constructor fixes that.
- **Differentiability** (Enzyme 0.13.220, Enzyme_jll 0.0.310): `parameters.jl` passes, including "tracer_active is
  inactive". In `forward_mode.jl` the BarotropicModel case passes (9 of 9), and so does the BarotropicModel case
  with a tracer (10 of 10), the forward/reverse check of the tracers. The ShallowWaterModel and PrimitiveWetModel cases
  error with an `IllegalTypeAnalysisException` in `bernoulli_grid_potential!` (the ShallowWater and the primitive
  equation method), code the tracers don't touch. With the same environment at part 1's HEAD (174afe29) both cases
  give the same error at the same line, so it is not caused by this change; part 1's 49 of 49 were with an
  environment that is gone, presumably with other versions of indirect dependencies.
- **Docs:** all 9 code blocks of `docs/src/tracers.md` run (with the plotting lines removed).
- **Formatting:** Runic reports no changes.

## Documentation changes

- `docs/src/tracers.md`: rewrite "Add/delete tracers" as "Define tracers" with the constructor keyword (several
  tracers, the key as name, a new model to change the set of tracers), "(De)activate tracers" with `Symbol`s and
  `Tracer(active = false)`, and the visualisation example. Note that `OutputVariables.Tracer(:abc)` (output) still
  takes a name, unlike the model's `Tracer()`.
- Docstrings of `Tracer`, `activate!`, `deactivate!`, `tracer_active` and `check_tracers`.

## Known limitations

- Breaking for all user code that defines or (de)activates tracers. Old code fails at `Tracer(:name)` with an
  error that shows the new way, so nothing changes silently.
- Changing the set of tracers needs a new model and simulation. Before, adding needed a new simulation anyway, but
  deleting during a simulation was possible (the tracer was frozen); deactivating has the same effect.
- Each set of tracer names gives a distinct model type. The variables type already contains the tracer names, so this
  does not add compilation in practice.
- Two different `Tracer`s: the model's (no name) and `OutputVariables.Tracer(:abc)` (with a name). The docs say so.
- `tracers = nothing` is stored as `nothing`, not converted to `(;)` (see section 2), so `model.tracers` is `nothing`
  for such a model. Everything that uses the tracers treats both the same.
- Part 1's limitation that `tracer_active` asserts a `Bool` no longer applies: `active` is always a `Bool`.
- Unrelated to this change, but found while verifying it: with Enzyme 0.13.220 and the current indirect dependencies,
  the ShallowWaterModel and PrimitiveWetModel cases of `test/differentiability/forward_mode.jl` fail locally with an
  `IllegalTypeAnalysisException` in `bernoulli_grid_potential!`, at part 1's HEAD as well (see Results).

## Open questions

Resolved 2026-11-17 (see revision log): remove `B`, keep the Enzyme inactive rules, accept `tracers = nothing` as
`(;)`.
