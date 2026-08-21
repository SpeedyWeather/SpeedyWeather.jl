# Reducing compile time of `initialize!(model)` and `run!(simulation)`

> Status: **in progress**. Changes 1 and 2 implemented and verified: total `initialize!` 34.0 s → 22.6 s. Change 3 (precompile workload) still outstanding and currently out of scope. `run!` profiled 2026-08-21: no structural fix found, see "Profiling `run!`".

Date of initial draft: 2026-08-20

Base revision: f89bce984b028f4aa9f62622accc32aba157d6a8

## Originating prompt

> We want to reduce the compile time to set up a simulation with SpeedyWeather. Creating a model is fast enough for now, but `initialize!(model)` runs essentially through 3 steps: initialize components, allocate variables and initialize variables to the initial conditions. For the primitive wet model at default resolution can you profile how long each step and substeps within take and propose some changes that could reduce compile time here?

## Revision log

- 2026-08-20: initial draft.
- 2026-08-20: implemented Change 1 (`set!` type erasure). Step 3 22.7 s → 13.8 s;
  total `initialize!` 31.4 s → 22.6 s. Required a correction to the original design: erasing
  the kwargs into a `Vector` of pairs was **not** sufficient on its own — see below.
- 2026-08-20: implemented Change 2 (variable collection via `Vector` instead of `Tuple`).
  `Variables(model)` 7.62 s → 4.63 s (−39 %); total `initialize!` 34.0 s → 31.4 s.
  Verified all 113 allocated variable arrays are bitwise identical before/after.
  Changes 1 and 3 deliberately left for separate commits.
- 2026-08-21: profiled `run!(simulation)` on request (no changes made). Found no dominant
  mechanism analogous to Mechanisms A and B — see "Profiling `run!`" below. Two small,
  isolated candidates identified (~1.4 s combined of 28.8 s); precompilation explicitly
  ruled out of scope for now.
- 2026-08-21: simplified the Change 1 machinery after review feedback that it read as
  convoluted. Removed `SET_OPTIONS`, its sync `@assert`, the duplicated second
  `_set_options` method and the `_set_step` helper (net −19 lines) with no compile-time
  cost. See "Simplification of the Change 1 machinery" below.

## Problem description

`initialize!(model)` for a default `PrimitiveWetModel` takes **~34 s** in a fresh session,
of which `@time` attributes **99.2 % to compilation**. Actual work is ~0.2 s: a second
model built from the same types initializes in 0.8 s, and re-running step 3 on an already
warm session takes 0.024 s.

Measured breakdown (Julia 1.12.7, macOS aarch64, default `SpectralGrid()`, `NF=Float32`):

| Step | Time | Share |
|------|------|-------|
| 1. initialize components | 3.0 s | 9 % |
| 2. `Variables(model)` | 7.6 s | 22 % |
| 3. `initialize!(variables, model)` (initial conditions) | 22.7 s | 67 % |
| **total** | **34.0 s** | |

Within step 1, only three components are non-trivial: `orography` 2.15 s,
`land_sea_mask` 0.37 s, `land` 0.14 s. The rest are ≲0.1 s each.

Within step 3, the cost is spread evenly across the seven initial-condition calls:

| Substep | Time |
|---------|------|
| `initialize!(vars, model.ocean, model)` | 1.99 s |
| `initialize!(vars, model.land, model)` | 3.06 s |
| `ZonalWind` (vordiv) | 4.96 s |
| `PressureOnOrography` | 2.69 s |
| `JablonowskiTemperature` | 2.55 s |
| `ConstantRelativeHumidity` | 3.63 s |

That flat profile — every initial condition costing 2–5 s regardless of what it computes —
is the signature of a shared code path being re-inferred once per call site rather than of
any one expensive computation.

## Background

`@snoop_inference` (SnoopCompileCore) over `initialize!(model)` attributes inference time
by method. The top entries, aggregated over all specializations:

```
6.397 s   8 specializations   set!(kwargs, ::typeof(set!), vars::Variables, model::AbstractModel)   src/variables/set.jl:393
4.063 s   2                   initialize!(vars, initial_conditions::NamedTuple, model)              src/dynamics/initial_conditions.jl:58
1.741 s   7                   set!(group, kwargs, ::typeof(set!), variables::Variables, args...)     src/variables/set.jl:76
1.026 s  11                   set!(step, add, spectral_transform, ..., kwargs, vars::NamedTuple, geometry::Geometry)  src/variables/set.jl:29
0.984 s   1                   variables(model::PrimitiveWet)                                        src/models/primitive_wet.jl:144
0.699 s   7                   filter(f, t::Tuple)                                                   Base tuple.jl:526
0.594 s   2                   build_fuse_parents(all_vars, model)                                   src/variables/variables.jl:398
0.450 s   1                   variables(::Type{<:PrimitiveWet}, nsteps)
0.438 s   1                   variables(::Type{<:PrimitiveDry}, nsteps)
0.413 s  16                   grids_match(A::AbstractField, Bs::AbstractField...)
0.398 s   8                   filter_variables(vars, VariableType)                                  src/variables/variables.jl:613
```

By module: SpeedyWeather 23.5 s, Base 5.5 s, SpeedyTransforms 1.8 s,
KernelAbstractions 1.2 s, RingGrids 1.0 s.

Two independent mechanisms account for almost all of it.

### Mechanism A — `set!` re-specializes per keyword signature

`set!(vars, model; kwargs...)` at [set.jl:393](../../../SpeedyWeather/src/variables/set.jl)
is the single entry point every initial condition uses. Its body loops
`for varname in keys(kwargs)`, which Julia unrolls at compile time against the concrete
`Base.Pairs{...,NamedTuple{names,types}}` kwargs type. Every distinct *(set of keyword
names, tuple of value types)* therefore produces a fresh specialization that re-infers the
whole downstream chain: `set!#221` → `set!#220` → the per-variable `set!` methods →
`transform!` → the Legendre/FFT kernels → KernelAbstractions launch code.

There are exactly eight such call sites in `initial_conditions.jl`, which matches the
eight specializations SnoopCompile reports:

```
initial_conditions.jl:147   set!(vars, model; vorticity = ξ)
initial_conditions.jl:236   set!(vars, model; vorticity = ξks)
initial_conditions.jl:451   set!(vars, model; vorticity = vor_ic, divergence = div_ic, static_func = true)
initial_conditions.jl:653   set!(vars, model; temperature = temp_grid)
initial_conditions.jl:793   set!(vars, model; pressure = lnp_grid)
initial_conditions.jl:813   set!(vars, model; pressure = log(...))      # scalar, not a field
initial_conditions.jl:854   set!(vars, model; humidity = humid_grid)
initial_conditions.jl:926   set!(vars, model; η = η)
```

Confirmed directly — each *new* signature costs a fresh 2.4–3.5 s, a repeat is free:

```
set! temperature=field    3.53 s   (99.99% compilation)
set! humidity=field       2.59 s   (99.99% compilation)
set! pressure=field       2.80 s   (99.99% compilation)
set! pressure=scalar      2.36 s   (100.0% compilation)
set! temperature=field    0.00012 s          <-- repeat signature, free
```

Warming just those four signatures cut step 3 from 22.7 s to 15.4 s.

Note that `pressure=field` and `pressure=scalar` are separate specializations: the value
type is part of the kwargs type, so a scalar and a `Field` for the same keyword do not
share inference work.

### Mechanism B — variable definitions are carried in a 60-element heterogeneous `Tuple`

`variables(model)` returns a `Tuple` with 60 distinct element types. `filter_variables`
then calls `filter` and `unique` on it once per variable type (8 types), and
`build_fuse_parents` / `copyto!` iterate it as well. Inference over a 60-element
heterogeneous tuple type is superlinear and does not scale.

Measured on the real 60-element tuple vs. the same data as a `Vector{AbstractVariable}`:

| Operation | on `Tuple` | on `Vector` |
|-----------|-----------|-------------|
| `filter(x -> x isa PrognosticVariable, v)` | 0.387 s | 0.035 s |
| `unique(identifier, v)` | 0.455 s | 0.046 s |

Roughly **10×**, and these run 8× inside `filter_variables`.

### Not a mechanism: orography

`initialize!(model.orography, model)` (2.15 s) has no hotspot — snooping it shows a flat
tail of NetCDF/TOML reading, `grids_match`, and the transform stack, all first-use costs
shared with the rest of the model. There is no structural fix here; it is a
precompilation candidate only.

### Specialization is per-`NF`

Changing `NF` from `Float32` to `Float64` (and grid and `nlayers`) costs a further
**24.9 s** — almost the full price again. Changing only truncation/`nlayers` at fixed `NF`
is nearly free (0.5 s). This matters for the design: a `@compile_workload` only pays for
the `NF` it covers, so it should target `Float32` (the default), and the structural fixes
below are worth more than precompilation alone because they shrink the cost *per* `NF`.

There is currently **no `PrecompileTools` workload anywhere in the monorepo** —
`grep` for `@compile_workload` / `precompile(` returns nothing in all five packages.

## Summary of changes

Three changes, in decreasing benefit-to-risk order.

### Change 1 — type-erase the `set!` keyword loop (largest win, ~9 s) — **DONE**

Give `set!` a non-kwarg, non-specializing core that takes the variables to set as a
runtime `Vector{Pair{Symbol,Any}}` instead of a compile-time-unrolled `Pairs` type. The
public kwarg API stays exactly as it is; it just funnels into the shared core:

```julia
function set_core!(
        vars::NamedTuple, geometry, pairs::Vector{Pair{Symbol, Any}};
        step = 1, add::Bool = false, spectral_transform = nothing,
        coslat_scaling_included::Bool = false, static_func::Bool = true, namespace = nothing,
    )
    # ... same logic, but `for (varname, value) in pairs` is a runtime loop
end

set!(vars::Variables, model::AbstractModel; group::Symbol = :prognostic, kwargs...) =
    set_core!(
        getfield(vars, group), model.geometry,
        Pair{Symbol, Any}[k => v for (k, v) in kwargs];
        spectral_transform = model.spectral_transform,
    )
```

Prototyped and measured against the current implementation:

| signature | current | prototype |
|-----------|---------|-----------|
| `temperature = field` | 3.53 s | 1.95 s |
| `humidity = field` | 2.59 s | 0.79 s |
| `pressure = field` | 2.80 s | 1.07 s |
| `pressure = scalar` | 2.36 s | 0.85 s |
| `vorticity = field` | — | 0.79 s |
| **total (5 signatures)** | **14.4 s** (4 sigs) | **5.5 s** (5 sigs) |

**Implemented, with a correction to the design.** Erasing the variables into a
`Vector{Pair{Symbol,Any}}` was necessary but *not sufficient*. Two further steps were
needed, each found by re-measuring after the previous one failed to move the number:

1. **The options must be passed positionally, not as keywords.** The first attempt split
   the kwargs into variables (erased to a vector) and options (splatted back out as
   keywords). That left `set_variables!` with one specialization per distinct *option*
   NamedTuple type, and the per-signature cost was unchanged at ~2.4 s. Collecting the
   options into a concrete `SetOptions` struct passed positionally gives `set_variables!`
   exactly one signature (confirmed: `Base.specializations` went from one-per-call-site
   to 2).

2. **`@nospecialize` on the kwargs at the outer boundary.** Even with (1), the cost stayed
   at ~2.4 s. Snooping a single `set!` call showed 2.12 s sitting in the *kwargs wrapper of
   the entry point itself* (`#set!#264`), whose body is three lines: Julia was inferring
   `_set_pairs`/`_set_options` against the concrete kwargs type and constant-propagating
   through them, so the erasure happened too late to help. Adding `@nospecialize` to the
   `kwargs` of the three entry points and to the two helpers stops the specialization at
   the boundary.

Measured after all three steps (repeat signatures, i.e. the marginal cost of one more
distinct `set!` call site):

| signature | before | after |
|-----------|--------|-------|
| `temperature = field` (first call, compiles shared chain) | 3.53 s | 4.89 s |
| `humidity = field` | 2.59 s | **0.43 s** |
| `pressure = field` | 2.80 s | **0.73 s** |
| `pressure = scalar` | 2.36 s | **0.51 s** |
| `vorticity = field` | 2.25 s | **0.44 s** |

The first call is now *more* expensive because it compiles the whole shared chain once;
every subsequent distinct signature is ~5× cheaper. That is the intended trade: the cost is
paid once rather than once per call site.

Independent confirmation: the `set!` unit tests (`test/dynamics/set.jl`, 49 tests) dropped
from **2m29s to 49s**.

An `@nospecialize` note for future readers: it is load-bearing on the *entry points*, not
just the helpers. Removing it from `set!(vars, model; kwargs...)` alone restores most of the
original cost.

Trade-off: the inner loop becomes dynamically dispatched. This is irrelevant here —
`set!` is called ~10 times during setup and never in the time-stepping loop. Runtime after
warmup was 0.024 s for all of step 3; the dispatch overhead is far below that.

#### Simplification of the Change 1 machinery

The first implementation was reviewed as convoluted. Three of its four moving parts turned
out to be incidental rather than load-bearing, and were removed on 2026-08-21 with **no
measurable compile-time cost** (net −19 lines):

- **`SET_OPTIONS` and its `@assert` — deleted.** The const duplicated
  `fieldnames(SetOptions)` and needed an assertion to keep the two in sync. Replaced by
  `is_option(name) = name in fieldnames(SetOptions)`. This removes the hazard recorded in
  Known limitations below: a new option is now added in **one** place, the struct.
- **The two `_set_options` methods — collapsed to one.** They differed only in whether
  `spectral_transform` was taken from the model. Making `SetOptions` a `Base.@kwdef` struct
  moves the defaults onto the struct, and a `defaults...` keyword handles the override:

  ```julia
  function _set_options(@nospecialize(kwargs); defaults...)
      options = merge(NamedTuple(defaults), NamedTuple(k => v for (k, v) in kwargs if is_option(k)))
      return SetOptions(; options...)
  end
  ```

  Twelve lines of hand-written `get(kwargs, :field, default)` became four.
- **`_set_step` — deleted.** `@kwdef`'s implicit `convert` to the declared
  `Union{Nothing, Int}` field type does the same job. Verified equivalent for `2`, `2.0`,
  `UInt8(2)`, `Int32(2)` and `nothing`.
- **The u/v special case** dropped `varnames = map(first, ...)`, two `findfirst` closures
  and two `has_*` flags in favour of one named `value_or_nothing` helper.

**What was deliberately kept:** the `@nospecialize` on the entry points and the positional
`SetOptions` struct. These are the two non-obvious pieces — removing either restores most of
the original cost — so the comment explaining them stays in the source. They are the part a
future reader is most likely to "simplify" away and silently lose ~9 s.

Verification after simplification:

| check | result |
|-------|--------|
| `test/dynamics/set.jl` | 49/49 pass, 48.7 s (vs 49 s before, 2m29s pre-Change-1) |
| `test/dynamics/initial_conditions.jl` | all pass |
| `Base.specializations(set_variables!)` | **3** — erasure intact, not one per call site |
| marginal cost per new `set!` signature | 0.40–0.49 s (was 0.43–0.73 s) |
| `initialize!` total | 22.0–22.2 s vs 21.9 s baseline — within noise |
| Runic | clean |

### Change 2 — collect variable definitions into a `Vector`, not a `Tuple` (~3–4 s) — **DONE**

In [variables.jl](../../../SpeedyWeather/src/variables/variables.jl), have
`variables(model)` / `all_variables` produce a `Vector{AbstractVariable}` and make
`filter_variables`, `build_fuse_parents` and `_allocate_namespace` operate on it. The
element type is already effectively erased downstream (`_allocate_namespace` builds a
`Pair{Symbol,Any}[]`, `build_fuse_parents` uses `Dict{...,Vector{AbstractVariable}}`), so
this mostly removes a type-level detour that buys nothing.

Note `variables(nt::NamedTuple, model)` at variables.jl:~605 ends in `|> Tuple` — that is
the place the tuple is formed.

The allocated `Variables` container itself stays a concretely typed `NamedTuple`; only the
*definition* objects flowing through the collection pipeline become a vector. This does
not affect runtime access to `vars.prognostic.vorticity` etc.

**Implemented.** Two edits in `src/variables/variables.jl`:

- `all_variables(model)` now `append!`s into an `AbstractVariable[]` instead of growing a
  tuple type via `t = tuple(t..., vars...)` in a loop over ~30 components. That loop was
  the main offender: it re-inferred an ever-larger heterogeneous tuple type at each of the
  ~30 iterations.
- `filter_variables` keeps vectors throughout, and `group` values become
  `Vector{AbstractVariable}` instead of `Tuple`. `get_namespaces` correspondingly takes a
  collection rather than varargs (the varargs splat `get_namespaces(vars...)` was itself a
  specialization hazard on a 60-element collection).

The component-facing extension API is untouched: components still define
`variables(::Component, ::Model)` returning tuples, and `append!` accepts those directly.

Measured: `Variables(model)` **7.62 s → 4.63 s** (−3.0 s, −39 %); total `initialize!`
**34.0 s → 31.4 s**. Steps 1 and 3 unchanged, as expected.

Verified bitwise-identical: all 113 allocated arrays across all seven variable groups
match the pre-change state exactly (size, sum of absolute values, and leading elements).

Snooping `Variables(model)` afterwards confirms `filter`, `unique`, `build_fuse_parents`
and `filter_variables` have dropped off the profile entirely (together ~1.7 s before).
The remaining ~1.9 s is in the `variables(::Type{<:PrimitiveWet}, nsteps)` /
`variables(model::PrimitiveWet)` methods *constructing* their large literal tuples of
variable definitions — a separate concern from collecting them, and a candidate for
future work rather than part of this change.

### Change 3 — add a `PrecompileTools` workload (absorbs most of the remainder)

With Changes 1 and 2 in, add to `SpeedyWeather.jl`:

```julia
using PrecompileTools
@setup_workload begin
    @compile_workload begin
        spectral_grid = SpectralGrid(truncation = 31, nlayers = 4)   # NF = Float32
        model = PrimitiveWetModel(spectral_grid)
        simulation = initialize!(model)
    end
end
```

Use a small truncation and few layers — the compiled code is resolution-independent
(a different truncation at the same `NF` costs 0.5 s), so the workload should be as cheap
to *run* as possible while covering the same methods.

Caveats to weigh before adopting:

- It only covers `NF = Float32`. `Float64` users pay ~25 s regardless. Covering both
  roughly doubles precompile time and package image size.
- It shifts cost from every session start to every `Pkg.precompile`. The package already
  takes ~22 s to precompile; this would add meaningfully to that.
- Orography initialization reads NetCDF assets from disk. Running it inside a workload
  either requires those assets at precompile time or means excluding that path. Prefer
  a `NoOrography`/`ZonalRidge` variant in the workload, or accept that the NetCDF read
  stays uncached, to avoid a precompile-time dependency on downloadable assets.

Recommendation: land Changes 1 and 2 first and re-measure. They are pure wins with no
precompile-time or image-size cost. Decide on Change 3 afterwards against the reduced
baseline.

### Result (measured)

| | before | after 2 | after 1 + 2 |
|---|---|---|---|
| 1. components | 3.0 s | 3.4 s | 3.6 s |
| 2. `Variables(model)` | 7.6 s | 4.6 s | 4.5 s |
| 3. initial conditions | 22.7 s | 23.1 s | **13.8 s** |
| **total `initialize!`** | **34.0 s** | 31.4 s | **22.6 s** |

Close to the 20–22 s estimate. All 113 allocated variable arrays remain bitwise identical
to the pre-change baseline after both changes.

Change 3 on top could take the default-`NF` path to a few seconds, at the cost of a
longer `Pkg.precompile` and a larger image.

**Measurement caveat:** compile-time numbers on this machine are very sensitive to CPU
contention — an intermediate measurement taken while the test suite was running showed
every step roughly doubled, including steps not touched by the change. Always re-measure on
an idle machine before concluding anything.

## Profiling `run!(simulation)`

Profiled 2026-08-21 on the branch state after Changes 1 and 2 (i.e. with `initialize!`
already reduced to ~22 s). Investigation only — **no changes were made**, and
precompilation was explicitly ruled out of scope.

### Baseline

Default `PrimitiveWetModel`, `truncation = 32`, `nlayers = 8`, `NF = Float32`,
Julia 1.12.7, macOS aarch64, fresh session:

| Step | Time | Compilation share |
|------|------|-------------------|
| `initialize!(model)` | 21.9 s | 98.8 % |
| **first `run!(simulation, period = Day(1))`** | **28.8 s** | **99.4 %** |
| second `run!` (same session) | 0.15 s | — |

So after Changes 1 and 2, `run!` costs *more* than `initialize!` and is now the larger
half of a cold-start simulation. Actual work in `run!` is 0.15 s; the rest is compilation.

Note `run!` defaults to `output = false`, so the 28.8 s is the core path only. NetCDF
output is a separate cost — see below.

### Structural breakdown

`@snoop_inference` over the first `run!`. Absolute numbers are inflated roughly 2.5× by
snooping overhead (75.9 s attributed vs 28.8 s wall clock), so read them as proportions:

```
kwcall(run!)                                     75.9 s
└─ #run!#1224                                    71.3 s
   ├─ time_stepping!                             44.4 s
   │  └─ time_step!(simulation, time_stepping)   39.9 s
   │     ├─ time_step!(vars, ts, ::PrimitiveEquation)  30.8 s
   │     │  ├─ dynamics_tendencies!              17.4 s
   │     │  ├─ reset_tendencies!                  3.1 s
   │     │  ├─ parameterization_tendencies!       2.1 s
   │     │  ├─ update_prognostic!                 1.4 s
   │     │  ├─ diffusion_and_implicit!            1.3 s
   │     │  └─ land_timestep!                     1.1 s
   │     ├─ progress!(feedback, vars, model)      3.0 s
   │     └─ reinitialize!(::PrimitiveWetModel)    2.0 s
   ├─ initialize!(simulation)                    20.2 s
   └─ finalize!(simulation)                       2.3 s
```

By module: SpeedyWeather 18.8 s, SpeedyTransforms 2.6 s, KernelAbstractions 1.6 s,
Base 1.6 s, Base.Broadcast 1.4 s.

### Conclusion: no dominant mechanism

Unlike `initialize!`, the `run!` profile is **flat**. There is no `set!`-style
keyword-signature explosion and no oversized heterogeneous tuple. Specifically:

- `dynamics_tendencies!` ([tendencies.jl:119](../../../SpeedyWeather/src/dynamics/tendencies.jl))
  is the single largest entry at 17.4 s inclusive, of which **1.39 s is exclusive
  self-inference of one method body**. That body is a long straight-line sequence of ~15
  calls (pressure gradient, geopotential, vertical velocity/advection, grid tendencies,
  batched transform, spectral tendencies, tracer advection). The cost is the size of the
  dynamical core itself, not a structural defect. No refactor identified.
- The column parameterizations are entirely flat — no individual `parameterization!`
  method exceeds 0.1 s (largest: `BulkRichardsonDiffusion` 0.09 s, `OneBandLongwave`
  0.06 s). Nothing to target.
- `reset_tendencies!` (3.1 s) is a deliberate recursive tuple unroll, commented in
  [tendencies.jl:1269](../../../SpeedyWeather/src/dynamics/tendencies.jl) as avoiding
  Union-typed iteration that Enzyme cannot differentiate. Converting it to a runtime loop
  the way Change 2 did would trade compile time against differentiability — not a free win,
  and out of scope here.
- The `Base.Iterators.Filter` entry (0.89 s over 8 specializations) that initially looked
  like a Mechanism-B repeat turned out to be diffuse Base machinery reached via
  `broadcast_unalias` → `unaliascopy(::LowerTriangularArray)`, not a SpeedyWeather
  collection pattern.

### Two small candidates (not implemented)

Both are real and independently measured, but together worth only ~1.4 s of 28.8 s.

**1. `ProgressMeter.next!` on a disabled meter — 0.40 s.**
`Feedback.verbose` defaults to `isinteractive()`, so in scripts and CI the progress meter
is disabled — yet `progress!(feedback::Feedback)` at
[feedback.jl:91](../../../SpeedyWeather/src/output/feedback.jl) still calls
`ProgressMeter.next!`, whose `enabled` check happens at *runtime*, inside the method. The
kwarg method is compiled purely to bump a counter nobody displays.

Measured in isolation: `ProgressMeter.next!` on a disabled `Progress` costs 0.40 s,
100 % compilation. An A/B of the guarded version over the full `run!`:

| | first `run!` |
|---|---|
| baseline | 29.19 s, 28.76 s |
| with guard | 28.74 s, 28.64 s |

i.e. ~0.3 s, consistent but marginal.

**Correctness trap:** the obvious guard (`enabled || return nothing`) is **wrong**.
`progress_meter.counter` is read by `progress!(feedback, vars, model)` to schedule
`nan_detection!` and debug output ([feedback.jl:99–132](../../../SpeedyWeather/src/output/feedback.jl)),
so it must still be incremented when the meter is off. Any fix has to bump the counter
directly and skip only the display work.

**2. `initialize!(implicit, ...)` — 1.05 s.**
Reached via `reinitialize!(::PrimitiveWetModel, vars)` → `reinitialize!(::ImplicitPrimitiveEquation, ...)`
on the first timestep. This is *not* redundant with `initialize!(model)`: the implicit
solver is built lazily and guarded by `implicit.Δt[] == Δt`, so `initialize!(model)` never
compiles it. It is genuine first-use cost with no obvious structural fix.

### The NetCDF output path is a separate, larger cost

`run!(..., output = true)` costs an **additional ~10.8 s** on top of the core path, and it
is fully disjoint: a session that has already run `run!` without output still pays the
full amount on the first output run.

The top exclusive entry in the entire `run!` profile belongs to this path:
`initialize!(::NetCDFOutput, vars, model)` at
[netcdf_output.jl:147](../../../SpeedyWeather/src/output/writers/netcdf_output.jl),
**2.46 s of exclusive self-inference** — again a long linear body, here of `defDim`/`defVar`
calls into NCDatasets. Same shape as `dynamics_tendencies!`: big body, no hotspot.

This is worth recording because any future precompilation work scoped to `initialize!` and
a plain `run!` would leave this 10.8 s entirely uncovered for every user who writes output.

## Testing and verification

- `set!` is widely used and user-facing; the existing `set!` tests must pass unchanged.
  Its keyword API is deliberately not altered by Change 1.
- Add a test that `set!` with a scalar, a `Field`, a `LowerTriangularArray` and a function
  all still dispatch correctly for both a namespaced (`:ocean`) and a non-namespaced
  variable, and that the `u`/`v` → vorticity/divergence special case still fires.
- Verify unchanged initial state: run `initialize!` before and after, and compare
  `vars.prognostic` arrays bitwise. The changes are meant to be purely a compile-time
  refactor, so this should be exact.
- Re-run the profiling scripts to confirm the numbers, and check
  `SpeedyWeather/benchmark/` for time-stepping regressions (none expected — `set!` and the
  variable-collection pipeline are setup-only paths).
- Test with `--check-bounds=yes` per project convention.

## Documentation changes

None required for Changes 1 and 2 — the public API is unchanged. If Change 3 lands, note
the precompilation behaviour and the `NF=Float32` coverage in the docs.

## Known limitations

- Compile time remains proportional to the number of distinct `NF` values used in a
  session; none of these changes make specializations shared across number formats.
- Change 1 trades a small amount of runtime dispatch for compile time. Correct here
  because `set!` is a setup-path function, but the same trade would be wrong inside the
  time-stepping loop.
- The ~0.4–0.7 s per-signature residual after Change 1 is inherent to first-inference of the
  transform stack and can only be removed by precompilation.
- After Changes 1 and 2, `run!` (28.8 s) costs more than `initialize!` (21.9 s) in a cold
  session, and no structural fix was found for it — the cost is spread across the
  dynamical core rather than concentrated in a re-specialization pattern.
- `SetOptions` fixes the set of `set!` options at the struct definition: a keyword that is
  not one of its fields is treated as a variable to set, so a typo'd option produces a
  "not defined" warning rather than an error. Adding a new option is a one-line change to
  the struct since the 2026-08-21 simplification.

## Future work

- The `Float64` path costs the same ~25 s and is untouched by a `Float32` workload.
- `initialize!(model.orography, model)` (2.15 s) is a diffuse first-use cost across the
  NetCDF/TOML/transform stack; worth revisiting if it dominates after these changes.
- `grids_match` shows 16 specializations for 0.41 s — a small but easy candidate for
  `@nospecialize` on its varargs.
- Guard `progress!(feedback::Feedback)` so a disabled progress meter does not compile
  `ProgressMeter.next!` (~0.3–0.4 s). Must keep incrementing `progress_meter.counter` —
  see the correctness trap noted under "Profiling `run!`".
- The NetCDF output path costs ~10.8 s on first use, disjoint from the core `run!` path,
  with 2.46 s of exclusive self-inference in `initialize!(::NetCDFOutput, ...)`. Untouched
  by anything in this document and the largest single remaining block.
- `reset_tendencies!` (3.1 s) could in principle move from a recursive tuple unroll to a
  runtime loop as Change 2 did, but the unroll is deliberate — it avoids Union-typed
  iteration that Enzyme cannot differentiate. Would need a differentiability check first.
