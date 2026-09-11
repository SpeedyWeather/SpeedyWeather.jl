# One `radiation` model component bundling shortwave and longwave

> Status: **in progress**. Implemented on `mg/numericalradiation`; unit tests and the whole-model
> bitwise comparison against the base revision pass on CPU. Awaiting human review/sign-off, GPU test
> run, PR number for the changelog, and the benchmark check.

Date of initial draft: 2026-09-11

Base revision: `1ef2a0e2` (`mg/numericalradiation`, on top of v0.22.1)

## Originating prompt

> I would like to run the ecCKD model that is defined in this package [NumericalRadiation.jl]
> with SpeedyWeather. Scan the package and SpeedyWeather.jl and make a plan how to do that.

> We might have to adjust SpeedyWeather upstream. This is fine and should be included in the plan.

> What do you think about in SpeedyWeather just having one radiation, `Radiation` field that
> bundles both shortwave and longwave. I am confident we can also make this work with the
> existing ones first, this would also need to be proven and implemented first upstream. In this
> case we want absolute identical results, but with the flexibility to add NumericalRadiation
> later.

> start with the SpeedyWeather upstream work in
> /Users/max/Nextcloud/SpeedyWeather/alt-version/SpeedyWeather.jl, also write a plan in this repo
> in accordance with its CLAUDE.md

## Revision log

- **2026-09-11, tests run.** `test/parameterizations/{radiation,longwave_radiation,shortwave_radiation,stochastic_physics}.jl`
  with `--check-bounds=yes --depwarn=yes`: bit-identity holds for the full 4×5 scheme matrix in
  Float32 and the default pair in Float64. One assertion of mine was wrong, not the code: a
  transparent shortwave with `longwave = nothing` correctly leaves the temperature tendency at
  zero; the "non-trivial tendency" check now only applies when a longwave scheme is present.
- **2026-09-11, "just keep those definitions in variables(::AbstractLongwave) and
  variables(::AbstractShortwave)".** Removed the `shortwave_variables()`/`longwave_variables()`
  helpers and the `variables(::AbstractRadiation)` union default again; the corresponding test
  case was dropped.
- **2026-09-11, whole-model comparison done.** `PrimitiveWetModel` and `PrimitiveDryModel`, T31 L8,
  default initial conditions, 2 days from 2000-01-01, Float32, CPU: all 19 dumped arrays (vorticity,
  divergence, temperature, pressure, humidity, grid temperature tendency, `outgoing_longwave`,
  `outgoing_shortwave`, `surface_longwave_down`, `surface_shortwave_down`) are `==` between a
  pristine worktree of `1ef2a0e2` and this branch.
- **2026-09-11, during implementation.** Dropped the `GreyRadiation` helper from the draft; the dry
  model spells out its default `Radiation(...; shortwave = OneBandGreyShortwave, longwave = OneBandGreyLongwave)` inline, like the other dry defaults.
- **2026-09-11, initial draft.** Written together with the companion plans in
  NumericalRadiation.jl (`docs/plans/ecckd_speedyweather.md`, `docs/plans/speedyweather_upstream.md`,
  item U1). Implementation started immediately after drafting.

## Problem description

`PrimitiveWetModel` and `PrimitiveDryModel` carry two independent radiation components,
`shortwave_radiation <: AbstractShortwave` and `longwave_radiation <: AbstractLongwave`, called
back to back inside the fused column kernel (`column_parameterizations!`). A correlated-k scheme
such as ecCKD (NumericalRadiation.jl) computes gas optical properties once and then solves both
streams. With two separate components that work is either done twice, or the shortwave component
has to smuggle its optics to the longwave one through scratch arrays plus an assumption on the
call order. Neither is acceptable for a scheme that will dominate the cost of the physics.

The requirement for this change is that it is **bit-identical** with the existing schemes: a
model built with `Radiation(shortwave, longwave)` must produce exactly the same tendencies and
diagnostics as one with the two separate components. Only after that is proven does an external
radiation scheme plug in as a single `radiation = ...` component.

## Background

- `AbstractRadiation <: AbstractParameterization` already exists as the common supertype of
  `AbstractShortwave` and `AbstractLongwave` (`parameterizations/radiation/shortwave_radiation.jl`).
- Composite parameterizations with `@component` sub-fields are an established pattern:
  `OneBandLongwave{T, R}` (transmissivity + radiative transfer), `OneBandShortwave{C, T, R}`,
  `OceanLandAlbedo{Ocean, Land}`. `OceanLandAlbedo` also shows how a composite collects the
  variables of its sub-components in `variables`.
- `parameterization!(ij, vars, p, model)` is dispatched per column in
  `column_parameterizations!`; the order is the `parameterizations` tuple of the model, which
  today lists `:shortwave_radiation, :longwave_radiation` consecutively. Wrapping the two calls
  in one function in the same order changes nothing about the floating-point operations, hence
  identical results are expected by construction.
- Duplicate `ParameterizationVariable`s are removed by identifier in `filter_variables`, so a
  composite can concatenate the variables of its parts without de-duplicating itself.
- `variables(component, model)` (two-argument form) falls back to `variables(component)`.
- Fallbacks for `nothing` exist for `parameterization!`, `initialize!`, `variables`, so either
  stream can stay `nothing`.
- No output writer or surface model reads `model.shortwave_radiation` / `model.longwave_radiation`;
  they read `vars.parameterizations.*`, which does not change. References to the two fields in
  the repository (at the base revision) are confined to `src/models/primitive_{wet,dry}.jl`,
  three test files (`test/parameterizations/{longwave,shortwave}_radiation.jl`,
  `stochastic_physics.jl`) and eight lines in `docs/src/radiation.md` and
  `docs/src/parameterizations.md`.
- GPU: `Adapt.adapt_structure(to, model)` only adapts `core_components`; the parameterizations
  reach the kernel through `get_parameterizations(model)`, which is `@generated` from the
  `parameterizations` tuple. A `Radiation` struct with `Adapt.@adapt_structure` adapts like
  `OneBandLongwave` does.

## Summary of changes

### New type `Radiation` (`src/parameterizations/radiation/radiation.jl`, new file)

```julia
export Radiation

@parameterized @kwdef struct Radiation{SW, LW} <: AbstractRadiation
    @component shortwave::SW
    @component longwave::LW
end

Radiation(SG::SpectralGrid; shortwave = OneBandShortwave(SG), longwave = OneBandLongwave(SG)) =
    Radiation(shortwave, longwave)

Adapt.@adapt_structure Radiation
Base.show(io::IO, R::Radiation) = show(io, R, values = false)

initialize!(radiation::Radiation, model::PrimitiveEquation)   # initialize! shortwave then longwave
variables(radiation::Radiation, model::AbstractModel)          # concatenation of both sub-schemes' variables

@propagate_inbounds function parameterization!(ij, vars, radiation::Radiation, model)
    parameterization!(ij, vars, radiation.shortwave, model)    # same order as today
    parameterization!(ij, vars, radiation.longwave, model)
    return nothing
end
```

### Diagnostics of third-party radiation schemes

`variables(::AbstractShortwave)` and `variables(::AbstractLongwave)` stay as they are. A scheme
that subtypes `AbstractRadiation` directly (one component computing both streams) declares its
own variables, e.g. by concatenating `variables(OneBandShortwave(SG))` and
`variables(OneBandLongwave(SG))`, so that the ocean, land and output code find the standard
fields. (An earlier draft factored the two method bodies into helpers with a union default on
`AbstractRadiation`; dropped at the user's request, see revision log.)

### Model structs (`src/models/primitive_wet.jl`, `src/models/primitive_dry.jl`)

- Type parameters `SW, LW` → `RA`.
- Fields `shortwave_radiation`, `longwave_radiation` → `@component radiation::RA = Radiation(spectral_grid)`
  (dry model: `Radiation(spectral_grid; shortwave = OneBandGreyShortwave(spectral_grid), longwave = OneBandGreyLongwave(spectral_grid))`).
- `parameterizations` tuple: `:shortwave_radiation, :longwave_radiation` → `:radiation`.
- `initialize!`: two calls → one.
- Backwards compatibility: the positional constructors `PrimitiveWetModel(spectral_grid; kwargs...)`
  and `PrimitiveDryModel(spectral_grid; kwargs...)` accept `shortwave_radiation` and/or
  `longwave_radiation`, wrap them into `Radiation` (missing stream from the default), and emit a
  deprecation warning once. A user-supplied `parameterizations` tuple containing
  `:shortwave_radiation` or `:longwave_radiation` has those entries replaced by a single
  `:radiation` at the position of the first, with a deprecation warning. Reading
  `model.longwave_radiation` is not shimmed; it errors with the usual `FieldError`, whose
  message is enough to find `model.radiation.longwave`.

### Tests

- `test/parameterizations/radiation.jl` (new): the bit-identity test, see below.
- `test/parameterizations/longwave_radiation.jl`, `shortwave_radiation.jl`,
  `stochastic_physics.jl`: construct with `radiation = Radiation(spectral_grid; ...)` and read
  `model.radiation.longwave` / `.shortwave`. One case each keeps the deprecated keyword to test
  the shim.

### Versioning and changelog

Breaking change of the model struct → `SpeedyWeather/Project.toml` version `0.22.1+DEV` →
`0.23.0-DEV`. CHANGELOG entry as first bullet under `## Unreleased`.

## Testing and verification

1. **Bit-identity unit test** (`test/parameterizations/radiation.jl`). For every combination of
   shortwave `{TransparentShortwave, OneBandShortwave, OneBandGreyShortwave, nothing}` and
   longwave `{UniformCooling, JeevanjeeRadiation, OneBandLongwave, OneBandGreyLongwave, nothing}`
   at `truncation = 32, nlayers = 8`: build the model with `Radiation(shortwave, longwave)`,
   allocate `Variables`, set a non-trivial state (temperature, humidity, pressures, albedos,
   `cos_zenith` via the zenith parameterization), `deepcopy` it, run
   `parameterization!(ij, vars_a, model.radiation.shortwave, model)` followed by the longwave
   call on one copy and `parameterization!(ij, vars_b, model.radiation, model)` on the other for
   all `ij`, and assert `==` (not `≈`) on the temperature tendency and every
   `vars.parameterizations` field (including `ocean` and `land` namespaces). Full matrix in
   `Float32`, the default pair also in `Float64`. Also checks that the variables of `Radiation`
   are the union of the two stream-specific sets.
2. **Whole-model bitwise comparison against the base revision** (one-off, not kept in the test
   suite). A script runs `PrimitiveWetModel` and `PrimitiveDryModel` at T31 L8 for two days
   from the default (deterministic) initial conditions on a pristine worktree of `1ef2a0e2`
   and on the branch, serialising prognostic arrays and radiation diagnostics; the two dumps
   must be `==`. Result recorded in the revision log.
3. **Existing tests** in `test/parameterizations/` run with `--check-bounds=yes`; `all_parameterizations`
   and the model-construction tests cover the struct change.
4. **GPU**: the change is dispatch-only; `test/GPU/` should be run by whoever has a device before
   merging (not available in the drafting environment).
5. **Performance**: `parameterization_tendencies!` benchmark before/after on CPU; expected within
   noise since the extra call level is `@propagate_inbounds` and inlines.

## Documentation changes

- `docs/src/radiation.md`: model construction examples use
  `radiation = Radiation(spectral_grid; longwave = ...)` / `shortwave = ...`; a short new section
  "Radiation as one component" explaining `Radiation`, that either stream may be `nothing`, and
  that a scheme subtyping `AbstractRadiation` can implement both streams in one
  `parameterization!`.
- `docs/src/parameterizations.md`: the custom-parameterization example's `parameterizations`
  tuple uses `:radiation`; the sentence about call order refers to `:radiation`.
- Docstrings of `PrimitiveWetModel` / `PrimitiveDryModel` update automatically via `TYPEDFIELDS`.

## Known limitations

- `model.longwave_radiation` / `model.shortwave_radiation` field access is not shimmed; only the
  constructor keywords and the `parameterizations` tuple are.
- `Radiation` calls shortwave then longwave unconditionally; a scheme that needs a different
  order or interleaving should subtype `AbstractRadiation` directly.
- Sub-stepping radiation (calling it every N steps) is out of scope here; see the companion
  upstream plan (item U3) in NumericalRadiation.jl.

## Future work

- NumericalRadiation.jl defines `EcCKDRadiation <: SpeedyWeather.AbstractRadiation` and is
  passed as `radiation = EcCKDRadiation(...)`; mixed setups such as
  `Radiation(spectral_grid; longwave = EcCKDLongwave(...))` remain possible.
- Prescribed ozone, radiation call frequency and a `(npoints, n, nlayers)` variable dimension are
  tracked as U2–U4 in NumericalRadiation.jl's `docs/plans/speedyweather_upstream.md`.
