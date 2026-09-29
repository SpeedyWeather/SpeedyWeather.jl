# A NumericalRadiation extension of SpeedyWeather

> Status: **implemented** on `mg/numericalradiation-extension` (decisions in the revision log);
> the extension tests are part of the normal test suite. This document stands on its own: the
> NumericalRadiation side is summarised below and linked by PR, not by its plan files.

Date of initial draft: 2026-09-22

Base revision: `e524ff10` (`mg/numericalradiation`, on top of v0.22.1; carries the `Radiation`
bundle of `radiation-bundle.md`). Working branch: `mg/numericalradiation-extension`.

## Originating prompt

> I want to move the extension upstream to SpeedyWeather. For this purpose we need to delete
> the SpeedyWeather extension here in NumericalRadiation and at the same time create one in
> SpeedyWeather. For SpeedyWeather use `/Users/max/Nextcloud/SpeedyWeather/alt-version/SpeedyWeather.jl`
> as a working directory. Base your SpeedyWeather branch on `mg/numericalradiation` and call it
> `mg/numericalradiation-extension`. First make a new plan for this. Then let me review it.
> Then make the changes but don't commit yet.

## Revision log

- **2026-09-29, NumericalRadiation registered.** Version 0.1.0 is in General, but it is
  NumericalRadiation's `main` before its PR #16 and has neither `ClearSkyEcCKDRadiation` nor the
  per-array types the extension needs, so with 0.1.0 the extension would fail to load. The git
  source in `test/Project.toml` and the `Pkg.add` line of the 1.10 CI step stay until PR #16 is
  released as 0.1.1; `[compat] NumericalRadiation = "0.1.1"` in `Project.toml` excludes 0.1.0
  already, and the release then makes the git source unnecessary.
- **2026-09-29, tests in the normal suite.** The extension test is
  `test/parameterizations/numericalradiation.jl`, run by `Pkg.test` on Julia 1.10 and 1.13
  like every other test: NumericalRadiation is a test dependency (`test/Project.toml`, git
  source until it is registered; the 1.10 CI step adds it with `Pkg.add`). The dedicated
  environment and workflow are gone. This plan no longer refers to NumericalRadiation's plan
  documents.
- **2026-09-24, NumericalRadiation side collapsed.** The three-PR stack there became one PR
  (NumericalEarth/NumericalRadiation.jl#16, branch `mg/speedy-update`, merged with its `main`);
  the test environment pulls that branch (`main` once the PR is merged, the registered package
  once it exists).
- **2026-09-22, naming.** NumericalRadiation's scheme type is `ClearSkyEcCKDRadiation`
  (the solver pair is part of its definition); SpeedyWeather-side mentions renamed.
- **2026-09-22, third review.** Same for `EcCKDRadiation`: the configured clear-sky ecCKD scheme
  (gas optics, prescribed mole fractions with `o3` and `co2` defaults, surface emissivity)
  becomes a host-neutral type of NumericalRadiation (`src/ecckd_radiation.jl` there); this
  extension adds methods and constructors from a `SpectralGrid` only. SpeedyWeather's `src`
  is untouched apart from `Project.toml`: no exported type, no stub, no `src` file. The
  ocean/land emissivity pair collapses to one `surface_emissivity`; without a CO₂ component
  the scheme's prescribed `mole_fractions.co2` (280 ppm by default) is used.
- **2026-09-22, second review.** No wrapper type for the longwave: nothing in SpeedyWeather
  requires a longwave to subtype `AbstractLongwave` (the `Radiation` bundle has unconstrained
  type parameters), so the extension adds `variables`, `initialize!` and `parameterization!`
  methods to NumericalRadiation's `AnalyticBandLongwave` itself plus a constructor from a
  `SpectralGrid`; the time-step selection, which does need a model component, is defined by
  the extension for the scheme types; `WilliamsLongwave` dropped. The wrapper's only state, the CO₂ default, becomes
  a fallback to 280 ppm (NumericalRadiation's own `AtmosphereProfile` default) when the model
  has no `greenhouse_gases.co2`. `EcCKDRadiation` stays a SpeedyWeather-side type because it
  bundles the gas optics with configuration. The #1261 follow-up is no longer needed for this.
- **2026-09-22, review.** Decisions: (1) `WilliamsLongwave` for now, to be revisited with the
  `Parameterization` stub of [#1261](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/1261),
  which would let the extension wrap NumericalRadiation's own scheme types without new names
  (TODO noted in `src/parameterizations/radiation/numericalradiation.jl`); (2) examples and
  validation scripts stay in NumericalRadiation, together with the full unit tests of the
  coupling; SpeedyWeather gets only basic functionality tests; (3) a documentation section
  that states which schemes are built in and which come from NumericalRadiation, with a usage
  example (non-executed code block; an executed example needs NumericalRadiation in the docs
  environment and the ecCKD tables); (4) on the NumericalRadiation side the removal is squashed
  into the bottom PR of its stack (#16) and propagated upward; (5) NumericalRadiation will be
  registered within days by its maintainers, after which the dedicated test environment and
  workflow collapse into `test/Project.toml` (done 2026-09-29, see above).
- **2026-09-22, initial draft.**

## Problem description

The coupling of NumericalRadiation.jl's radiation schemes to SpeedyWeather currently lives in
NumericalRadiation as the package extension `NumericalRadiationSpeedyWeatherExt`
(branch `mg/ecckd-speedy`, PR #19 there): a `SpeedyAnalyticBandLongwave <: AbstractLongwave`
wrapping the analytic-band longwave scheme, and `EcCKDRadiation <: AbstractRadiation`, a
clear-sky ecCKD scheme computing both streams from one gas-optics evaluation. It is to move
here, into a SpeedyWeather extension triggered by loading NumericalRadiation, and be deleted
over there.

Reasons: the extension is SpeedyWeather-shaped code (parameterization variables, column
kernel, `Radiation` bundle, host constants), it changes whenever SpeedyWeather's
parameterization interface changes, and SpeedyWeather already hosts its coupling code to
other packages (Terrarium) the same way. NumericalRadiation stays a host-neutral column
radiation library and keeps only the small package-side changes made for the coupling
(per-array types of `ColumnAtmosphere`, an optional caller-owned scratch of the clear-sky
shortwave solver, a two-stream robustness fix), which any host benefits from.

## Background

- **Extensions here.** `SpeedyWeather/ext/` holds eight extensions, two as directories
  (`SpeedyWeatherTerrariumExt/`, `SpeedyWeatherZarrExt/`). Optional dependencies sit in
  `[weakdeps]` with a `[compat]` entry; the docs environment lists the ones it needs.
- **Exposing extension functionality.** Two patterns exist. `TerrariumOutput` is a function
  declared in `src/` (`function TerrariumOutput end`, exported) with its methods in the
  extension; `TerrariumLand` is a type defined inside the extension and reached through
  `Base.get_extension`. Radiation schemes are model components that users pass by name and
  test with `isa`, so this plan puts the *types* in `src/` and only the methods that need
  NumericalRadiation in the extension.
- **Tests.** `test/runtests.jl` autodiscovers `test/**/*.jl` and removes the
  `GPU/`, `differentiability/` and `reactant/` prefixes, which have their own environments
  and workflows (`CI_Enzyme.yml` develops the monorepo packages into
  `test/differentiability/` and runs a script there).
- **NumericalRadiation is not registered yet.** It can be a weak dependency (only a UUID is
  needed) but any environment that installs it needs a `[sources]` git entry, which requires
  Julia ≥ 1.11, or `Pkg.add(url = ...)`. `CI_SpeedyWeather.yml` runs Julia 1.10 and 1.13 and
  has a 1.10-only step that `Pkg.develop`s the monorepo packages because `[sources]` is
  ignored there; the same step adds NumericalRadiation from git (checked: `Pkg.add` of a weak
  dependency from a git URL works on 1.10, moves it to `[deps]` of that throwaway checkout,
  and the extension still loads). Registration removes both workarounds.
- **NCDatasets** is a hard dependency of SpeedyWeather, so NumericalRadiation's NetCDF reader
  extension (needed to load the ecCKD tables) is active whenever both packages are loaded;
  no extra dependency is needed for `ClearSkyEcCKDRadiation(spectral_grid, "32x32")`.
- **What the extension needs from NumericalRadiation**, all exported there:
  `AtmosphereProfile`, `ColumnGrid`, `SurfaceState`, `PhysicalConstants`,
  `LongwaveDiagnostics`, `solve_longwave!`, `AnalyticBandLongwave` (analytic-band adapter);
  `EcCKDTabulatedGasOpticsModel`, `ColumnAtmosphere`, `RadiativeFluxes`, `LongwaveOptics`,
  `ShortwaveOptics`, `CloudlessLongwave`, `CloudlessShortwave`, `ShortwaveColumnScratch`,
  `TabulatedSurfaceEmission`, `LongwaveBoundaryConditions`, `ShortwaveBoundaryConditions`,
  `optical_properties!`, `radiative_fluxes!`, `read_reference_ecckd_gas_optics`, `gas_names`
  (ecCKD). The code to move is ~500 lines in three files plus a ~270-line test file, a
  ~110-line docs example and three validation scripts (`git ls-tree mg/ecckd-speedy` in
  NumericalRadiation).

## Summary of changes

### Project.toml

```toml
[weakdeps]
NumericalRadiation = "cd8119b0-1744-44d6-9ede-6ad1ad750b26"

[extensions]
SpeedyWeatherNumericalRadiationExt = "NumericalRadiation"

[compat]
NumericalRadiation = "0.1.1"
```

Version stays `0.23.0-DEV` (already bumped on the base branch for the `Radiation` bundle;
the extension is additive).

### No `src` type

Both scheme types are NumericalRadiation's own: `AnalyticBandLongwave` (already there) and
`ClearSkyEcCKDRadiation` (added there on 2026-09-22 as the configured clear-sky ecCKD scheme: tabulated
gas optics, prescribed mole fractions of the gases the host does not carry with `o3` and `co2`
defaults, surface emissivity). Nothing in SpeedyWeather requires a radiation scheme to subtype
`AbstractRadiation`: the `Radiation` bundle has unconstrained type parameters and the model's
`radiation` component is dispatched on by the methods the extension defines. The one place
that needs a SpeedyWeather component is the time-step selection
(`get_prognostic_step(var, time_stepping, component::AbstractModelComponent)`); the extension
defines these methods for both scheme types, forwarding to the steps any radiation
parameterization gets (a zero-field `RadiationStandIn <: AbstractRadiation`), so the schemes
work both as `model.radiation` and inside a `Radiation` bundle. Earlier drafts had a
`WilliamsLongwave` wrapper and a SpeedyWeather-side `EcCKDRadiation` type with an
`ArgumentError` stub; see the revision log.

### The extension `ext/SpeedyWeatherNumericalRadiationExt/`

- `SpeedyWeatherNumericalRadiationExt.jl`: module, `using SpeedyWeather, NumericalRadiation`,
  `using DocStringExtensions` (a dependency of SpeedyWeather), imports of the
  SpeedyWeather internals it extends (`variables`, `initialize!`, `parameterization!`,
  `get_prognostic_step`, `get_tendency_step`, `pressure`, `pressure_half`,
  `flux_to_tendency`, `ParameterizationVariable`, `GridXYZ`, `Grid3D`, `Grid4D`) and of
  NumericalRadiation's staged API listed above; includes the two files below.
- `analytic_band_longwave.jl`: the analytic-band adapter as it is on `mg/ecckd-speedy` (constants
  from the host, stepped prognostics), as methods on NumericalRadiation's `AnalyticBandLongwave`
  (`variables`, `initialize!`, `parameterization!`, constructor from a `SpectralGrid`), used as
  `Radiation(SG; longwave = AnalyticBandLongwave(SG))`.
- `ecckd_radiation.jl`: the `ClearSkyEcCKDRadiation` kernel as it is on `mg/ecckd-speedy` after the merge with
  NumericalRadiation's `main` (lazy `TabulatedSurfaceEmission` surface sources, host
  `PhysicalConstants` through `ColumnAtmosphere.constants`, `ShortwaveColumnScratch` from
  seven column views, work arrays in the `:ecckd` namespace, stages
  `ecckd_surface_state`, `ecckd_column_atmosphere!`, `ecckd_column_optics`,
  `ecckd_longwave!`, `ecckd_shortwave!`, `ecckd_heating!`), as methods on NumericalRadiation's
  type plus constructors from a `SpectralGrid`. Identifier style follows this repository
  (`variables`, `radiation`, `Nz`).

The code is moved, not rewritten; the diff against `mg/ecckd-speedy` should be renames,
import lines and the split of the struct definitions into `src/`.

### Tests: `test/parameterizations/numericalradiation.jl`

Part of the normal suite (autodiscovered by `test/runtests.jl`, run by `Pkg.test` on Julia
1.10 and 1.13). The file starts with `using NumericalRadiation`; NumericalRadiation is a test
dependency in `test/Project.toml` with a `[sources]` git entry (`rev` = the NumericalRadiation
branch of PR #16 until merged/registered), and the 1.10 step of `CI_SpeedyWeather.yml` adds
it with `Pkg.add(url, rev)` before `Pkg.test`. The `ecrad_data` artifact (~30 MB) is fetched
by NumericalRadiation on first use and cached by `julia-actions/cache`.

What it checks (basic functionality only): the extension is active; `AnalyticBandLongwave`
constructs from a `SpectralGrid`, a longwave-only model runs one column update with finite
tendencies and positive OLR, and without a CO₂ component gives the same OLR as with one at
280 ppm; `ClearSkyEcCKDRadiation` constructs from the reference tables in the grid's number
format, its `:ecckd` work arrays have the expected sizes, one column update gives finite
tendencies, positive OLR and outgoing shortwave below the surface flux; a default wet model
with `ClearSkyEcCKDRadiation` runs two time steps. The full tests of the coupling (budget and
energy-conservation checks, a cross-check of one column against NumericalRadiation's staged
API, CO₂ forcing, night, 4× CO₂) live in NumericalRadiation.

Once NumericalRadiation 0.1.1 is registered, the `[sources]` entry and the `Pkg.add` line go
and it becomes an ordinary test dependency; the docs environment can then also execute the
example.

### Scripts and the full tests

Decided at review: the example, the validation scripts and the full unit tests of the
coupling stay in NumericalRadiation (its `examples/speedyweather_ecckd.jl`, `validation/`
and `test/speedyweather/`, whose environment pulls this branch of SpeedyWeather).
SpeedyWeather tests basic functionality only, see above.

### Documentation

- `docs/src/radiation.md`: a section "Radiation schemes from NumericalRadiation.jl":
  what the schemes are, how to construct them, the `:ecckd` work
  arrays, the clear-sky caveat and the ozone default. As a non-executed code block: the docs
  environment would otherwise need NumericalRadiation (git source; docs build on Julia 1.12,
  so possible) and the ecCKD tables, and an executed example adds minutes to the build.
  Decided: plain code block for now.
- No docstrings on this side: both scheme types are documented in NumericalRadiation, so the
  section names them in code font and links NumericalRadiation's documentation instead of
  `@ref`.
- CHANGELOG entry under `## Unreleased`.

### The NumericalRadiation side

NumericalEarth/NumericalRadiation.jl#16 (branch `mg/speedy-update`; a former stack of three
PRs, collapsed on 2026-09-24) carries everything this extension needs and nothing
SpeedyWeather-specific:

- `ClearSkyEcCKDRadiation` and `default_ozone_profile` (new `src/ecckd_radiation.jl`): the
  configured clear-sky scheme this extension dispatches on.
- `AtmosphereProfile` and `ColumnAtmosphere` with one array-type parameter per array, so the
  extension can pass views of different `SubArray` types (stepped prognostic arrays, 2D work
  arrays, interface arrays) without copies.
- An optional caller-owned `ShortwaveColumnScratch` argument of the clear-sky shortwave
  solver, so the per-column call in the fused kernel does not allocate.
- Its former extension `NumericalRadiationSpeedyWeatherExt` deleted, no `SpeedyWeather` weak
  dependency; the full coupling tests, the example and the validation scripts kept there,
  testing this branch of SpeedyWeather.

Its `main` gained an independent fix of the shortwave two-stream pole (`λμ₀ = 1`, 54b6a75
there) that this coupling had also guarded; `main`'s version stands.

## Testing and verification

1. `Pkg.test("SpeedyWeather")` includes `parameterizations/numericalradiation` (21 tests,
   ~1.5 min after precompilation); locally also
   `julia --project=SpeedyWeather/test SpeedyWeather/test/runtests.jl parameterizations/numericalradiation`.
2. The rest of the suite is unchanged by the extension: `src/` is untouched apart from
   `Project.toml`, and without NumericalRadiation loaded neither scheme type exists.
3. NumericalRadiation's own coupling tests (54 tests: budgets, energy conservation, a
   column cross-checked against its staged API, CO₂ forcing, night) run there against this
   branch, in its `test/speedyweather/` environment.
4. CI: `CI_SpeedyWeather.yml` on Julia 1.10 (git dependency added by the 1.10 step) and 1.13.

## Documentation changes

As listed above: `docs/src/radiation.md` section, docstrings, CHANGELOG, this plan.

## Known limitations

- Until NumericalRadiation 0.1.1 is released (0.1.0 predates PR #16), users install it from
  that branch, the test environment pulls it from git, and the docs example is not executed.
- The ecCKD scheme is clear-sky; ozone comes from an analytic default profile
  (`default_ozone_profile` in NumericalRadiation, a Chapman-layer shape peaking at 8 ppmv near
  30 hPa) because SpeedyWeather has no ozone field; CO₂ is the only gas taken from the model.
- GPU: untested (Metal fails earlier in SpeedyWeather; needs a CUDA machine).

## Future work

- After NumericalRadiation 0.1.1 is released: drop the git source in `test/Project.toml` and
  the `Pkg.add` line of the 1.10 CI step, and execute the docs example.
- Prescribed ozone, as a SpeedyWeather component modelled on `greenhouse_gases` (a zonal-mean,
  pressure-dependent climatology filled at `initialize!`) or as a prescribed tracer; the
  scheme then reads a per-column ozone profile instead of `mole_fractions.o3`.
- A radiation call frequency: parameterizations that run every `N` steps and hold their
  tendency in between; a 32×32 g-point column per step dominates cost at climate resolution.
- Clouds for the ecCKD scheme (a cloudy counterpart of `ClearSkyEcCKDRadiation` in
  NumericalRadiation, which already has cloud-overlap solvers) once SpeedyWeather has a cloud
  state to feed them.
