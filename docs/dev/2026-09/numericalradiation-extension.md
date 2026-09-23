# A NumericalRadiation extension of SpeedyWeather

> Status: **in progress**. Reviewed 2026-09-22 (decisions in the revision log); implemented on
> `mg/numericalradiation-extension`, uncommitted, local tests pending. Companion plan on the
> NumericalRadiation side: `docs/plans/extension_upstream.md` there.

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
  workflow collapse into `test/Project.toml` (TODO noted in `test/numericalradiation/Project.toml`).
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
- **NumericalRadiation is not registered.** It can be a weak dependency (only a UUID is
  needed) but any environment that installs it needs a `[sources]` git entry, which requires
  Julia ≥ 1.11, or `Pkg.add(url = ...)`. `CI_SpeedyWeather.yml` runs Julia 1.10 and 1.13 and
  already has a 1.10-only step working around missing `[sources]` support. Listing an
  unregistered package in `test/Project.toml` would break `Pkg.test` on 1.10.
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
NumericalRadiation = "0.1"
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

### Tests: `test/numericalradiation/` (own environment)

Because NumericalRadiation cannot be a plain test dependency (unregistered, see Background):

- `test/numericalradiation/Project.toml`: `NumericalRadiation` (git source, URL of
  NumericalEarth/NumericalRadiation.jl, `rev` = the branch carrying the package-side changes
  until they are on `main`), `SpeedyWeather` and the monorepo siblings by path, `Statistics`,
  `Test`. Needs Julia ≥ 1.11 for the git source.
- `test/numericalradiation/runtests.jl` includes `basic.jl` (see below).
- `test/runtests.jl`: add `numericalradiation/` to the prefix filter so `Pkg.test` does not
  pick these files up.
- `.github/workflows/CI_NumericalRadiation.yml`: modelled on `CI_Enzyme.yml`; Julia 1.11,
  instantiate `test/numericalradiation`, run its `runtests.jl`; skip labels as the other
  workflows. The `ecrad_data` artifact (~30 MB download) is fetched by NumericalRadiation
  on first use and cached by `julia-actions/cache`.

Alternative (**decision for review**): register NumericalRadiation in General, after which it
becomes an ordinary entry of `test/Project.toml` and `docs/Project.toml`, the test file moves
to `test/parameterizations/`, and the dedicated workflow and environment go away. The
dedicated environment is the right shape until then and is a small change to undo.

### Scripts and the full tests

Decided at review: the example, the three validation scripts and the full unit tests of the
coupling stay in NumericalRadiation (`examples/`, `validation/`, `test/speedyweather/`, whose
environment pulls this branch of SpeedyWeather). SpeedyWeather tests basic functionality only
(`test/numericalradiation/basic.jl`: both schemes construct, models build, one column update
and two time steps give finite output, the `:ecckd` work arrays exist) and, in the main suite,
that the constructors throw the explanatory error without NumericalRadiation.

### Documentation

- `docs/src/radiation.md`: a section "ecCKD and analytic-band radiation from
  NumericalRadiation.jl": what the schemes are, how to construct them, the `:ecckd` work
  arrays, the clear-sky caveat and the ozone default. As a non-executed code block: the docs
  environment would otherwise need NumericalRadiation (git source; docs build on Julia 1.12,
  so possible) and the ecCKD tables, and an executed example adds minutes to the build.
  Decided: plain code block for now.
- Docstrings of the two `src/` types are picked up by the API page through `@autodocs`.
- CHANGELOG entry under `## Unreleased`.

### On the NumericalRadiation side (its own plan, `docs/plans/extension_upstream.md`)

Delete the extension, its test and environment, the example and the validation scripts;
drop the SpeedyWeather weak dependency, extension entry, compat and CI job; rewrite the
README section to point here. Everything in `src/` stays.

## Testing and verification

1. `test/numericalradiation/runtests.jl` locally against the NumericalRadiation branch, one
   Julia process (the suite takes ~2 min after precompilation; 53 tests on the source
   branch).
2. `Pkg.test("SpeedyWeather")` unaffected: no new test dependency, the new files are filtered
   out of autodiscovery, and the `src/` types add only two struct definitions and two error
   methods (checked by the existing radiation tests, which construct `Radiation` and the
   primitive models).
3. Without NumericalRadiation loaded neither scheme type exists; SpeedyWeather's own test
   suite is unchanged by the extension.
4. The dedicated workflow on the PR.

## Documentation changes

As listed above: `docs/src/radiation.md` section, docstrings, CHANGELOG, this plan.

## Known limitations

- Until NumericalRadiation is registered, users install it with
  `Pkg.add(url = "https://github.com/NumericalEarth/NumericalRadiation.jl")` and the tests run
  only on Julia ≥ 1.11 in their own workflow.
- The ecCKD scheme is clear-sky and takes ozone from an analytic default profile; both are
  tracked on the NumericalRadiation side (their plan, items U2 and clouds) and unchanged by the
  move.
- GPU: untested (Metal fails earlier in SpeedyWeather; needs a CUDA machine).

## Future work

- Register NumericalRadiation and collapse the dedicated test environment.
- Prescribed ozone (a SpeedyWeather component or a tracer) and a radiation call frequency;
  both were planned as upstream items U2 and U3 in NumericalRadiation's plan and belong here
  now.
- Clouds for the ecCKD scheme (a cloudy counterpart of `ClearSkyEcCKDRadiation` in
  NumericalRadiation) once SpeedyWeather has a cloud state to feed them.
