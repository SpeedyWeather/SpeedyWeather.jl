# Sea-ice-aware surface longwave emission

> Status: **completed**. #1406 is option 1 plus the review changes; the energy conservation of
> `ThermodynamicSeaIce` (option 2) moves to a separate PR, as agreed on 2026-11-17.

Date of initial draft: 2026-11-15

Base revision: 7a76967fecd8308dbcfe0db52e35ef3a85ce9bde (main merged into #1406 by the maintainers on 2026-11-16;
first draft on dd988ccb9c3e6c87d6a0a8c0dde3e783484b5cf7)

## Originating prompt

> Would you consider any of these changes due to bugs in SpeedyWeather that have been unreported?

> Test it first [...] it is important that we are careful and check this thoroughly.

> Rerun the tests and do the Earth runs for all four

> Please remove option 3 from the draft. I do not think it reasonable to suggest his climate model break
> conservation of energy.

> The PR for Speedy Weather, it appears, should have a verification section in the plan. I should compare
> nonlinear radiative blending values using the code's linear area-weighted effective temperature and true
> grid-cell average flux. I should also show that the method used to calculate Bulk Richardson Drag by means
> of liquid seawater temperature underneath the ice instead of the cold ice skin temperature is of little
> consequence due to the insulation factor. I should unit test mathematical invariants for open-water limit,
> cold-air sea ice limit, warm air cap and land isolation.

> Another part of verification can be refactoring of code. In this case, I think a subagent tasked with
> reviewing the code in the new PR once it is drafted would be appropriate [...]

> (2026-11-17, after the #1414 finding) Yes, split it. [...] Go ahead and continue the split. Let's try and
> wrap up our end of the work tonight.

> (2026-11-17) Yes, include the Richardson check and the Jensen tests

## Revision log

- 2026-11-14: #1403 reported. The surface longwave over sea ice is emitted at the sea surface temperature,
  the water under the ice (about −1.8 °C), however cold the air above it.
- 2026-11-14: the fix was audited across all four ocean-emitting schemes (`OneBandLongwave`,
  `JeevanjeeRadiation`, and `AnalyticBandLongwave` and `ClearSkyEcCKDRadiation` from the
  NumericalRadiation extension). There is no separate sea-ice temperature, as asked for in #1403.
  #1406 was opened.
- 2026-11-14: the Earth A/B runs show a side effect with `SlabOcean` + `ThermodynamicSeaIce`: summer sea ice
  vanishes. #1406 offered two options: as is (1), or with a conductive heat loss through the ice (2). An option
  that does not conserve energy was considered and rejected.
- 2026-11-15: review of #1406. Option 2 was asked for in the same PR, plus:
  - dispatch the emission temperature on `model.sea_ice`;
  - rely on the parameterization order instead of the `> 0` check;
  - rename the helper to `ocean_surface_emission_temperature`;
  - name all four schemes in the CHANGELOG;
  - add this plan;
  - comment on the test SST.

  The version is already `0.25.0-DEV` (#1262).
- 2026-11-15 to 2026-11-17: option 2 was designed, signed off (2026-11-16) and implemented: ice volume as the
  prognostic variable, energy-conserving freezing and melt, conduction through the ice from a diagnostic ice
  surface temperature, and a local radiation split. The 3-year Earth A/B runs with it restore the summer ice.
- 2026-11-17: #1414's review found that `ThermodynamicSeaIce` on main is an energy source: once ℵ = 1 the heat
  the water loses is thrown away with the clamped ice, and melting takes no heat from the water. The maintainers
  proposed keeping #1406 as option 1, physically right on the atmosphere side, and fixing the sea-ice conservation
  in a separate PR. Agreed. The option 2 work moves to that PR with its own plan. #1406 keeps the review changes.
  Making `sea_ice` a core component (below) stays in #1406, as the dispatch needs it on the GPU.

## Problem description

`OneBandLongwave` and the other ocean-emitting schemes compute the upward longwave over the ocean part of a
grid point as ε σ SST⁴. Under sea ice the SST is the water below the ice, held near freezing. A winter ice
surface under −30 °C air therefore radiates as if it were at −1.8 °C. This overstates the upward longwave by
about 100 W/m² (with ε = 1, at −30 °C against −1.8 °C).

## Background

- **Emission temperature:** `T_eff = (1 − ℵ)·SST + ℵ·min(SST, T_sa)`.
  - ℵ is the sea ice concentration.
  - `T_sa` is `vars.parameterizations.surface_air_temperature`, the lowermost layer temperature lowered
    dry-adiabatically to the surface (`surface_condition.jl`, in `:boundary_layer`, which runs before `:radiation`).
  - The ice never emits warmer than the water under it.
  - Over open water, and without sea ice, `T_eff = SST` exactly.
- **`SlabOcean`** subtracts the same upward longwave, so the surface exchange conserves energy.
- **`ThermodynamicSeaIce`** on main is not energy-conserving itself (#1414), which is why the sea ice climate
  changes with this PR. That is fixed separately.

## Summary of changes

- `ocean_surface_emission_temperature(ij, vars, sst, scheme, model)` (was `ocean_skin_temperature`), used by
  all four schemes, dispatches on `model.sea_ice`. `PrescribedSeaIce`, `ThermodynamicSeaIce` and `nothing` all use
  the blend above. With `nothing` the concentration can still be set manually (`docs/src/sea_ice.md`), and the
  albedo and the surface fluxes react to it, so the longwave does too. Without an allocated concentration the
  blend is the SST. A future sea ice model with its own ice surface temperature only needs a new method.
- The value check `ifelse(surface_air_temperature > 0, ...)` is removed. The order (`:boundary_layer` before
  `:radiation`) is documented in the docstring and in `radiation.md`, and tested for the default tuples. The
  compile-time fallback to the SST when `surface_air_temperature` is not allocated stays.
- `:sea_ice` is a core component of `PrimitiveWetModel` and `PrimitiveDryModel`, because the column kernels on the
  GPU only see the core components. `ThermodynamicSeaIce` is therefore immutable (nothing mutated it) and
  `Adapt`ed. `PrescribedSeaIce` is a singleton.
- CHANGELOG names all four schemes. `docs/src/radiation.md` documents the emission temperature and the dispatch.

## Testing and verification

All tests run on the CPU with Julia 1.13 and `--check-bounds=yes`. All pass.

### 1. Unit tests: invariants of the emission temperature (`longwave_radiation.jl`)

For `OneBandLongwave` and `JeevanjeeRadiation`:

| Invariant | Test |
|---|---|
| Open-water limit | ℵ = 0 gives exactly εσSST⁴, for any surface air temperature (`==`) |
| Cold-air sea ice limit | ℵ = 1 and T_sa < SST give εσT_sa⁴ |
| Warm-air cap | T_sa > SST gives εσSST⁴ (ice never emits warmer than the water) |
| Half-covered cell | ℵ = 0.5 gives the linear blend |
| Bounds | `min(SST, T_sa) ≤ T_eff ≤ SST` over a grid of ℵ, T_sa and SST |
| Monotonicity | more ice is not warmer; warmer air is not colder, and warms T_eff by at most ℵ·ΔT_sa |
| Number format | Float32 in gives Float32 out |
| Land isolation | the land emission is identical with and without sea ice, at every grid point |
| Energy closure | at every grid point, the upward longwave the atmosphere receives is `lf·land + (1 − lf)·ocean`, the fluxes the land and the slab ocean subtract |
| Jensen gap | the blend is below the true cell average `εσ[(1 − ℵ)SST⁴ + ℵT_sa⁴]`, and the gap matches `6εσT̄²ℵ(1 − ℵ)ΔT²` within 5 % (four cases, ΔT 10–60 K) |
| Dispatch | `nothing` (with a set concentration), `PrescribedSeaIce` (instance and a model, no time step dimension) and `ThermodynamicSeaIce` |
| Order | `:boundary_layer` before `:radiation` in `PrimitiveWetModel` and `PrimitiveDryModel` (replaces the `> 0` check) |

`numericalradiation.jl`: for `AnalyticBandLongwave` and `ClearSkyEcCKDRadiation`, ice emits less than open water
and land is unaffected (test SST 295 K everywhere, commented).

The related test files pass too:
- `ocean_sea_ice.jl`, `radiation.jl`, `albedo.jl`, `surface_fluxes.jl`, `surface_condition.jl`;
- `surface_roughness.jl`, `ocean_land_time_stepping.jl`, `all_parametrizations.jl`, `variables.jl`.

`radiation.jl` needs its own module when it runs after `numericalradiation.jl`, because NumericalRadiation exports a
`OneBandShortwave` of its own.

### 2. Jensen gap and surface Richardson number in an Earth run

Default `PrimitiveWetModel` (`OneBandLongwave`, `SlabOcean`, `ThermodynamicSeaIce`), T32 L8, this PR, 3 years.
Sampled once a day in the last year, over the ocean part of every point with sea ice. The script, `earth_verify.jl`,
is outside the repository.

- **Jensen gap of the linear blend** (true cell-average flux minus the emitted flux):
  - the monthly mean over the ice points is 0.02–1.7 W/m², largest in the Arctic in January. That is 0.01–0.7 % of
    the emitted flux;
  - as a global mean it is at most 0.06 W/m²;
  - the largest at any point on any day is 10.8 W/m²;
  - the second-order formula reproduces it to within 0.03 W/m² on average.

  It stays a known limitation. A flux-equivalent temperature, `((1 − ℵ)SST⁴ + ℵT_ice⁴)^¼`, would remove it. That
  belongs with a sea ice model that has its own ice temperature.
- **Surface Richardson number over ice.** The model uses the SST, and the saturation humidity at the SST, for the
  ocean surface, also under ice. The alternative uses the ice surface temperature `min(SST, T_sa)` on the ice
  fraction ℵ:
  - the ice skin would make the surface layer stable at 1–43 % of the ice points (monthly), but the drag
    coefficients barely change (momentum 1.211 → 1.210 × 10⁻³ in the largest case);
  - wind stress changes by at most 0.06 % in the monthly mean, and at most 0.0027 N/m² at any point;
  - the sensible heat flux changes by at most 0.003 W/m² in the monthly mean, and at most 0.19 W/m² at any point.
    Without the insulation the monthly change would be up to 0.26 W/m², so the insulation removes about 99 % of it;
  - the reconstructed SST-based drag matches the model's to within 0.6 %.

  So the SST-based Richardson number matters little here. For heat that is because of the insulation, and for
  momentum because the stability correction is small at these wind speeds.

### 3. GPU

`test/GPU/primitive_wet.jl`, "GPU surface longwave over sea ice matches CPU" (run by Buildkite): with
`OneBandLongwave` and `JeevanjeeRadiation`, each with `PrescribedSeaIce` and `ThermodynamicSeaIce`, sea ice on every
other point, the ocean surface longwave on GPU equals the CPU's within `rtol = 1e-5`. Passes locally (RTX 4050, 8/8).
`get_sea_ice` is tested directly in `longwave_radiation.jl`.

Check by hand beforehand:

NVIDIA RTX 4050, CUDA, with random SST, surface air temperature and sea ice on half the points:

- The ocean surface longwave differs between GPU and CPU by at most 1.8 × 10⁻⁴ W/m² over ice and 3.1 × 10⁻⁵ W/m²
  over open water, at most 6.3 × 10⁻⁷ relative. That is Float32 round-off: open water, where the code is
  unchanged, differs too.
- 20 full-model GPU time steps stay finite, both with `OneBandLongwave` and with `JeevanjeeRadiation`.

### 4. Earth A/B

Default `PrimitiveWetModel` (`OneBandLongwave`), T32 L8, last of 3 years, base 7a76967, ice in 10⁶ km²:

| | main | this PR |
|---|---|---|
| Arctic sea ice Mar / Sep | 17.2 / 0.0 | 16.5 / 0.0 |
| Antarctic sea ice Feb / Sep | 8.2 / 15.7 | 0.3 / 15.1 |
| Winter air poleward of 60° N / S [°C] | −17.3 / −18.2 | −19.9 / −20.5 |

The loss of Antarctic summer ice comes from `ThermodynamicSeaIce`, where melting costs no heat (see #1414).

### 5. Benchmark

The repository's regression check (`benchmark/regression/regression.jl`, CPU and GPU) ran on 2026-11-15 for
dd988cc. It found no regression: 0.996 of main on CPU and 1.014 on GPU. The changes since add no per-step work: the
dispatch is resolved at compile time, and a value check was removed. So it was not rerun.

### 6. Independent code review

A sub-agent without the development context reviewed the diff and found no correctness bugs. Applied:
- `nothing` also blends a manually set concentration;
- the order requirement for custom `parameterizations` tuples is documented;
- the argument is named `scheme`;
- a `PrescribedSeaIce` model test;
- shorter lines;
- a shorter CHANGELOG entry;
- a cross-link from `sea_ice.md`.

A GPU testset was added afterwards too (see section 3).

### Not done

- Enzyme gradients at the `min` kink (T_sa = SST). The blend itself has no branches: a `min` and a linear
  combination.
- The full test suite, and Reactant.

## Documentation changes

- `docs/src/radiation.md`: the emission temperature over sea ice and its dispatch on the sea ice model.
- CHANGELOG.

## Known limitations

- The emission is blended linearly in temperature. Because T⁴ is convex, εσT_eff⁴ underestimates the true
  cell-average flux `ε σ [(1 − ℵ) SST⁴ + ℵ T_ice⁴]` by about `6 ε σ T̄² ℵ(1 − ℵ) ΔT²`, about 5ε W/m² for
  ℵ = 0.5 and ΔT = 30 K. In the Earth run it is at most 1.7 W/m² in the monthly mean over the ice (testing,
  section 2).
- `T_sa` comes from a dry-adiabatic extrapolation. It is warm-biased in surface inversions, so the correction is
  conservative.
- A custom `parameterizations` tuple without `:boundary_layer` before `:radiation` leaves the surface air
  temperature unset (zero), so ice would emit almost nothing. This is documented, not checked.
- The surface Richardson number over ice still uses the SST. In the Earth run that changes the wind stress by
  under 0.1 % and the sensible heat flux by at most 0.003 W/m² in the monthly mean (testing, section 2).
- `ThermodynamicSeaIce` is not energy-conserving (#1414).

## Future work

- Separate PR: energy-conserving `ThermodynamicSeaIce` (melting takes latent heat from the ocean, no ice thrown
  away at full cover), then conduction through the ice from a diagnostic ice surface temperature, which would add
  a method of `ocean_surface_emission_temperature` for it.
- An ice-skin surface Richardson number over sea ice.
