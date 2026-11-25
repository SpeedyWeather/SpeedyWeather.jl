# Swappable saturation vapor pressure formulations

> Status: **completed**. Four formulations with `grad` implemented as a field of `EarthAtmosphere` in #1441; large-scale condensation uses `grad`, convection not yet.

Date of initial draft: 2026-11-24

Base revision: 75bf89d9aeac994eb8857a06d136a5ceafec0099

## Originating prompt

> Can you have a look at our saturation vapour pressure function. I know there's several
> Clausius-Clapeyron approximations out there, which one are we using and is that actually the best?

> Would you suggest to hardcode the derivative of CC in the condensation or revisit what I had
> started for Tetens to have a more general grad function that could be inlined depending on the
> actual CC approximation used?

> And the svp type could live in `EarthAtmosphere`?

> Can you create locally a branch, write a plan, then write the 4 implementations Constant-L,
> linear L (Romps), IFS and the Bolton version of Magnus/Tetens each with a `grad` function for it.
> Don't adapt the large-scale condensation yet to use this, at the moment the
> `saturation_vapor_pressure` function interface should stay the same and no changes to the
> parameterizations for now. Do remove the `saturation_vapor_pressure` field but only do +DEV,
> don't consider this one field being a public interface. Don't push that branch yet, just keep it
> local, write some tests and find a place in the documentation to write a summary about this
> down + user interface.

## Revision log

> - Naming by assumption: "Clausius-Clapeyron" is the exact phase-change relation, our formula is its
>   constant-latent-heat integration, hence `ClausiusClapeyronConstantL` / `ClausiusClapeyronLinearL`.
> - `EarthAtmosphere.saturation_vapor_pressure` (e₀) removed without deprecation, it moves into
>   `ClausiusClapeyronConstantL`; treated as non-public, so only a `+DEV` bump.
> - Code moved from `dynamics/` to `parameterizations/thermodynamics.jl` (user suggestion), still
>   included before `dynamics/atmosphere.jl`; tests in `test/parameterizations/thermodynamics.jl`.
> - Docs: new page `docs/src/saturation_vapor_pressure.md` under Physics.
> - Condensation initially kept its hardcoded `dq_sat/dT = q_sat Lᵥ/(Rᵥ T²)`.
> - 2026-11-25, after merging #1431 (reevaporation): "update the condensation code to use the `grad`
>   function". New `grad_saturation_humidity(T, p, A) = ε/p deₛ/dT`, used for both the implicit
>   condensation divisor and the reevaporation limiter. Convection's moist adiabat still implicitly
>   assumes constant L (left for a follow-up).

## Problem description

`saturation_vapor_pressure(T, ::EarthAtmosphere)` integrates Clausius-Clapeyron with a constant
latent heat of condensation. Compared with reference data (Murphy and Koop, 2005) it overestimates
eₛ by ~2.4% at 30˚C and ~4% at 40˚C (because Lᵥ decreases with temperature), and it has no ice
phase (+23% at -20˚C, +58% at -40˚C relative to saturation over ice). The formula is not
swappable, and the derivative used by the implicit condensation is hardcoded for this one formula.
An unused `TetensEquation` exists whose `grad` has a sign error in the denominator.

## Background

Exact Clapeyron: `deₛ/dT = L / (T (v_vapor - v_liquid))`. Neglecting the condensed volume, assuming
an ideal gas and a constant L gives `eₛ = e₀ exp(L/Rᵥ (1/T₀ - 1/T))`. Letting L vary linearly with
temperature (constant heat capacities, Kirchhoff) gives the Romps (2017) form. Magnus/Tetens-type
formulas `e₀ exp(a (T - T₀)/(T - b))` are empirical fits (Bolton 1980; IFS uses separate water and
ice coefficients with a quadratic blend between -23˚C and 0˚C).

## Summary of changes

- `abstract type AbstractSaturationVaporPressure` with
  `saturation_vapor_pressure(svp, T, atmosphere)` and `grad(svp, T, atmosphere)` = deₛ/dT [Pa/K].
  The atmosphere argument lets a formulation take consistent constants (Lᵥ, Rᵥ) from it.
- Four implementations (all isbits, `NF`-parametric, GPU-safe):
  1. `ClausiusClapeyronConstantL` (default, e₀ = 610.78 Pa, T₀ = 273.15 K, Lᵥ and Rᵥ from the
     atmosphere) — bit-identical to the current formula.
  2. `ClausiusClapeyronLinearL` (Romps 2017, over liquid, own constants at the triple point).
  3. `IFSMixedPhase` (Tetens over water and ice, quadratic blend, ECMWF IFS Cy47r3 Part IV).
  4. `MagnusBolton` (Bolton 1980, over liquid).
- `EarthAtmosphere{NF}` → `EarthAtmosphere{NF, SVP}` with new field `saturation::SVP`;
  field `saturation_vapor_pressure` (e₀) removed. Convenience constructors default to
  `ClausiusClapeyronConstantL(NF)`.
- `saturation_vapor_pressure(T, A)` and `saturation_humidity(T, p, A)` keep their signatures and
  forward to `A.saturation`. New `grad_saturation_humidity(T, p, A)` used by `ImplicitCondensation` (implicit divisor and reevaporation limiter) instead of the hardcoded Clausius-Clapeyron derivative.
- `TetensEquation` kept as is but its `grad` sign error fixed and its docstrings corrected.
- SpeedyWeather version `0.26.0-DEV` already carries a `-DEV` tag (unreleased minor bump), so no further bump is needed.

## Testing and verification

`test/parameterizations/thermodynamics.jl`:
- each formulation: value against reference values (Murphy and Koop 2005) within a tolerance,
  `grad` against centred finite differences in Float64 over -60…45˚C,
- default atmosphere bit-identical to the old formula, `saturation_humidity` forwards,
- continuity of `IFSMixedPhase` at the blend boundaries,
- Float32 type stability, `TetensEquation` grad fix.

## Documentation changes

New page "Saturation vapor pressure" under Physics with formulas, a comparison and the user
interface (`EarthAtmosphere(spectral_grid, saturation = MagnusBolton(spectral_grid))`).

## Known limitations

- `ClausiusClapeyronLinearL` and `MagnusBolton` are over liquid only.
- The formulations are not `@parameterized`, so they are not part of `parameters(model)`.
- Convection's moist-adiabatic lapse rate still implicitly uses the constant-L derivative.

## Future work

- Use `grad` in convection's moist adiabat; possibly a fused value-and-derivative to save an `exp`.
- Ice branch for Romps; decide on a new default.
- Consistent latent heat L(T) in the condensation heating.
