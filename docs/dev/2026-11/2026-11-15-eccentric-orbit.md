# Eccentric orbits: Kepler orbit, Sun–planet distance and precession

> Status: **completed**. Implemented as agreed in #1407; unit tests pass and daily insolation matches climlab to 0.4 W/m² RMS.

Date of initial draft: 2026-11-15

Base revision: caed82d647c38b6fd793ff49ff674e62d437b8cc

## Originating prompt

> Full Milankovitch cycles: obliquity, eccentricity and precession, with the precession
> chosen as the month of perihelion.

(From Paleo Planet Sim, a worldbuilding app that runs SpeedyWeather for user-designed
planets. Prompted by teaching uses like Build Your Own Earth, Schultz et al. 2017, BAMS,
whose orbital experiments vary tilt, eccentricity and the orbit through the last glacial
cycle.)

## Revision log

- 2026-11-15: initial plan with flat `planet.eccentricity` / `planet.longitude_of_perihelion`
  parameters and a `KeplerOrbit` component.
- 2026-11-15 (review in #1407):
  > What about making planet.orbit by default CircularOrbit or alternatively KeplerOrbit with
  > eccentricity and perihelion? Then we can dispatch over the orbit in the zenith calculation
  > to compute the distance factor?

  and

  > I wouldn't put it into the radiation code as every grid cell would redo that ij-independent
  > orbit calculation. Maybe make it part of the zenith's global parameterization! call? And then
  > dispatch over planet.orbit? [...] Require the distance variable as ref you can then write
  > into and read from inside kernels.

  The design changed to `planet.orbit::AbstractOrbit` (`CircularOrbit` default, `KeplerOrbit`).
  The orbit position is computed once per time step in `parameterization!(vars, ::AbstractZenith,
  model)`, outside any kernel, and the distance factor is stored in a `ScalarDim` parameterization
  variable that the shortwave schemes read. The equation of time from the orbit is left as future
  work.
- 2026-11-15 (implementation):
  - `cos_zenith!` takes the declination as an optional last argument, defaulting to today's
    `S.solar_declination`.
  - `solar_declination(T, zenith, orbit_time, planet)` follows the orbit; seasonal ozone uses it.
  - The NumericalRadiation ecCKD shortwave boundary also gets the distance factor.
  - With e = 0, `KeplerOrbit` uses the exact sin δ = sin(tilt)·sin λ, while `SinSolarDeclination`
    uses δ = tilt·sin λ. They differ by up to ~0.004 rad, so `CircularOrbit` keeps the existing
    path rather than routing through Kepler.
  - The equinox phase uses the mean year (365.25 days), like `SinSolarDeclination`, so the same
    leap-year offset of up to about half a day applies.
- 2026-11-16 (review in #1407, Milan and milankl-claude):
  > - planet.orbit is responsible to define the orbit and needed variables for its computation
  > - model.solar_zenith is responsible to execute that computation before the zenith angle is computed
  > - model.radiation uses zenith + distance and falls back to what we currently have for CircularOrbit

  The structure changed accordingly:
  - **Variables:** `variables(planet::Earth)` forwards to `variables(planet.orbit)`. Only
    `KeplerOrbit` declares `solar_distance_factor`; `CircularOrbit` declares nothing.
  - **Zenith:** `orbit_position!(vars, T, orbit, zenith, planet, orbit_time)` returns the
    declination and writes the factor for `KeplerOrbit` only.
  - **Radiation:** `solar_distance_factor(vars)` returns the stored factor, or 1 when the variable
    doesn't exist. That's decided at compile time through `haskey` on the NamedTuple, so nothing
    is read inside the kernels for a circular orbit. It is used at all four sites.
  - **`KeplerOrbit`:**
    - it is `@parameterized`, with `[OPTION]` fields and `@param eccentricity (bounds = 0 .. 1,)`;
    - its defaults are the constants `EARTH_ECCENTRICITY` and `EARTH_LONGITUDE_OF_PERIHELION`;
    - its docstring notes that `KeplerOrbit()` differs from `CircularOrbit()` in insolation
      (±3.4%) and that the equation of time stays present-day Earth's;
    - it has a `[DERIVED]` `mean_anomaly_equinox`, computed by `mean_anomaly_at_equinox`, which
      also throws a `DomainError` unless 0 ≤ e < 1;
    - `KeplerOrbit(spectral_grid; ...)` is a constructor.
  - **Equinox phase:** `equinox_angle` is shared by `SinSolarDeclination` and `KeplerOrbit`.
  - **Warning:** a warning at initialization when a `KeplerOrbit` meets a non-default
    `solar_declination`, which it ignores.
  - **Tests:**
    - the ozone test compares `ozone_absorption!` for `CircularOrbit` and an eccentric orbit
      with perihelion at the December solstice. There the solstice falls on 23 January
      instead of 22 December, so on 1 February ozone is higher in the north and lower in the
      south;
    - variable allocation is tested per orbit;
    - the simulation tests use T21 L4 for 1 day.
  - The Newton iteration count stays a runtime keyword (5 iterations once per step, no
    measurable gain from compile time).
- 2026-11-17/18 (Reactant GPU failure in #1407): a traced value can't be written into the plain
  `Ref{NF}` of a `ScalarDim`, so the distance factor and CO₂ moved to a new `ArchitectureScalarDim`
  allocated with `scalar_ref` (moved to SpeedyWeatherInternals).
- 2026-11-19 (review in #1407):
  > We want to consolidate both of them just in `ScalarDim`

  `ArchitectureScalarDim` is merged back into `ScalarDim`, which now allocates with `scalar_ref`.
  With Reactant, `prognostic.scale` is then traced in the compiled step too, so `scale_prognostic!`,
  `scale!`, `unscale!` and the `NCycleLorenz` update take `scale::Number` instead of `Real`.
- 2026-11-19 (JET test in CI): with the `KeplerOrbit` default the `solar_distance_factor` Ref is one
  more variable type in `vars.parameterizations`, enough for `reset_variable!(vars, var::Symbol, ...)`
  to infer `getfield` as a widened union and dispatch at runtime. It now takes the name as `Val(var)`.

## Problem description

SpeedyWeather's orbit is circular:

- `SinSolarDeclination`: δ = tilt · sin(g − equinox), with the year angle g advancing
  uniformly with `orbit_time`. Seasons all have equal length.
- The shortwave schemes use `planet.solar_constant * cos_zenith` (three places in
  `shortwave_radiation.jl`), so TOA sunlight never depends on the Sun–planet distance.
- `SolarDeclination` (Spencer 1971) and `SolarTimeCorrection` carry Earth's real orbit
  implicitly, but only for Earth's tilt, eccentricity and perihelion date. They can't be
  changed.

So none of these can be run:
- Milankovitch experiments: mid-Holocene, the Eemian at 126 ka, glacial inception at 116 ka.
  PMIP-style orbital forcing is mostly eccentricity × precession.
- Mars (#838), where e = 0.093 gives a ±20% insolation swing.
- Eccentric exoplanets (#462).
- Earth's own ±3.4% perihelion–aphelion asymmetry.

## Background

- Berger (1978), *J. Atmos. Sci.* 35, 2362–2367: orbital elements for the last few Myr as a
  short trigonometric series. Berger & Loutre (1991) for the insolation quantities.
- The PMIP protocol (Otto-Bliesner et al. 2017, *GMD*): the vernal equinox is fixed at
  21 March and the orbit is set by (e, obliquity, longitude of perihelion ϖ measured from the
  vernal equinox). Present day: e ≈ 0.0167, with perihelion in early January (ϖ − 180° ≈ 102°
  in Berger's convention).
- climlab (`climlab.solar.orbital`, `daily_insolation`; MIT licence): a reference
  implementation to test against.

## Summary of changes

1. **Orbit types** (`dynamics/planet.jl`):

   ```julia
   abstract type AbstractOrbit end
   struct CircularOrbit <: AbstractOrbit end
   @kwdef struct KeplerOrbit{NF} <: AbstractOrbit
       eccentricity::NF = 0.0167
       longitude_of_perihelion::NF = 282.9    # [˚] from the vernal equinox
   end
   ```

   `Earth` gains `orbit::O = CircularOrbit()`. The vernal equinox (`planet.equinox`) fixes the
   calendar, as in PMIP, and the year angle of `orbit_time` is the mean anomaly clock.
2. **Orbit position** (`radiation/zenith.jl`):
   `orbit_position(T, orbit, zenith, planet, orbit_time) -> (δ, distance_factor)`.
   - `CircularOrbit`: δ from the zenith's `solar_declination` component as now, with a distance
     factor of 1.
   - `KeplerOrbit`:
     - the mean anomaly is M = M₀ + (g − g_equinox), where M₀ is the mean anomaly at the vernal
       equinox (true anomaly ν₀ = −ϖ);
     - Kepler's equation E − e·sin E = M is solved with 5 fixed Newton iterations;
     - then ν = 2 atan(√(1+e) sin(E/2), √(1−e) cos(E/2)) and λ = ν + ϖ;
     - sin δ = sin(tilt) sin λ and (a/r)² = 1/(1 − e cos E)².

   `solar_declination(T, zenith, orbit_time)` gains a `planet` method so that seasonal
   ozone follows the same orbit.
3. **Distance factor:** `ParameterizationVariable(:solar_distance_factor, ScalarDim(1))` is
   declared by the zenith. It is written in the zenith's global `parameterization!` and read by
   the three `solar_constant * cos_zenith` sites. `cos_zenith` stays the geometric cosine.

`CircularOrbit` keeps exactly today's path, so default simulations are unchanged.

## Testing and verification

- **Unit tests** (`test/parameterizations/orbit.jl`, 908 pass):
  - the Kepler residual |E − e sin E − M| is below 1e-10 for e ∈ [0, 0.8] (Float64). Five
    iterations reach round-off up to e = 0.8; e = 0.9 needs six;
  - `CircularOrbit` returns `S.solar_declination` and a distance factor of exactly 1;
  - `KeplerOrbit` with e = 0 gives δ = asin(sin(tilt) sin(g − g_equinox)) and a distance factor
    of 1;
  - the annual mean of the distance factor is 1/√(1 − e²);
  - the solstices and equinoxes fall at λ = 0, 90°, 180°, 270°;
  - perihelion falls where ϖ puts it;
  - a short `PrimitiveWetModel` run with `KeplerOrbit` stays finite, and the stored distance
    factor varies between 1/(1 + e)² and 1/(1 − e)²;
  - the default model is unchanged (`CircularOrbit`, distance factor 1).
- **Existing tests still pass:** zenith (72), shortwave (3267), radiation (82), variables
  (13), all parameterizations (63), greenhouse gases (12), longwave (10).
- **Insolation against climlab** (offline, reported in the PR rather than as a unit test):
  - setup: daily-mean TOA insolation, 5° latitudes × 365 days, for the Berger & Loutre (1991)
    orbits in climlab's table at 0, 6, 116 and 127 ka, with S₀ = 1365.2 and the equinox on
    21 March;
  - agreement: 0.36–0.38 W/m² RMS and 1.6–1.9 W/m² at most. The remainder is the 365- vs
    365.2422-day calendar; the best calendar offset is 0;
  - June 21 at 65°N:

    | Orbit | This PR | climlab |
    |---|---|---|
    | 0 ka | 478.9 | 478.9 |
    | 6 ka | 506.5 | 506.5 |
    | 116 ka | 439.9 | 439.9 |
    | 127 ka | 547.5 | 547.6 |

  - global annual means agree to 0.02 W/m².

## Documentation changes

An "Orbit" section in `radiation.md`: `CircularOrbit` / `KeplerOrbit`, the PMIP equinox
convention, the calendar-effect caveat, and an example setting a 126 ka-like orbit.

## Known limitations

- The orbit is fixed within a run (time-slice experiments only).
- The equation of time (`SolarTimeCorrection`) stays Earth's; it isn't derived from the orbit yet.
- For very eccentric orbits (e ≳ 0.3) the perihelion peak can be far brighter than the mean.
  That's a stability question for the physics at those fluxes, not for the orbit code.

## Future work

- `BergerOrbit(date)`: e, tilt and ϖ from Berger (1978) for any date in the last few Myr.
- Secular (transient) orbits.
- The equation of time from the orbit.
