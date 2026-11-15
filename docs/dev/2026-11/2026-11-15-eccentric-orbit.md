# Eccentric orbits: Kepler orbit, Sun–planet distance and precession

> Status: **plan**. Nothing implemented yet; posted for discussion before any code.

Date of initial draft: 2026-11-15

Base revision: `main` as of 2026-11-15

## Originating prompt

> Full Milankovitch cycles: obliquity, eccentricity and precession, with the precession
> chosen as the month of perihelion.

(From Paleo Planet Sim, a worldbuilding app that runs SpeedyWeather for user-designed
planets. Prompted by teaching uses like Build Your Own Earth, Schultz et al. 2017, BAMS,
whose orbital experiments vary tilt, eccentricity and the orbit through the last glacial
cycle.)

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
  in Berger's convention). The exact values come from the PMIP tables when implementing.
- climlab (`climlab.solar.orbital`, `daily_insolation`; MIT licence): a reference
  implementation to test against.

## Proposed changes

1. **Planet parameters** (`planet.jl`):
   - `eccentricity::NF = 0` (bounds 0 ≤ e < 1);
   - `longitude_of_perihelion::NF` [°, from the vernal equinox]. The default is irrelevant at
     e = 0; Earth's present value would be a natural choice.

   The existing `equinox` field fixes the calendar (PMIP convention), and `length_of_year`
   is the orbital period.
2. **An orbit component** used by `SolarZenith`/`SolarZenithSeason` in place of the
   declination functor. Something like:

   ```julia
   abstract type AbstractOrbit <: AbstractModelComponent end
   struct KeplerOrbit <: AbstractOrbit end   # reads e, ϖ, tilt, equinox from planet
   # returns (δ, distance_factor = (a/r)^2) for a given orbit_time
   ```

   Steps:
   - the mean anomaly comes from orbit_time, anchored so that the true solar longitude
     λ = 0 at `planet.equinox`;
   - Kepler's equation E − e·sin E = M is solved by Newton's method (a fixed 5 iterations, so
     there is no data-dependent control flow for Reactant/Enzyme);
   - then ν from E, λ = ν + ϖ, sin δ = sin(tilt)·sin λ, and (a/r)² = ((1 + e·cos ν)/(1 − e²))².

   At e = 0 this reduces exactly to `SinSolarDeclination` with a distance factor of 1, so
   existing simulations are unchanged. That case would be tested bitwise.
3. **Distance factor in the shortwave.** The zenith step stores the scalar (a/r)² for the
   current time step: a new 0-dimensional parameterization variable, or a field on the zenith
   component, whichever fits your design. The three `solar_constant * cos_zenith` sites
   multiply by it. `cos_zenith` stays the geometric cosine, so path lengths and
   zenith-dependent albedos are unaffected.
4. **Equation of time** (optional, second step): for a general orbit
   EoT = (M − ν) − (λ_apparent − α) can replace the Earth-only Spencer series in
   `SolarTimeCorrection`. It is small, so it could be left as future work.

## What it does not do

- No secular change of the orbit during a run: the elements are constants, as in PMIP time
  slices. An `orbital_time` option evaluating Berger (1978) at a given date could follow later.
- No change to the calendar/month definitions. PMIP's "calendar effect" on monthly means at
  high e gets a docs note instead.

## Testing and verification

- **Unit tests:**
  - e = 0 is bitwise equal to `SinSolarDeclination` and gives a distance factor of 1;
  - Kepler residual below 1e-6 across e ∈ [0, 0.5];
  - the distance factor integrates over a year to 1/√(1 − e²);
  - the solstices and equinoxes fall where λ says they should.
- **Insolation against climlab** (or Berger & Loutre tables): daily-mean TOA insolation by
  latitude × day of year for present day, 6 ka, 116 ka and 126 ka, within 1 W/m². That
  includes the classic 65°N June value.
- **Model runs** (T31 L8, slab ocean):
  - present-day orbit vs circular, where we expect small, Southern-Hemisphere-summer-weighted
    differences;
  - 6 ka vs present, where we expect more North African summer rain;
  - a Mars-like e = 0.09 at zero tilt, where both hemispheres should warm together near
    perihelion.

## Documentation changes

A short "Orbit" section in `radiation.md`: the parameters, the PMIP equinox convention, the
calendar-effect caveat, and an example setting a 126 ka orbit.

## Known limitations

- The orbit is fixed within a run (time-slice experiments only).
- For very eccentric orbits (e ≳ 0.3) the perihelion peak can be far brighter than the mean.
  That's a stability question for the physics at those fluxes, not for the orbit code.

## Future work

- `BergerOrbit(date)`: e, tilt and ϖ from Berger (1978) for any date in the last few Myr.
- Secular (transient) orbits.
- The equation of time from the orbit (item 4) if it isn't done here.
