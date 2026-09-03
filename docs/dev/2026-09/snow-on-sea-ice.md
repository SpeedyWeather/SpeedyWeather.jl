# Snow on sea ice

> Status: **planned**. Awaiting human sign-off before implementation.

Date of initial draft: 2026-09-03

Base revision: e7c569df

## Originating prompt

> Our SnowModel accumulates snow over land but not over the ocean where it would melt in the
> water. But at the moment we don't account for snow landing on sea ice where it should have an
> effect on albedo and flux insulation. When the sea ice melts the snow should disappear. Can you
> make changes to the code to implement this with a flag in the snow model to enable/disable this
> effect?

Design decisions taken from a follow-up clarification with the user:

- Snow depth over sea ice is stored **per unit ice area**, not as a grid-cell mean.
- Snow disappears both **on sea ice retreat** and by **melting when the surface air temperature
  exceeds the melting threshold**.
- All three couplings are wired in: **albedo**, **sensible heat flux**, **humidity flux**.

## Revision log

- 2026-09-03: initial draft.

## Problem description

`SnowModel` (`SpeedyWeather/src/parameterizations/land/snow.jl`) only accumulates snow where
`land_fraction[ij] > 0`; its kernel is guarded by that condition and its prognostic variable
`snow_depth` lives in the `:land` namespace. Over the ocean, falling snow is implicitly assumed to
melt on contact with the water, which is correct for open ocean but wrong over sea ice: in reality
snow accumulates on the ice, where it substantially raises the surface albedo and insulates the
surface turbulent fluxes. When the ice melts away, the snow it carried must go with it.

The model already has all the surrounding machinery:

- `ThermodynamicSeaIce` (`parameterizations/sea_ice.jl`) carries a prognostic sea ice concentration
  `sea_ice_concentration` (ℵ) in the `:ocean` namespace.
- `OceanSeaIceAlbedo` (`parameterizations/albedo.jl`) already blends ocean and ice albedo by ℵ, and
  `LandSnowAlbedo` already has snow-cover schemes (`LinearSnowCover`, `SaturatingSnowCover`) that
  map snow depth to a cover fraction.
- `SurfaceOceanHeatFlux` / `SurfaceOceanHumidityFlux` already damp their fluxes by a
  `sea_ice_insulation` factor, and their land counterparts already damp by a
  `snow_insulation_depth` factor — the exact form we want to reuse over ice.

What is missing is a snow reservoir over the ice and the wiring of that reservoir into these three
consumers.

## Background

The land snow budget in `land_snow_kernel!` is a bucket in equivalent liquid water depth [m]:
snow falls at `snow_rate` [m/s], melts at a rate limited by the energy available from the top soil
layer being above `melting_threshold`, is floored at 0 and capped at `snow_depth_cap`. The excess
melt that could not be drawn from the bucket is folded into `snow_melt_rate` for the soil moisture
model.

Over sea ice there is no prognostic ice surface temperature and no soil layer, so the land energy
argument does not carry over. The available thermal driver is
`vars.parameterizations.surface_air_temperature`, which is what this plan uses. This is
deliberately a simple degree-day-style melt rather than a closed surface energy budget — see
*Known limitations*.

Note also that the two snow sinks act on different quantities. Melting removes snow mass. Ice
retreat does not melt the snow *per unit ice area* — it removes the area that was carrying it.
Because we store depth per ice area, retreat is represented by rescaling the reservoir so that the
total snow mass `ℵ·S` is reduced in proportion to the lost ice area, which for a per-area depth `S`
means `S` is unchanged by retreat alone but the mass it represents drops with ℵ automatically. To
make the snow actually *disappear* when the ice does — as the prompt requires — we additionally
zero the reservoir once ℵ falls below a small threshold. See *Summary of changes*, item 3.

## Summary of changes

### 1. Flag and parameters on `SnowModel`

Add to `SnowModel{NF}` in `parameterizations/land/snow.jl`:

```julia
"[OPTION] Accumulate snow on sea ice, affecting albedo and surface fluxes"
sea_ice_snow::Bool = false

"[OPTION] Sea ice concentration below which snow on sea ice is removed entirely [1]"
@param sea_ice_threshold::NF = 0.01 (bounds = 0 .. 1,)
```

`sea_ice_snow` defaults to `false` so existing setups are bit-for-bit unchanged. It is a plain
`Bool` field, not a `@param`, since it is a structural switch rather than a calibratable number.
`snow_depth_cap` and `melting_threshold` are shared with the land budget.

### 2. New prognostic variable

Extend `variables(snow::SnowModel)` to conditionally declare

```julia
PrognosticVariable(:snow_depth, Grid2D(), namespace = :ocean, units = "m",
                   desc = "Snow depth on sea ice in equivalent liquid water height")
```

only when `snow.sea_ice_snow` is `true`, so the allocation is skipped when the feature is off.
This makes `variables` depend on the flag's value, which is already the pattern used elsewhere for
optional components. It lives in the `:ocean` namespace (accessed as
`vars.prognostic.ocean.snow_depth`) because its lifecycle is tied to sea ice, not to land.
`initialize!(vars, snow, model)` sets it to 0 alongside the land field.

All consumers access it defensively via
`haskey(vars.prognostic.ocean, :snow_depth)`, falling back to zero — matching how
`sea_ice_concentration` is already accessed throughout. This means no consumer needs to know
about the flag, and the feature degrades cleanly for models without sea ice.

### 3. Sea ice snow budget

Add `sea_ice_snow_kernel!` and launch it from `timestep!(vars, snow::SnowModel, model)` when
`snow.sea_ice_snow` is true and both `sea_ice_concentration` and the ocean `snow_depth` are
present. Per grid point, guarded by `land_fraction[ij] < 1` (at least partially ocean):

1. **Accumulation.** `snow_rate[ij]` [m/s] falls on the ice. Because `S` is a per-ice-area depth,
   accumulation is *not* weighted by ℵ.
2. **Melt.** `δT = max(surface_air_temperature[ij] - melting_threshold, 0)`, giving a melt rate
   `melt_rate = melt_factor * δT` with a new
   `@param sea_ice_melt_factor::NF = 1e-7 (bounds = Nonnegative,)` [m/s/K] on `SnowModel`. As on
   land, melt is limited to the snow actually available, `min(S/Δt, melt_rate)`, so the bucket
   cannot go negative.
3. **Removal on ice retreat.** After the Euler-forward step,
   `S = ifelse(ℵ[ij] < sea_ice_threshold, 0, S)` — when the ice is gone, so is the snow. ℵ is read
   at the current prognostic step via `get_prognostic_step(..., model.sea_ice)`.
4. **Cap.** `S = min(S, snow_depth_cap)`, as on land.

The melt water is *not* routed anywhere: over ocean it joins the ocean, which the model does not
track as a freshwater budget. This is a genuine (small) mass-conservation gap, noted below.

Ordering: `timestep!(vars, land.snow, model)` runs inside the land timestep, and the sea ice
timestep runs separately in the main loop. The sea ice snow update reads ℵ from the current
prognostic step, so it sees the ice state from the previous sea-ice update. This one-step lag is
consistent with how other cross-component couplings in the model already work.

### 4. Albedo coupling

Give `OceanSeaIceAlbedo` the same snow fields `LandSnowAlbedo` has:

```julia
@param albedo_snow::NF = 0.4 (bounds = 0 .. 1,)
@param snow_depth_scale::NF = 0.05 (bounds = Positive,)
@param snow_cover::Scheme = SaturatingSnowCover() (group = :snow_cover,)
```

and in `albedo!` blend the snow contribution over the ice fraction only:

```julia
albedo[ij] = albedo_ocean + ℵ * (albedo_ice - albedo_ocean)
# snow sits on the ice, so weight its contribution by the ice fraction
albedo[ij] += ℵ * snow_cover * albedo_snow
```

with `snow_cover = snow_cover_scheme(snow_depth_ocean[ij], snow_depth_scale)`. The result is
clamped to 1 so that `albedo_ice + albedo_snow > 1` cannot produce an unphysical albedo. This adds
a type parameter `Scheme` to `OceanSeaIceAlbedo`, which is a breaking change to that struct's type
signature (keyword construction is unaffected).

### 5. Surface flux coupling

Add `snow_insulation_depth::NF = 0.05` to `SurfaceOceanHeatFlux` and `SurfaceOceanHumidityFlux`,
and apply it after the existing sea ice insulation:

```julia
flux_ocean /= 1 + sea_ice_concentration / heat_flux.sea_ice_insulation
# snow on sea ice insulates further, weighted by the ice fraction it sits on
flux_ocean /= 1 + sea_ice_concentration * snow_depth / heat_flux.snow_insulation_depth
```

The `sea_ice_concentration` weight is what makes this vanish over open ocean. The two divisions
compose the same way the land scheme composes its own damping factors.

## Testing and verification

New file `SpeedyWeather/test/parameterizations/snow_on_sea_ice.jl`, kept short per the repo
convention, covering:

1. **Flag off is inert.** With `sea_ice_snow = false`, `vars.prognostic.ocean` has no
   `snow_depth` key, and a short `PrimitiveWetModel` run reproduces the unmodified result exactly.
2. **Accumulation.** With the flag on, prescribed ℵ = 1 and a prescribed `snow_rate`, snow depth
   over ice grows by `snow_rate * Δt` per step and stays 0 over land-free open ocean.
3. **Melt.** With surface air temperature above `melting_threshold` and no snowfall, snow depth
   decreases monotonically and is floored at 0, never negative.
4. **Disappearance on ice melt.** Starting from nonzero snow depth, setting ℵ to 0 zeroes the snow
   depth on the next timestep.
5. **Cap.** Snow depth never exceeds `snow_depth_cap`.
6. **Albedo response.** `OceanSeaIceAlbedo` returns a strictly larger albedo with snow than
   without at ℵ = 1, equal albedo at ℵ = 0, and stays within [0, 1] even for
   `albedo_ice + albedo_snow > 1`.
7. **Flux insulation.** `SurfaceOceanHeatFlux` and `SurfaceOceanHumidityFlux` magnitudes decrease
   with snow depth at ℵ = 1 and are unchanged at ℵ = 0.

Run with `--check-bounds=yes` per the repo convention. A short `PrimitiveWetModel` run with the
flag enabled will also be checked for stability (no NaNs) on both `Float32` and `Float64`.

GPU correctness is inherited from the kernel form — the new kernel uses only `ifelse`/`min`/`max`
and no branching on data — but is not separately tested here.

## Documentation changes

- Docstrings on all new fields, following the existing `"[OPTION] ..."` convention.
- The land/snow section of the documentation gains a short subsection on snow over sea ice,
  documenting the flag, the budget, and the three couplings.
- `CHANGELOG.md` entry under `## Unreleased`.
- Version bump for `SpeedyWeather` (new public option and a changed struct signature on
  `OceanSeaIceAlbedo`, so a minor bump with `-DEV`).

## Known limitations

- **No surface energy budget over ice.** Melt is a tuned degree-day function of surface air
  temperature, not a closed energy balance. The `sea_ice_melt_factor` default of 1e-7 m/s/K is a
  plausible order of magnitude, not a calibrated value.
- **Melt water is not conserved.** Snow melting on sea ice is dropped rather than added to an
  ocean freshwater budget, which the model does not currently track.
- **Snow does not insulate the ice itself.** The insulation acts on the atmosphere-surface fluxes
  only; it does not slow the `ThermodynamicSeaIce` growth/melt, so snow cannot yet preserve ice
  through a warm season.
- **No blowing snow, no snow-ice formation, no flooding** of thick snow depressing the ice below
  sea level.
- **One-step lag** between the sea ice update and the snow-on-ice update, as described above.

## Future work

- Couple snow insulation into the `ThermodynamicSeaIce` budget so snow cover feeds back on ice
  survival.
- Replace the degree-day melt with a surface energy balance once a prognostic ice surface
  temperature exists.
- Track snow-on-ice mass in a conservation diagnostic alongside the land snow budget.
