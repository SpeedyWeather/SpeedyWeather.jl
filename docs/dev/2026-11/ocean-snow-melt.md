# Snow falling onto the slab ocean cools it by melting

> Status: **completed**. Implemented in this PR: the slab ocean is cooled by snowfall and by snow left over 100% ocean cells, and that snow is removed.

Date of initial draft: 2026-11-06

Base revision: c96ec4807fa6ca3d97e1cadb62f608837e8744dd

## Originating prompt

> Brandon responded, can you implement the melting + warming the meltwater in a PR? In that same PR can you fix that snow depth over 100% ocean cells should be set to 0, melting that too? Because it's possible to start with snow_depth = 0.2 globally for example and that snow melts over land but doesn't over the ocean ...

## Revision log

> (none)

## Problem description

Issue #1340: snow that falls into the ocean carries no latent heat of fusion. On land it takes
`Lᵢ` from the soil when it melts, but over the ocean it just disappears. Separately, `SnowModel`
only time steps cells with `land_fraction > 0`. Snow set over 100% ocean cells, e.g. from a
global initial `snow_depth = 0.2`, therefore never changes.

## Background

The estimate in #1340 shows that melting the snow (`Lᵢ`) plus warming the meltwater from 0˚C to the
SST (`c_w (SST - 0˚C)`) makes up 87–99% of the cooling that snow falling into the ocean causes. Warming the snow up to
0˚C is left out because it would need a snow temperature.

## Summary of changes

- `SlabOcean`: new sink `-ρ_w s (Lᵢ + c_w (SST - T₀))` [W/m²] in the SST tendency, with `ρ` the slab ocean `density`, `s` the
  snow rate [m/s] (`vars.parameterizations.snow_rate`, if any precipitation scheme defines it),
  `c_w = specific_heat_capacity`, `T₀ = 273.15 K` (new field `melting_temperature`).
- Over 100% ocean cells `s` additionally includes snow on the ground, `snow_depth / Δt_ocean`, if
  a land snow model exists. So this snow cools the ocean when it melts.
- `SnowModel` `filter!`: snow depth over 100% ocean cells is set to 0 after every time step.
  It is removed (an adjustment, not a tendency), consistent with the ocean having taken its heat
  during the same step.

## Testing and verification

- Unit tests in `test/parameterizations/ocean_sea_ice.jl`: SST tendency with/without snow matches
  the analytical sink, snow over 100% ocean cells is zero after a time step and cools the SST.

## Documentation changes

Docstrings of `SlabOcean` and `SnowModel`.

## Known limitations

- Warming the snow to 0˚C is ignored (≤ ~13% of the cooling, see #1340).
- Over partially-land cells, snowfall goes both into the land snow bucket and into the slab
  ocean (each per m² of its own fraction), as before for rain.
- Snow falling onto sea ice also cools the ocean underneath instead of the ice.
- Prescribed oceans just let the snow over ocean cells disappear.

## Future work

Convective snow should release `Lᵢ` into the atmosphere (the other half of #1340).
