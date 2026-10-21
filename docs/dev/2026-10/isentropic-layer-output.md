# Output on isentropic layers

> Status: **planned**. All decisions are settled (revisions 1 and 2); awaiting sign-off to
> implement.

Date of initial draft: 2026-10-21

Base revision: `adf0d61d` (fork `main`, branch `isentropic_layers`)

## Originating prompt

> For a new branch of the updated fork, make a plan for adding isentropic surfaces, surfaces
> of constant potential temperature, as output layers to be interpolated onto

## Revision log

- **2026-10-21, initial draft.** Branch `isentropic_layers` created from fork `main`
  (`adf0d61d`). Upstream `main` is one commit ahead (`d6e3286d`, #1330, an Enzyme extension
  fix that does not touch the output code), so the base is current for this work.
- **2026-10-21, separate clone.** To avoid interfering with the PV work, this work moved to
  its own folder `isentropic/` (next to `pv/`) with a separate clone of the fork. The
  branch `isentropic_layers` is there; the `pv/` clone is back on `pv`. The PV check script
  lives in `isentropic/` and only reads `pv/potential_vorticity_output.jl`.
- **2026-10-21, revision 1 (decisions).**

  > ln p and potential temperature should vary linearly, make sure that the interpolator
  > does that.
  > 1. Make the extrapolation routine consistent with the pressure level practice.
  > 2. Θ should be linear with ln p
  > 3. Is this really a problem on sigma or sigma-hybrid?
  > 4. Fine. 5. Fine. 6. Use dry potential temperature. 7. Yes, follow CLAUDE.md. 8. Yes.

  - The interpolator now assumes θ linear in ln p between two model layers. It solves for
    the isentrope's pressure p* and interpolates every variable linearly in ln p at p*.
    The `interpolation` option and the two θ interpolation types are removed.
  - Extrapolation now works as for pressure layers. Each variable's own `extrapolation`
    (constant by default, dry-adiabatic for temperature, optional `SubsurfaceMask`) is
    applied at p*. The layers' `extrapolation` field and `MissingValueExtrapolation` are
    removed.
  - Decision 3 is answered below. The search direction stays open.
  - Decisions 4–8 are settled as proposed.
- **2026-10-21, revision 2 (decisions 3 and 1a).**

  > First crossing from the top is OK. Search the web for conventions for 1a and list them

  > 1a: use the Lorenz condition, p* = pₛ below (and p* = p_top above), as proposed.

  - Decision 3: the first crossing from the top.
  - Conventions for underground isentropes, from a web search, are listed under 1a.
  - Decision 1a: the Lorenz condition for an unstable outermost layer pair.

## Problem description

Output writers can write 3D atmospheric variables on the model's layers (`ModelLayers`) or
on pressure layers (`PressureLayers`, #1256/#1257). Many dynamical analyses want them on
**isentropic layers**, surfaces of constant dry potential temperature θ = T (p₀/p)^κ,
instead: PV on 315/330/350 K, flow along isentropes, the tropopause as the 2-PVU contour
on an isentrope. Today this has to be done in post-processing from model-layer output that
carries T, pₛ and the coordinate definition along.

The goal is a third subtype of `AbstractOutputLayers`,

```julia
output = NetCDFOutput(spectral_grid, PrimitiveDry, layers = IsentropicLayers(spectral_grid, [300, 315, 330, 350]))
```

that writes every 3D atmospheric output variable on the given θ layers [K], through the
same hook as `PressureLayers`, on CPU and GPU.

## Background

How `PressureLayers` works (`SpeedyWeather/src/output/writers/output_layers.jl`,
`SpeedyWeather/src/output/vertical_interpolation.jl`):

- `layers` is a parametric field of `NetCDFOutput`, `ZarrOutput` and `HEALPixOutput`.
- The generic `output!` (`writers/general.jl`) calls
  `interpolate_layers!(output_layers(output), field, variable, simulation)` after selecting
  the time step and before the horizontal interpolation. The vertical interpolation runs on
  the model grid and the model's architecture.
- Everything else dispatches on the layers type: `get_nlayers` (sizes the writer's 3D
  scratch field), `vertical_dimension_name`, `define_vertical_coordinate!`, `show`.
- `interpolate_pressure_layers!` launches one kernel over `(ij, k_out)`. In each column it
  finds the bracketing model layers in pressure (`pressure(k, pₛ, coordinates)`, sigma or
  hybrid) and blends them with `interpolation_weight`.
- Outside the model layers, the variable's `extrapolation` decides, via
  `output_extrapolation(variable)`:
  - `ConstantExtrapolation` is the default;
  - `DryAdiabaticExtrapolation` is used for `TemperatureOutput`, with κ synced from the
    atmosphere by `sync_extrapolations!`;
  - `SubsurfaceMask` is optional.

  Below the lowest layer this is `extrapolate_below(ξ_bottom, p, p_bottom, pₛ, E)`. Above
  the top it is `extrapolate_above(ξ_top, p, p_top, E)`, which is constant for every method.

Because all hooks dispatch on the layers type, isentropic layers are **additive**: new
methods for a new type, no change to existing functions or to the writers. The PV user
script calls `interpolate_layers!(output_layers(output), pv, variable, simulation)`, so PV
gets isentropic output without changes.

What is different for θ:

1. **The vertical coordinate depends on the state.** p(k) needs only pₛ; θ(k) also needs
   the temperature. `vars.grid.temperature` is absolute temperature [K] with a step
   dimension, so it is taken at the same step as the field being written. κ and p₀ come
   from `model.atmosphere.κ` and `model.atmosphere.reference_pressure` (1e5 Pa). θ is the
   dry potential temperature, also in the wet model (decision 6).
2. **θ is not always monotonic in a column.** See decision 3 below.
3. **Isentropes intersect the ground as a rule, not an exception.** For example, 290 K is
   below ground in the tropics and well above ground at the poles in winter. The
   extrapolation methods take a pressure, so the isentrope's pressure p* is extrapolated
   first (see 3 below).

## Summary of changes

All in `SpeedyWeather/src/output/` unless noted.

### 1. `IsentropicLayers` (`writers/output_layers.jl`)

```julia
export ModelLayers, PressureLayers, IsentropicLayers

"""Potential temperature [K] of the layers that `IsentropicLayers` interpolates onto by default."""
const DEFAULT_ISENTROPIC_LAYERS = Float32[285, 300, 315, 330, 350]

struct IsentropicLayers{V, F} <: AbstractOutputLayers
    "[OPTION] dry potential temperature [K] of the layers to interpolate onto, monotonic"
    potential_temperature::V
    "[DERIVED] scratch field (npoints, nθ) on the model grid"
    scratch::F
end

IsentropicLayers(SG::SpectralGrid; potential_temperature = DEFAULT_ISENTROPIC_LAYERS)
IsentropicLayers(SG::SpectralGrid, potential_temperature::AbstractVector)
```

The constructor follows `PressureLayers`: monotonic and positive θ are asserted, and θ and
the scratch field are allocated on `SG.architecture` in `SG.NF`. There is no
`interpolation` option, because θ linear in ln p is fixed (decision 2), and no
`extrapolation` option, because extrapolation is per variable as for pressure layers
(decision 1).

Methods:

- `get_nlayers(layers::IsentropicLayers, ::SpectralGrid) = length(layers.potential_temperature)`
- `vertical_dimension_name(::IsentropicLayers) = "potential_temperature"`
- `define_vertical_coordinate!(dest, ::IsentropicLayers, model)`: values in K,
  `units = "K"`, `long_name = "potential temperature"`,
  `standard_name = "air_potential_temperature"`, `positive = "up"`
- `Base.show` like `PressureLayers` (θ in K)
- `interpolate_layers!(layers::IsentropicLayers, field, variable, simulation)`, see 4.

### 2. Interpolation: θ linear in ln p

Between two neighbouring model layers k−1 and k in a column, θ is assumed linear in ln p:

    θ(ln p) = θ₁ + (θ₂ − θ₁)·(ln p − ln p₁)/(ln p₂ − ln p₁)

Here θ₁, p₁ belong to layer k−1 and θ₂, p₂ to layer k. For a target θ between θ₁ and θ₂,
the isentrope's pressure follows in closed form:

    ln p* = ln p₁ + w·(ln p₂ − ln p₁),   w = (θ − θ₁)/(θ₂ − θ₁)

Every output variable is then interpolated linearly in ln p at p*, as `PressureLayers` does
by default:

    ξ(θ) = (1 − w)·ξ₁ + w·ξ₂,   with w = interpolation_weight(p*, p₁, p₂, LinearInLogPressure())

Both weights are the same number: the ln p weight at p* equals the θ weight. The kernel
computes w once, and p* only for the extrapolation. A field linear in ln p is reproduced
exactly, and θ itself, interpolated this way, returns the target θ.

### 3. Extrapolation, as for pressure layers

The variable's own extrapolation (`output_extrapolation(variable)`) is used through the
existing `extrapolate_below`/`extrapolate_above`. These take a pressure, so p* is
extrapolated with the same rule as inside the column (θ linear in ln p), using the slope of
the outermost layer pair:

- **Below the lowest layer** (θ < θ(nlayers), no crossing):
  `ln p* = ln p_bottom + (θ − θ_bottom)·(ln p_bottom − ln p_{nlayers−1})/(θ_bottom − θ_{nlayers−1})`,
  so p* > p_bottom in a stably stratified bottom pair. Then
  `extrapolate_below(ξ_bottom, p*, p_bottom, pₛ, E)`, which behaves as follows:
  - `ConstantExtrapolation` (default) → the lowest layer's value;
  - `DryAdiabaticExtrapolation` (temperature) → T_bottom·(p*/p_bottom)^κ;
  - `SubsurfaceMask` → `missing_value` where p* > pₛ (an isentrope underground), and its
    `above_surface` method between the lowest layer and the surface.
- **Above the top layer** (θ > θ(1), no crossing): p* from the top pair in the same way,
  and `extrapolate_above(ξ_top, p*, p_top, E)`. For every method this is currently the
  top layer's value, exactly as for pressure layers.
- **Unstable outermost pair** (θ_bottom ≥ θ_{nlayers−1}, or θ_top ≤ θ_2): no p* beyond the
  layer follows from the slope. The Lorenz condition applies: p* = pₛ below and
  p* = p_top above (decision 1a).

Because the variables' extrapolation is reused, `sync_extrapolations!` (κ, number format)
already covers it. Nothing new is needed at `initialize!`.

### 4. The kernel and the hook

`vertical_interpolation.jl`:

```julia
interpolate_isentropic_layers!(
    out_field,          # OUTPUT: (horizontal, nθ), externally allocated
    in_field,           # INPUT: (horizontal, nlayers) on model layers
    temperature,        # INPUT: (horizontal, nlayers) absolute temperature [K]
    surface_pressure,   # INPUT: (horizontal,) [Pa]
    θ,                  # isentropic layers [K]
    coordinates,        # model.geometry.vertical_coordinates
    κ, p₀,              # from model.atmosphere
    extrapolation = ConstantExtrapolation(),
)
```

Like `interpolate_pressure_layers!`, this has no `Variables` or model argument, checks
grids and sizes (`DimensionMismatch`), and launches one kernel over `(ij, k_out)` with
`RingGridWorkOrder`. The column function:

1. computes θ(k) = T[ij, k]·(p₀/pressure(k, pₛ, coordinates))^κ on the fly, for sigma and
   hybrid coordinates alike;
2. searches for the bracketing pair (decision 3), with θ between θ(k−1) and θ(k) in either
   order. If θ(k−1) = θ(k) = θ, it returns the value of layer k−1, avoiding 0/0;
3. blends `(1 − w)·ξ[k−1] + w·ξ[k]` with w as in 2;
4. without a crossing, extrapolates as in 3.

Cost: nθ·nlayers powers per column per variable. That is small compared with the
horizontal interpolation, so θ is recomputed per variable rather than cached per output
step.

`writers/output_layers.jl`:

```julia
function interpolate_layers!(layers::IsentropicLayers, field::AbstractField,
        variable::AbstractOutputVariable, simulation::AbstractSimulation)
    (is3D(variable) && !is_land(variable)) || return field
    (; model, variables) = simulation
    temperature = get_prognostic_step(variables.grid.temperature, model.time_stepping, layers)
    interpolate_isentropic_layers!(
        layers.scratch, field, temperature, variables.parameterizations.surface_pressure,
        layers.potential_temperature, model.geometry.vertical_coordinates,
        model.atmosphere.κ, model.atmosphere.reference_pressure,
        output_extrapolation(variable))
    return layers.scratch
end
```

The step: the hook has no `output` argument, so `layers` (an `AbstractModelComponent`)
selects the step. Neither writers nor layers have a `which_prognostic_step` specialisation,
so both fall back to step 1, the step the writer passes in `field`. A test pins this. The
hook's signature stays, so the PV script's call keeps working.

Only models with temperature (`PrimitiveDry`, `PrimitiveWet`) support isentropic layers.
Constructing a writer with `IsentropicLayers` for Barotropic or ShallowWater raises a clear
error at `initialize!`.

### 5. Other

- `SpeedyWeather/Project.toml`: `0.23.0+DEV` → `0.24.0-DEV`. This follows CLAUDE.md for a
  new public type (decision 7).
- `CHANGELOG.md`: one bullet under `## Unreleased`. There is no PR number because nothing is
  pushed, so it is a placeholder.
- Exports: `IsentropicLayers` only, like `PressureLayers`. `interpolate_isentropic_layers!`
  stays unexported and is used qualified in the docs.

## Testing and verification

New `SpeedyWeather/test/output/isentropic_layer_output.jl`, auto-discovered:

- **Kernel, exact cases** for `SigmaCoordinates` and `SigmaPressureCoordinates` with a
  horizontally varying pₛ. T on the layers is built from a θ profile that is linear in ln p
  in every layer, T = θ·(p/p₀)^κ. In Float64, `rtol = 1e-10`:
  - a field linear in ln p is reproduced exactly at every isentrope;
  - ln p itself, written as a field, gives ln p* from the closed form;
  - θ itself, written as a field, gives the target θ;
  - isentropes equal to a layer's θ return that layer's value.
- **Extrapolation, as for pressure layers**:
  - `ConstantExtrapolation` gives the outermost layer's value above and below;
  - `DryAdiabaticExtrapolation` gives T_bottom·(p*/p_bottom)^κ with p* from the bottom
    pair;
  - `SubsurfaceMask` masks exactly where p* > pₛ and not between the lowest layer and the
    surface;
  - the unstable outermost pair: p* = pₛ below and p* = p_top above (decision 1a). With an
    unstable bottom pair, `SubsurfaceMask` does not mask and `DryAdiabaticExtrapolation`
    gives T_bottom·(pₛ/p_bottom)^κ.
- **Non-monotonic column**: a column with two crossings returns the crossing chosen in
  decision 3.
- **Constructor**: number format and architecture of θ and scratch, `get_nlayers`,
  assertions for non-monotonic and non-positive θ, `DimensionMismatch` for wrong sizes.
- **Step**: `get_prognostic_step(var, time_stepping, layers)` and the writer's
  `get_prognostic_step(var, time_stepping, output)` select the same step.
- **NetCDF round trip**, as for pressure layers: `PrimitiveDryModel` on a
  `FullGaussianGrid`, so the horizontal interpolation is a copy. After a short run, the file
  is compared with a direct `interpolate_isentropic_layers!` on the same state
  (`rtol = 1e-5`). It also checks:
  - the `potential_temperature` coordinate (values, `units = "K"`);
  - that `layer` is absent;
  - that 2D variables (`mslp`) stay 2D;
  - that temperature uses the dry-adiabatic extrapolation.
- **Zarr and HEALPix**: one short write test each, appended to their existing test files,
  for the same Julia 1.13 GC reason as recorded in the pressure-layer plan.
- **GPU**: a CPU/GPU comparison of `interpolate_isentropic_layers!` added to
  `test/GPU/output.jl` next to the pressure-layer one. Run locally with Metal if the
  `MetalGPU` test environment instantiates; otherwise it is left to CI and reported as not
  run.

Run with `--check-bounds=yes`: the new test file, plus `vertical_interpolation.jl` and
`pressure_layer_output.jl` as regression tests (a few minutes each).

**PV check in `isentropic/` (decision 8).** It needs `isentropic/Project.toml` with
SpeedyWeather from this clone via `[sources]`, NCDatasets and CairoMakie, and leaves
`pv/potential_vorticity_output.jl` unchanged (only `include`d).
`isentropic/check_isentropic_potential_vorticity.jl` writes PV of the JW run (T42 L30,
Γ = 5 K/km, as in check 4) on 315, 330 and 350 K. It checks the following:

- PV is positive in the NH and negative in the SH, wherever the isentrope is above ground;
- the 2-PVU contour on 330 K lies in the midlatitudes;
- values match PV on model layers interpolated by hand in a few columns.

About 5 min. A figure is out of scope unless requested.

## Documentation changes

- `docs/src/vertical_interpolation.md`: a section "Isentropic layers". It covers:
  - θ on model layers;
  - θ linear in ln p and the closed form for p*;
  - the search for the crossing;
  - extrapolation through the variables' own methods at p*;
  - a runnable example that writes temperature and u on 300/330 K.
- `docs/src/output.md`, section Output layers: one paragraph and one code line for
  `IsentropicLayers`.

## Decisions

| # | Question | Decision |
|---|---|---|
| 1 | Extrapolation | As for pressure layers: the variable's own extrapolation, at p* extrapolated with θ linear in ln p |
| 1a | Unstable outermost pair: no p* beyond the layer | Lorenz condition: p* = pₛ below, p* = p_top above |
| 2 | Interpolation | θ linear in ln p; variables linear in ln p at p* |
| 3 | Several crossings in a column | first crossing from the top |
| 4 | Default layers | 285, 300, 315, 330, 350 K |
| 5 | Dimension name | `potential_temperature` |
| 6 | θ in the wet model | dry θ |
| 7 | Version bump | `0.24.0-DEV`, following CLAUDE.md |
| 8 | PV check in `isentropic/` | yes |

### 3. Is non-monotonic θ really a problem on sigma or hybrid coordinates?

The coordinate does not cause it. For both `SigmaCoordinates` (p = σ·pₛ) and
`SigmaPressureCoordinates` (p = A·p_ref + B·pₛ), pressure increases strictly with the layer
index in every column, which the pressure-layer search already relies on. θ is
non-monotonic only where the **state** has θ decreasing with height, that is, where it is
statically unstable (superadiabatic). That happens on any vertical coordinate.

Whether SpeedyWeather states have it, checked on the existing PV runs (T42 L30, θ from T
on the 100–850 hPa output layers, layers above the surface only):

| run | day 0 | day 6 | day 12 |
|---|---|---|---|
| JW, Γ = 5 K/km (`run_pv_jw_0003`) | 0 % of columns | 0 % | 0 % |
| JW, Γ = 9.8 K/km (`run_pv_jw_0002`) | 44 % of columns | (blew up on day 3) | |

So a stable dry run has no such columns between 850 and 100 hPa. The unstable initial
state of the default `JablonowskiTemperature` in PrimitiveDry has them in almost half the
columns. Not covered: layers below 850 hPa, and full-physics runs (PrimitiveWet with
surface fluxes). There, superadiabatic lowest layers over heated land are expected,
transiently, until convection and vertical diffusion remove them.

So the search direction matters only in unstable columns, but it has to be defined so that
the kernel gives a unique answer. Decided: the first crossing from the top.

### 1a. Unstable outermost layer pair

If θ is below every θ in the column and the bottom pair is unstable
(θ_bottom ≥ θ_{nlayers−1}), θ linear in ln p gives no p* below the lowest layer. The same
holds at the top. Decided (the Lorenz condition): p* = pₛ below (the isentrope is put at the surface; constant
and dry-adiabatic extrapolation are then well defined, and `SubsurfaceMask` does not mask),
and p* = p_top above (irrelevant today, as `extrapolate_above` is constant). The
alternative below is p* = ∞, so that `SubsurfaceMask` masks it as underground; but then
`DryAdiabaticExtrapolation` returns ∞ and would need a guard.

Conventions elsewhere for an isentrope below the lowest level or the ground (web search,
2026-10-21). None of them extrapolates θ below the lowest level; each assigns a fixed
value:

| convention | value of an underground isentrope | used by |
|---|---|---|
| missing | missing / NaN | MetPy `isentropic_interpolation` (no bracketing pair, or p* beyond the deepest pressure → NaN); GEMPAK `GDVINT` without `/L` |
| lowest level | value of the lowest model level, held constant | wrf-python `vinterp` with a θ target (extrapolation only for pressure/height targets, otherwise the lowest model level); NCEP reanalysis isentropic fields (held constant outside the domain) |
| Lorenz condition / massless layer | p(θ) = pₛ, z(θ) = zₛ, surface values, zero thickness | Lorenz (1955); Held & Schneider (1999, JAS 56); Koh & Plumb (2004); GEMPAK `GDVINT` with `/L` |

The proposal p* = pₛ is the Lorenz condition. With the default `ConstantExtrapolation` it
gives the lowest layer's value, which is the "lowest level" convention. With
`DryAdiabaticExtrapolation` it gives temperature on the adiabat at the surface. With
`SubsurfaceMask` it does not mask, so the "missing" convention is not reachable for
unstable bottom pairs.

## Known limitations

- **Masking dilates** by the horizontal interpolation's stencil (NaN·0 = NaN) with
  `SubsurfaceMask`, as for pressure layers.
- **Unstable columns**: only one crossing is written (decision 3).
- **No pressure on isentropes.** There is no 3D pressure output variable, so the file does
  not hold p*, or the Montgomery potential usually wanted with it. See Future work.
- θ is recomputed for every variable written. That is cheap, but redundant.

## Future work

- A 3D pressure output variable (p on model layers). It would give p* on isentropes for free
  through the same hook (exactly, since it is interpolated linearly in ln p), and the
  Montgomery streamfunction M = cₚT + Φ.
- Masking on the output grid instead of the model grid, as already listed for pressure
  layers.
- Caching θ (or the bracketing index and weight) once per output step, if a benchmark shows
  the recomputation matters.
