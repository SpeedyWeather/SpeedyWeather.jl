# Output on pressure layers

> Status: **in review** ([#1257](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/1257)).
> The vertical interpolation from the model's layers onto pressure layers is merged
> ([#1256](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/1256)); this plan covers
> wiring it into the output writers behind a `layers` keyword argument, its documentation
> and tests.

Date of initial draft: 2026-09-16

Base revision: `af19cf29` (`main`, the merge of #1256)

## Originating prompt

The work started as a review of @henrikauestad's draft in
[#1256](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/1256):

> Henrik opened a pull request towards output on pressure levels. Can you review this pull
> request by commenting on it, and providing code suggestions to guide him? Guidelines are:
> The functionality should eventually be usable with all our output writers probably as a
> specific interpolator. The interpolation from sigma to pressure levels should work on both
> CPU/GPU with a kernel launched over the horizontal as the operation is parallelisable
> across columns. The interface should not depend on Variables but some externally allocated
> output field. Boundary conditions can be hardcoded for now but indicate where the code
> needs changes, e.g. to extrapolate the temperature dry-adiabatically down towards 1000hPa.
> Eventually this should also work with hybrid sigma-pressure coordinates so just indicate
> where changes would be needed in the future to accommodate for this.

and continued with the follow-up in
[#1257](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/1257):

> Can you branch off from main now and write a plan for the remaining tasks here into a new
> pull request? Don't make output on pressure levels yet the default but I liked the
> interface you designed with it being a keyword argument to the output writers. We'll need
> some test that writing output in a e.g. netcdf file with pressure levels actually works
> though.

> Subsurface default should be constant for all variables except temperature which should get
> a dry adiabatic descent. I prefer "pressure" over "plev" as it's more expressive [...] Do
> not export the vertical interpolation function or types for now, maybe we'll make this
> public interface at some point but not yet.

## Revision log

- **2026-09-15, review of #1256.** Design review of the draft interpolation, posted as a PR
  review with inline suggestions. Two bugs found (swapped interpolation weights, an
  out-of-bounds access at the top-most layer) plus the four design points from the prompt.
- **2026-09-16, #1256 numerics.** The reworked interpolation was contributed as a PR into
  Henrik's branch and merged with #1256.
- **2026-09-16, plan.** #1257 opened as a draft with this plan in its description, plus the
  documentation #1256 never had.
- **2026-09-16/17, implementation.** The `layers` keyword, the writer integration, tests and
  docs. Three decisions settled by the user: constant subsurface extrapolation except
  dry-adiabatic for temperature, the dimension named `pressure` rather than `plev`, and the
  interpolation function and its option types kept unexported.
- **2026-09-17, review of #1257.** Four changes from the review: *layers* instead of *levels*
  throughout (`PressureLayers`, `ModelLayers`, `interpolate_pressure_layers!`), matching
  `nlayers` and avoiding the ambiguity of ECMWF's "level"; `PressureLayers` made an immutable
  struct with concrete field types and a `::SpectralGrid` constructor that allocates, which
  removed its `initialize!`; `vertical_dimension` renamed to `vertical_dimension_name` as it
  returns the dimension's name; and the `extrapolation` field of `TemperatureOutput` made
  parametric, as is the `layers` field of the three writers.

## Problem description

SpeedyWeather integrates on terrain-following vertical coordinates, `SigmaCoordinates` or
`SigmaPressureCoordinates`. A model layer therefore sits at a different pressure in every
column, bending around the orography. Almost all analysis — comparison against reanalyses,
against observations, or against other models — wants the data on fixed pressure layers
instead, so this conversion otherwise has to be done by every user in post-processing, from
output that carries the surface pressure and the coordinate definition along.

## Summary of changes

### The interpolation (#1256, merged)

`SpeedyWeather/src/output/vertical_interpolation.jl`:

```julia
interpolate_pressure_layers!(
    out_field,          # OUTPUT: (horizontal, npressure), externally allocated
    in_field,           # INPUT: (horizontal, nlayers) on model layers
    surface_pressure,   # INPUT: (horizontal,) [Pa]
    p,                  # pressure layers [Pa]
    coordinates,        # model.geometry.vertical_coordinates
    interpolation,      # LinearInLogPressure() by default
    extrapolation,      # ConstantExtrapolation() by default
)
```

Four properties are worth recording, because they are what the review asked for:

- **No `Variables`, no `AbstractModel`.** The output field is allocated by the caller, which
  is what decouples the number of pressure layers from the number of model layers. The draft
  borrowed `vars.scratch.grid.a`, which both capped the layer count at `nlayers` and, because
  `getindex(::AbstractField, ::Colon, k...)` allocates a copy rather than a view, allocated a
  fresh array on every call.
- **The search runs in pressure, not in sigma.** `pressure(k, pₛ, coordinates)` dispatches
  over the coordinate type — `σ[k]·pₛ` for sigma, `A[k]·p_ref + B[k]·pₛ` for hybrid — so
  hybrid sigma-pressure coordinates are covered with no further changes. Interpolating in σ
  would have been wrong under hybrid coordinates, and log-in-σ only equals log-in-p because
  `pₛ` cancels for pure sigma.
- **One kernel over `(ij, k)`.** Columns and target layers are independent, so
  `launch!(arch, RingGridWorkOrder, size(out_field), ...)` covers CPU and GPU with a single
  code path. The pressure layers must live on the same architecture as the fields.
- **Extrapolation is a trait, not a branch.** `ConstantExtrapolation`,
  `DryAdiabaticExtrapolation(κ)` and `SubsurfaceMask(above_surface = ...)` extend
  `extrapolate_below`, which takes `pₛ` so that "below the lowest full model layer" and
  "below ground" stay distinguishable — a 1000 hPa layer under a lowest full layer at
  ~985 hPa is still above ground and wants extrapolating, not masking.

Two bugs in the draft were fixed on the way: the interpolation weights were swapped, mirroring
the result inside the layer (with `σ = [0.2, 0.6]`, values `[10, 20]` and a target at `0.3` it
returned `17.5` instead of `12.5`), and a pressure layer coinciding with the top-most model
layer fell into the interpolation branch with `k = 1` and indexed `σ[0]`.

### Output integration (#1257)

The vertical output grid is a property of the file being written, like the horizontal
`output_grid` and `interpolator` already are — not of the individual output variables. So it
attaches to the writer:

```julia
output = NetCDFOutput(spectral_grid, PrimitiveWet, layers = PressureLayers(spectral_grid, [850, 500, 200] .* 100))
```

- `AbstractOutputLayers` in `SpeedyWeather/src/output/writers/output_layers.jl`, with
  `ModelLayers()` (default, a no-op) and `PressureLayers(spectral_grid, pressure;
  interpolation)`, an immutable struct with concrete field types constructed from a
  `SpectralGrid` like the other model components. It owns the layers on the model's
  architecture and one scratch field on the model grid, both allocated at construction.
  Two buffers for the whole file, shared by every output variable regardless of how many
  are written.
- `layers` keyword on `NetCDFOutput`, `ZarrOutput` and `HEALPixOutput`, which also sizes their
  3D scratch field on the output grid via `get_nlayers(layers, spectral_grid)`.
- **One hook** in the generic `output!` (`writers/general.jl`), between selecting the time
  step and the copy to CPU, so every writer and every 3D atmospheric variable is covered at
  once:

  ```julia
  ori = has_step ? get_prognostic_step(ori, ts, output) : ori
  ori = interpolate_layers!(output_layers(output), ori, variable, simulation)   # no-op on ModelLayers
  raw = on_architecture(CPU(), ori)
  ```

- `initialize!(::OutputWriterCore, output, model)`, which every writer already calls, only
  syncs κ from the model's atmosphere; the buffers exist by then. `output_layers(output)`
  falls back to `ModelLayers()` so that writers without a `layers` field keep working.
- `vertical_dimension_name(output, variable)` gives variables on the model's atmospheric layers the
  writer's dimension while leaving custom ones (`soil_layer`, Terrarium's `soil_depth`)
  alone. The hardcoded sigma coordinate in each writer's `initialize!` became
  `define_vertical_coordinate!(dest, layers, model)`, backend-agnostic through the existing
  `define_coordinate!`.

### Why vertical first, then horizontal

The vertical interpolation deliberately runs on the model grid and the model's architecture,
*before* the copy to CPU and the horizontal interpolation. Three reasons, from the
[design discussion](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/1256#issuecomment-5679019257):

1. Each column is converted with **its own** surface pressure. Interpolating horizontally
   first would need `pₛ` on the output grid, which smooths the orography and makes the
   below-ground extrapolation inconsistent with the model's own surface.
2. It is usually less work: with fewer pressure layers than model layers, both the GPU→CPU
   copy and the horizontal interpolation shrink. (The reverse case — many pressure layers out
   of few model layers — makes this order do *more* horizontal work.)
3. It keeps the GPU-capable part on the GPU.

A fused 3D interpolator was considered and deferred. It would be *numerically identical* —
both stages are linear and the horizontal weights are layer-independent, so fusing changes
nothing — and it cannot be a single precomputed weight matrix, because the vertical weights
depend on `pₛ` and change every output step. It would save the intermediate buffer and the
per-layer kernel launches, at the cost of interpolating each source column roughly four times
(once per output point referencing it). Worth revisiting with a benchmark, behind the same
interface.

### Extrapolation defaults

Per variable, through an optional `extrapolation` field read with the same `hasproperty`
pattern that `output!` already uses for `transform` and `unscale`:

| variable | extrapolation below the lowest model layer |
| --- | --- |
| `TemperatureOutput` | `DryAdiabaticExtrapolation()`, `T(p) = T·(p/p_bottom)^κ` |
| everything else | `ConstantExtrapolation()` |

κ = R_dry/cₚ is taken from the model's atmosphere at `initialize!` (`sync_extrapolations!`),
so the adiabat matches the model rather than a hardcoded 2/7; the field is marked `[DERIVED]`
for that reason. The sync keeps the number format of the extrapolation so that the parametric
`extrapolation` field of an output variable does not change type; `extrapolate_below`
converts κ to the number format of the field it works on. This is the same adiabat the mean
sea-level pressure output uses to bring the lowest model layer down to the surface.

`SubsurfaceMask(above_surface = ...)` masks below the surface with the variable's
`missing_value`, but is not a default anywhere. Note that it masks *before* the horizontal
interpolation, and `anvil_average` is a bilinear blend where `NaN * 0 == NaN`, so one masked
source point poisons an output point even at zero weight: the mask dilates by the full
four-point stencil. A crisp mask needs masking on the output grid afterwards — see Future work.

## Testing and verification

`SpeedyWeather/test/output/vertical_interpolation.jl` (26 tests, from #1256):

- a field linear in `p` is reproduced exactly by `LinearInPressure`, one linear in `log(p)` by
  `LinearInLogPressure`, for both coordinate types and with a spatially varying `pₛ`
- pressure layers coinciding with model layers return the model layer values
- all three extrapolations, above the top, below the lowest layer but above ground, and below
  ground
- dimension mismatches throw; more pressure layers than model layers works

Both bug fixes were checked by **reverting them**: restoring the swapped weights fails the
exactness tests, and restoring the strict `<` together with a search starting at `k = 1`
reproduces `BoundsError: attempt to access 8-element Vector{Float32} at index [0]`.

`SpeedyWeather/test/output/pressure_layer_output.jl` (38 tests):

- `PressureLayers` construction (including that pressure and scratch land in the model's
  number format), `get_nlayers`, the sizing of `field3D` against the untouched `field3Dland`,
  and the monotonic/positive assertions, which now fire at construction
- extrapolation defaults per variable and the κ sync, including through a `SubsurfaceMask`
- a **NetCDF round trip**: a `PrimitiveWetModel` on a `FullGaussianGrid`, so the horizontal
  interpolation degenerates to a copy and the file can be compared against a direct
  `interpolate_pressure_layers!` on the same state (`rtol = 1e-5`). It checks that `pressure`
  exists with the right values and units, that `layer` is gone, that temp/u/v/vor/humid are on
  `pressure` while `st` stays on `soil_layer` and `mslp` stays 2D, and that the values are
  finite and decrease with height
- the default path is unchanged: a writer without `layers` still writes `layer` with the sigma
  values

Zarr (5 tests) and HEALPix (4 tests) pressure-layer tests live in their existing test files
rather than with the NetCDF ones, for the reason below.

### A local test-environment caveat

On Julia 1.13.0 on the development machine, the output tests abort in the GC
(`GC error (probable corruption)`, segfaults in `gc_mark_outrefs`). This is **not** caused by
these changes, verified by stashing them:

- pristine `test/output/netcdf_output.jl` aborted in 3 of 3 runs; with the changes in the tree
  it passed 1 of 2, i.e. nondeterministic either way and uncorrelated with the change
- pristine `test/output/zarr_output.jl` aborts as well
- `include`ing `netcdf_output.jl` and `zarr_output.jl` in the *same* process aborts reliably,
  on a pristine tree too

`ParallelTestRunner` gives each test file its own worker, so CI never exercises the
combination — which is why the Zarr and HEALPix pressure-layer tests were appended to their
own files instead of being grouped with the NetCDF ones. Possibly specific to this machine or
its depot; worth knowing about now that CI runs on 1.13 ([#1258](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/1258)).

## Documentation changes

- `docs/src/vertical_interpolation.md`, a new page in the Advanced section: the interface, the
  interpolation methods with their weight formulas, the extrapolation methods, and a runnable
  example interpolating temperature onto 850/500/200 hPa and down to 1000 hPa with subsurface
  masking. #1256 shipped without any documentation; this fills that gap.
- An [Output layers](../../src/output.md) section in `docs/src/output.md`, next to Output
  grid, covering the `layers` keyword and the extrapolation choices.

## Known limitations

- **Above the model top** every extrapolation method holds the top-most layer constant. Marked
  with a TODO in `extrapolate_above`; alternatives are masking with the missing value or an
  isothermal continuation.
- **Masking dilates.** `SubsurfaceMask` runs before the horizontal interpolation, so the masked
  region grows by the anvil stencil, see above.
- **One interpolation method per file.** `PressureLayers.interpolation` applies to every
  variable; only the extrapolation is per-variable.
- **Not exported.** `interpolate_pressure_layers!` and the interpolation/extrapolation types
  are internal for now, deliberately: only `ModelLayers` and `PressureLayers` are exported, as
  the handle for the `layers` keyword. The docs use qualified names accordingly.
- No GPU test for the output path yet; the interpolation itself is architecture-agnostic and
  tested on CPU.

## Future work

- A `mask_subsurface` option on `PressureLayers` that masks on the **output** grid after the
  horizontal interpolation, using a horizontally interpolated `pₛ` (one extra `Field2D` and one
  `interpolate!` per output step). That gives a crisp mask, keeps NaNs out of the interpolator
  entirely, and makes the mask edge agree with the `pres` field written into the same file.
- Geopotential height on pressure layers, which wants the hypsometric equation with the
  extrapolated temperature — another `extrapolate_below` method.
- The fused 3D interpolator, behind the existing interface, if a benchmark justifies it.
- Moving the horizontal output interpolation onto the GPU. `_interpolate!` is already a
  KernelAbstractions kernel and `AnvilLocator` is `Adapt`-ed; output is CPU-bound only because
  the writers allocate their scratch fields on a CPU output grid. That is independent of this
  work and a bigger win than the vertical kernel alone.
