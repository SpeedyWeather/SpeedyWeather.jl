# HEALPix Zarr output writer

> Status: **completed**. `HEALPixOutput` writes flat HEALPix-RING Zarr stores, skipping interpolation when the simulation already runs on the target HEALPix grid; unit tests pass and the opt-in healpy interoperability check passes at nside=8 (its optional cuHPX section is unrun, no CUDA hardware).

Date of initial draft: 2026-08-31

Base revision: 20467269399f40212fee413a545496226de999ca

## Originating prompt

> Write a new output writer that writes out healpix grids.
>
> * It should work with all grids in the simulation in which case the data has to interpolated to healpix, however if the simulation is already in healpix format, the interpolation should be skipped
> * The data should be saved in zarr format.
> * Try to reuse some of the utilities for zarr that we have for `ZarrOutput`
> * The data should be saved in a flat array (as is done in all our `Fields`
> * Additionally we want to save one flat array with the latitude coordinates of all data points, one for the longitude coordinates and one for the index of the ring
> * write unit test that check the funtionaility
> * write an additonal test - not part of the unit test - that check for compatability with cuHPX.

## Revision log

- 2026-08-31: initial draft and implementation in one pass.
- 2026-08-31: the out-of-suite check was reframed from cuHPX-first to **healpy**-first at the
  user's request (`write a test script with healpy instead`), and the directory renamed
  `test/cuhpx/` -> `healpix_compat/`. healpy is the reference HEALPix implementation,
  runs on CPU anywhere and can therefore actually be executed; cuHPX is CUDA-only and could
  not be run on the development machine. cuHPX coverage is retained as an optional
  auto-skipping section of the same script.
- 2026-08-31: the check was moved out of the `test/` tree entirely, to
  `SpeedyWeather/healpix_compat/` beside `benchmark/` (`don't actually add the
  interoperability check to the unit tests, not even as a check`). `test/runtests.jl` needed
  an exclusion regex while the directory sat under `test/`; with the move it needs nothing and
  is unmodified.

## Problem description

Machine-learning weather models and the wider "earth-2"/DestinE tooling ecosystem consume
gridded data on HEALPix in its native **flat, RING-ordered** layout: one unravelled vector of
`npix = 12·nside²` cells per (level, time), plus per-cell coordinate vectors. Both existing
gridded writers ([`NetCDFOutput`](../../../SpeedyWeather/src/output/writers/netcdf_output.jl)
and `ZarrOutput`) interpolate onto a **full** (rectangular `lon`×`lat`) grid and write a
`(nlon, nlat, …)` array. There is no way to get HEALPix data out of SpeedyWeather without a
lossy round trip through a regular lat/lon grid, and no way to avoid interpolation entirely
when the simulation itself already runs on `HEALPixGrid`.

## Background

- `RingGrids.HEALPixGrid` is the standard HEALPix discretisation with `nlat_half = 2·nside`,
  `npoints = 3·nlat_half² = 12·nside²` and `nlat = 4·nside − 1` rings. Its ring order and
  within-ring longitude offsets (`get_latd`, `get_lond_per_ring`, `each_index_in_ring!`)
  reproduce the canonical HEALPix RING scheme of Górski et al. (2005) exactly, including the
  alternating half-pixel offset `s ∈ {1,2}` in the equatorial belt. This makes the flat index
  `ij` (1-based) equal to the HEALPix RING pixel index (0-based) plus one, so no reordering is
  needed on write and the store is directly consumable by healpy/cuHPX.
- A `Field`'s `parent` is already the flat `(npoints, nlayers…)` array, so the flat layout is
  the *natural* one: `ZarrOutput` currently relies on Zarr.jl reshaping that flat buffer into
  its `(nlon, nlat, …)` array. Writing flat removes that implicit reshape.
- `RingGrids.interpolate!` already short-circuits to `copyto!` when the source and destination
  fields match, but `RingGrids.interpolator(out, in)` still precomputes a full `AnvilLocator`
  (four stencil indices and three weight vectors per output point). Skipping interpolation must
  therefore happen at *construction* time, not just at write time.

## Summary of changes

### New public type

`HEALPixOutput` (`SpeedyWeather/src/output/writers/healpix_output.jl`, implementation in
`SpeedyWeather/ext/SpeedyWeatherZarrExt/healpix_output.jl`) — a Zarr-backed output writer whose
output grid is a `HEALPixGrid` or an `OctaHEALPixGrid` and whose arrays keep the flat
horizontal dimension:

| array | shape (Julia, column-major) | `_ARRAY_DIMENSIONS` (row-major, xarray order) |
|---|---|---|
| 3D variable | `(npix, nz, ntime)` | `["time", "layer", "cell"]` |
| 2D variable | `(npix, ntime)` | `["time", "cell"]` |
| `lat`, `lon`, `ring` | `(npix,)` | `["cell"]` |
| `cell` | `(npix,)` | `["cell"]` |

`lat`/`lon` are the per-cell coordinates in degrees, `ring` the 1-based ring index of every
cell (north to south) as taken from `grid.whichring`. `cell` is the flat index `1:npix` and
exists so `lat`/`lon`/`ring` are well-formed non-dimension coordinates for xarray.

Resolution is set with `nside` (or `nlat_half = 2·nside`), defaulting to the simulation's own
`nlat_half` rounded up to the nearest even number (`HEALPixGrid` requires even `nlat_half`;
`OctaHEALPixGrid` takes any).

### Two grids, one store convention

Both `HEALPixGrid` and `OctaHEALPixGrid` are equal-area and ring-ordered north to south, and
both are written with an identical store layout — same flat `cell` dimension, same flat
`lat`/`lon`/`ring` vectors, same shapes and `_ARRAY_DIMENSIONS`. They differ only in the
advertised metadata, resolved by `healpix_store_attributes(grid)`:

Both carry the *same* key set, so a reader needs no special-casing: `grid`, `npix`,
`nlat_half`, `nrings`, `cell_ordering`, `equal_area`, plus `healpix_nside`, `healpix_npix`,
`healpix_nfaces` and `healpix_order = "RING"`. `nside` is generalized to **the side length of
one square face** so that the invariant `healpix_npix == healpix_nfaces * healpix_nside^2`
holds for both — 12 faces and `nside = nlat_half ÷ 2` for `HEALPixGrid` (the usual
`npix = 12nside²`), 4 faces and `nside = nlat_half` for `OctaHEALPixGrid` (`npix = 4nside²`).
`healpix_order` is `"RING"` for both, since cells run ring by ring north to south either way.

`healpix_nfaces` is therefore the discriminator, and it carries real weight: healpy and cuHPX
implement only the 12-face tessellation, so a reader assuming 12 faces would compute
`12nside²` and silently misread an OctaHEALPix store (at `nlat_half = 8`: 256 actual points
vs. 192 assumed). OctaHEALPix stores additionally carry an explicit `note`, and
`check_healpy.py` aborts on `healpix_nfaces != 12` rather than on a missing `healpix_nside`,
which no longer discriminates.

Grid selection (`healpix_output_grid`) is, in order of precedence: an explicit `output_grid`;
`nside`, which always means `HEALPixGrid` since `nside` is a HEALPix-only parameter; otherwise
the *model's own* grid type when it already runs on one of the two, and `HEALPixGrid` for
everything else. That default is what lets both HEALPix and OctaHEALPix simulations skip the
interpolation without the user asking for it.

### Interpolation skip

`HEALPixOutput`'s constructor compares the output grid against the (CPU) model grid with
`RingGrids.grids_match`. On a match it stores `interpolator = nothing` and the writer copies
straight from the model field into the output field.

To make that reachable from the shared output path, the five direct
`RingGrids.interpolate!(dest, src, output.interpolator)` call sites are replaced by a new hook
`interpolate_output!(output, dest, src)` (`src/output/writers/general.jl`), which dispatches to
`copyto!` when `output.interpolator === nothing`. Because the interpolator's type is a struct
type parameter, the `isnothing` branch is resolved at compile time.

### Zarr utility reuse

`ZarrOutput` gains a supertype, `AbstractZarrOutput <: AbstractOutput`, shared with
`HEALPixOutput`. The extension's time-chunk buffering machinery (`time_buffered`,
`write_array!`, `write_slice!`, `flush_time_chunk!`, `flush_partial_time_chunks!`), the time
axis writer `output!(::AbstractZarrOutput, ::DateTime)`, `Base.close`, `write_coordinate!`,
`define_coordinate!`, `get_dimension_length` and `resolve_compressor` are all shared verbatim.
Only two hooks differ:

- `zarr_write_indices(output, i, variable)` — `(:, :, :, i)`-style indices with the optional
  ensemble slot for `ZarrOutput`, flat `(:, :, i)`-style indices for `HEALPixOutput` (new
  `get_flat_indices` in `src/output/variables/output_variables.jl`).
- `zarr_ensemble_index(output)` — `0` for writers without ensemble support.

Ensemble output is deliberately **not** offered by `HEALPixOutput` (see Known limitations).

The single-file extension `ext/SpeedyWeatherZarrExt.jl` becomes a directory
`ext/SpeedyWeatherZarrExt/` with `SpeedyWeatherZarrExt.jl`, `shared.jl`, `zarr_output.jl` and
`healpix_output.jl`, matching the layout of `SpeedyWeatherTerrariumExt`.

## Testing and verification

Unit tests in `SpeedyWeather/test/output/healpix_output.jl`:

1. **Type and defaults** — construction, `nside`/`nlat_half` resolution selection, odd
   `nlat_half` rounded up to even.
2. **Interpolating path** — `ShallowWaterModel` on the default octahedral Gaussian grid;
   checks store layout, flat shapes, coordinate arrays, time axis and finiteness.
3. **Skip-interpolation path** — a model constructed on `HEALPixGrid` at the output resolution;
   asserts `output.interpolator === nothing` and that the written data is *bit-identical* to
   the model field (`Float32` output of a `Float32` model), which is only true if no
   interpolation ran.
4. **Coordinates are canonical HEALPix RING** — the stored `lat`/`lon`/`ring` are compared
   against an independent in-test implementation of the Górski et al. (2005) `z`/`φ` pixel
   centre formulas (the canonical `z = 1 − i²/3nside²` form rather than RingGrids' `acosd`
   form), plus the RING invariants: `4nside−1` rings, cells grouped ring by ring, ring
   lengths `4i`/`4nside`/`4(4nside−i)`, monotonically decreasing ring latitude, no pixel
   at either pole, and the self-describing `healpix_*` global attributes.
5. **Time chunking** — `time_chunk > 1` buffering flushes a trailing partial chunk on
   `finalize!`, matching an unbuffered reference run.
6. **PrimitiveWet with soil layers** — the `soil_layer` vertical dimension and 2D-from-3D
   variables (`mslp`, `u10`, `tsurf`) that carry their own `output!` methods.

The interoperability check is **not** part of the unit tests, in any form: it needs a Python
environment with healpy. It lives in `SpeedyWeather/healpix_compat/` — *outside* the `test/`
tree, next to `benchmark/` — so the suite's `find_tests` discovery never sees it and
`test/runtests.jl` is left completely untouched by this work. It is a Julia writer script
plus a Python checker, `check_healpy.py`, verifying against **healpy**, the reference
HEALPix implementation:

- our `lat`/`lon`/`ring` equal `healpy.pix2ang` and `healpy.pix2ring` pixel for pixel;
- the inverse closes — `healpy.ang2pix` on our own coordinates returns each cell's own
  index, and the centres coincide as unit vectors. This is the check that would break under
  any off-by-one, flipped hemisphere or rotated ring;
- `healpy.reorder` permutes the written *data map* RING↔NEST, matching the `nest2ring`
  index map, so the data and not just the coordinates is RING-ordered;
- healpy consumes it as a map: `ud_grade`, `get_interp_val`, `anafast`,
  `get_all_neighbours` adjacency, and the equal-area pixel size.

An optional final section exercises cuHPX (`ring2nest` vs `healpy.reorder`, and
`ring2flat`/`flat2ring` round-trips over all 8 conventions) and is reported as SKIPPED when
cuHPX or a CUDA device is unavailable; `--strict` turns skips into failures for a GPU CI job.
cuHPX is GPU-only and its conversions are index permutations of exactly the RING convention
verified above, so healpy conformance is what establishes cuHPX compatibility.

**Verified:** checks 1-5 pass at nside=8 against healpy 1.20.0 (22 checks, exit 0). The
cuHPX section has *not* been executed — no CUDA hardware was available.

Two bugs were found only by actually running this, both now fixed: the writer script's
`mktempdir` needed `cleanup=false` (Julia deleted the store on exit, before the checker could
read it), and the checker assumed Julia's axis order — Zarr stores row-major and Zarr.jl
reverses at the metadata boundary, so the flat `cell` dimension is the **last** axis from
Python and the first from Julia. That asymmetry is now documented in `other_output.md`.

## Documentation changes

- `docs/src/other_output.md` gains a "HEALPix output" section documenting `HEALPixOutput`, the
  store layout, the RING-ordering guarantee (and how it is verified), the Julia/Python axis-order
  asymmetry, and the interpolation-skip behaviour.
- `CHANGELOG.md` entry under `## Unreleased`.
- `SpeedyWeather/Project.toml` version bumped `0.22.1+DEV` → `0.23.0-DEV` (new public type).

## Known limitations

- **No ensemble output.** `ZarrOutput`'s ensemble machinery (shared store, readiness marker,
  member validation) is not wired up for `HEALPixOutput`; `zarr_ensemble_index` returns 0.
  Adding it later is mechanical — the ensemble axis would append to the flat shape exactly as
  it does for the rectangular one.
- **`OctaHEALPixGrid` stores are not portable.** They are written on request, with the same
  store convention *and* the same metadata key names, but healpy and cuHPX cannot read them.
  Consumers must branch on `healpix_nfaces` (or `grid`); `check_healpy.py` does exactly that
  and aborts with an explanation rather than cascading through meaningless failures.
- **Even `nlat_half` only for `HEALPixGrid`**, inherited from `RingGrids.nside_healpix`; odd
  requests are rounded up with an `@info`. `OctaHEALPixGrid` accepts any positive `nlat_half`.
- The interpolation-skip check is on grid *type and resolution*, not architecture: a GPU
  simulation on `HEALPixGrid` still pays the `on_architecture(CPU(), …)` copy per output step,
  as every writer does.

## Future work

- Ensemble support, mirroring `ZarrOutput`.
- Optional NESTED-order output (a permutation of the flat index at write time) for consumers
  that prefer it over RING.
- A matching reader/`animate` path for flat HEALPix stores in the GeoMakie extension.
