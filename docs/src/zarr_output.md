
## Zarr Output

[Zarr](https://zarr.dev) is a chunked, compressed, cloud-friendly format for N-dimensional arrays. 

`ZarrOutput` is implemented as an **extension** that is only loaded once
[Zarr.jl](https://github.com/JuliaIO/Zarr.jl) is imported:

```@example zarr
using SpeedyWeather
using Zarr     # this loads SpeedyWeatherZarrExt and enables ZarrOutput

spectral_grid = SpectralGrid(truncation=32, nlayers=8)
output = ZarrOutput(spectral_grid, PrimitiveWet, interval=Hour(6))
model = PrimitiveWetModel(spectral_grid; output)
simulation = initialize!(model)
run!(simulation, period=Day(10), output=true)
nothing #hide
```

The constructor and option fields mirror [`NetCDFOutput`](@ref): `path`, `id`, `overwrite`,
`interval`, `variables`, `write_restart`, `write_parameters_txt`, `write_progress_txt` all
behave the same way, the run folder layout (`run_<id>_NNNN/`) is identical, and the same
`AbstractOutputVariable` types are used to declare which variables are written.
The on-disk layout differs:

- the Zarr store is a *directory* (`output.zarr/`), not a single file
- per-variable arrays live as subdirectories with `.zarray` and `.zattrs` metadata
- coordinates `lon`, `lat`, `layer`, `soil_layer`, `time` are stored as 1D arrays in
  the same group, tagged with the conventional `_ARRAY_DIMENSIONS` attribute so that Xarray-compatible readers can rebuild the dataset

Two extra options are specific to `ZarrOutput`:

| Option | Meaning |
|--------|---------|
| `time_chunk::Int` | Number of time steps per chunk along the time axis (default `1`). Larger values give bigger chunks and usually better compression at the cost of higher write latency. |
| `compressor` | Any `Zarr.Compressor` (e.g. `Zarr.BloscCompressor(clevel=3)`, `Zarr.ZlibCompressor()`); `nothing` (default) uses Blosc with the SpeedyWeather default compression level. |

```julia
using SpeedyWeather, Zarr

spectral_grid = SpectralGrid(truncation=32, nlayers=8)
output = ZarrOutput(spectral_grid, PrimitiveWet;
    interval = Hour(1),
    time_chunk = 24,                        # bundle one day per chunk on the time axis
    compressor = Zarr.BloscCompressor(clevel=5),
)
```

Reading back the data only needs Zarr.jl:

```@example zarr
using Zarr
g = Zarr.zopen(joinpath(output.run_path, output.filename))
g["time"][:]             # all stored hours since startdate
g["vor"][:, :, 1, :]     # vorticity, top layer, all time steps
nothing #hide
```

Custom output variables work exactly as with `NetCDFOutput`
(see [Output variables](@ref) and [Customizing netCDF output](@ref)): subtype
`AbstractOutputVariable`, implement `path(::MyOutputVariable, simulation)` to
return the `AbstractField` to write, and `add!(output, MyOutputVariable())` it to the
`ZarrOutput`. 

### Reading the store from Python with xarray

`ZarrOutput` writes the stores according to the conventions of `xarray` to ensure compability with it. Simulations can be easily opened in Python as in the following: 

```python
import xarray as xr
ds = xr.open_zarr("run_0001/output.zarr", consolidated=False)
print(ds)
# <xarray.Dataset>
# Dimensions:  (time: 41, layer: 8, lat: 32, lon: 64)
# Coordinates:
#   * time     (time) datetime64[ns]  2000-01-01 ... 2000-01-11
#   * layer    (layer) float32        0.06 0.19 ... 0.94
#   * lat      (lat) float64          85.76 ... -85.76
#   * lon      (lon) float64          0.0 ... 354.4
# Data variables:
#     vor      (time, layer, lat, lon) float32
#     u        (time, layer, lat, lon) float32
#     v        (time, layer, lat, lon) float32
#     temp     (time, layer, lat, lon) float32
#     humid    (time, layer, lat, lon) float32
#     mslp     (time, lat, lon) float32

ds["temp"].isel(time=-1, layer=0).plot()       # last step, top layer
ds["mslp"].mean(("lat", "lon")).plot.line(x="time")
```

`xarray` decodes the `time` axis to `datetime64` automatically (from the
CF-style `hours since <startdate>` units), so resampling and slicing by date
work without extra work:

```python
ds.sel(time=slice("2000-01-05", "2000-01-08"))["temp"].mean("time")
```

### Ensemble output

Several runs of an ensemble can be written into a *single* Zarr store along an
additional `ensemble` dimension. This is disabled by default and controlled by two
options:

| Option | Meaning |
|--------|---------|
| `ensemble_index::Int` | This writer's ensemble member index. `0` (default) disables ensemble output and reproduces the layout described above. A value `> 0` adds an `ensemble` dimension and makes this writer store its data into slot `ensemble_index`. Members are indexed `1..ensemble_size`. |
| `ensemble_size::Int` | Total number of ensemble members. It sizes the `ensemble` dimension and must be passed up front; it has to satisfy `ensemble_size ≥ ensemble_index`. |
| `ensemble_timeout::Int` | Seconds a member waits for member 1 to create the shared store before erroring (default `600`). |

The design targets **parallel ensemble members running as separate processes**, all
writing into one store. This is safe because the `ensemble` axis is chunked with size 1,
so member `e` only ever writes chunk files carrying its own index and no two members
touch the same chunk file. Each process constructs a `ZarrOutput` with the *same*
`path`, `id`, `run_number`, `filename` and `ensemble_size`, and a distinct `ensemble_index`,
so all members resolve to the same run folder and store:

```julia
using SpeedyWeather, Zarr

# in each process, `member` is this process' ensemble index (1..ensemble_size)
spectral_grid = SpectralGrid(truncation=32, nlayers=8)
output = ZarrOutput(spectral_grid, PrimitiveWet;
    ensemble_index = member,        # e.g. read from an environment variable / job array id
    ensemble_size = 10,
    interval = Hour(6),
)
model = PrimitiveWetModel(spectral_grid; output)
simulation = initialize!(model)
run!(simulation, period=Day(10), output=true)
```

The single shared artifacts (the group and coordinate metadata, and the `time` axis) are
written once by the *creator* — the member with `ensemble_index == 1`, which builds the
store schema and then signals readiness. The remaining *writer* members (`ensemble_index
> 1`) wait for that signal, open the existing store and write only their own ensemble
slice. This assumes that all members share the same time stepping so that one common
`time` axis is correct.

The resulting store gains an `ensemble` coordinate; the ensemble axis is the outermost
(slowest-varying) dimension, so an xarray-compatible reader reports e.g.

```python
import xarray as xr
ds = xr.open_zarr("run_0001/output.zarr", consolidated=False)
print(ds)
# Dimensions:  (ensemble: 10, time: 41, layer: 8, lat: 32, lon: 64)
# Data variables:
#     temp     (ensemble, time, layer, lat, lon) float32
ds["temp"].mean("ensemble")     # ensemble mean
```

Note that this ensemble layout is currently specific to `ZarrOutput`; the
[`NetCDFOutput`](@ref) does not support concurrent multi-process writes. 

## HEALPix Output

[HEALPix](https://healpix.sourceforge.io) is an equal-area discretization of the sphere. However with its non-rectangular coordinates, it's not straight forward to save HEALPix grids in standard NetCDF files. Instead `HEALPixOutput` writes
a Zarr store, keeping the horizontal dimension **flat** — one
unravelled vector of `npix = 12nside²` cells, exactly as a `Field` stores its data — rather
than interpolating onto a rectangular `lon`×`lat` grid the way [`NetCDFOutput`](@ref) and
`ZarrOutput` do. Use it with:

```@example healpix
using SpeedyWeather
using Zarr     # this loads SpeedyWeatherZarrExt and enables HEALPixOutput

spectral_grid = SpectralGrid(truncation=31, nlayers=8)
output = HEALPixOutput(spectral_grid, PrimitiveWet, nside=16, interval=Hour(6))
model = PrimitiveWetModel(spectral_grid; output)
simulation = initialize!(model)
run!(simulation, period=Day(1), output=true)
nothing #hide
```

The resolution is set with `nside` (as in most HEALPix implementation like cuHPX), 
equivalently `nlat_half = 2nside`, and defaults to the model grid's own `nlat_half` 
(rounded up to the nearest even number, which `HEALPixGrid`
requires). 

All the other file options of `ZarrOutput` (`path`, `id`, `overwrite`, `interval`,
`variables`, `time_chunk`, `compressor`, …) behave identically; `lon_chunk`/`lat_chunk` are
replaced by a single `cell_chunk` for the flat dimension.

### Store layout

| array | shape | `_ARRAY_DIMENSIONS` |
|---|---|---|
| 3D variable, e.g. `temp` | `(npix, nlayers, ntime)` | `["time", "layer", "cell"]` |
| 2D variable, e.g. `mslp` | `(npix, ntime)` | `["time", "cell"]` |
| `lat`, `lon`, `ring` | `(npix,)` | `["cell"]` |
| `cell` | `(npix,)` | `["cell"]` |

Alongside the data the store holds three flat vectors, one entry per cell, so that a consumer
can place every data point on the sphere without knowing anything about SpeedyWeather's grids:

- `lat`, `lon`: the coordinates of each cell in degrees north/east
- `ring`: the latitude ring index of each cell, 1 (north) to `4nside-1` (south)

```@example healpix
g = Zarr.zopen(joinpath(output.run_path, output.filename))
g["lat"][1:4], g["lon"][1:4], g["ring"][1:4]    # the 4 cells of the northernmost ring
```

### RING ordering and interoperability

Cells are written in standard **HEALPix RING order** of [cuHPX](https://github.com/NVlabs/cuHPX): 
cell `ij` (1-based, as in a `Field`) is RING pixel `ij-1` in the 0-based convention of healpy and
[cuHPX](https://github.com/NVlabs/cuHPX), so no reordering is needed on either side. The
store also carries `healpix_nside`, `healpix_npix` and `healpix_order` as global attributes:

```python
import xarray as xr, healpy as hp
ds = xr.open_zarr("run_0001/output_healpix.zarr", consolidated=True)
nside = ds.attrs["healpix_nside"]
hp.mollview(ds["temp"].isel(time=-1, layer=-1).values, nest=False)   # RING map, plots directly
```

These conventions should be compatabile with `healpy` and `cuHPX`.  Note that while `HEALPixOutput`
writes any even `nlat_half`, the NESTED and earth-2 flat layouts those tools convert to
require `nside` to be a power of two.

Our `lat` is latitude in degrees **from the equator** and `lon` is in `[0, 360)`, which is
exactly `healpy`'s `lonlat=True` convention. Mind that this is *not* healpy's default: with
`lonlat=False` (the default) healpy uses **colatitude in radians from the north pole**, so
convert with `theta = deg2rad(90 - lat)`, `phi = deg2rad(lon)` or use `lonlat=True`.

Mind the axis order when reading from Python: Zarr stores shape row-major while Julia is
column-major, so a variable written as `(cell, layer, time)` from Julia is seen as
`(time, layer, cell)` from Python — matching its `_ARRAY_DIMENSIONS` tag. A flat `(npix,)`
map to hand to `healpy` is therefore `arr.reshape(-1, npix)[-1]`, and the `hp.mollview` call
above works because `.isel(time=-1, layer=-1)` already reduces to the trailing `cell` axis.

### OctaHEALPix output

`HEALPixOutput` also writes the [`OctaHEALPixGrid`](@ref), a SpeedyWeather-specific member of
the HEALPix family (4 faces, `4nlat_half²` points, no equatorial belt). It is
written with **exactly the same store convention**: the same flat `cell` dimension, the same
flat `lat`/`lon`/`ring` coordinate vectors, the same array shapes and `_ARRAY_DIMENSIONS`.

```@example healpix
output = HEALPixOutput(spectral_grid, PrimitiveWet, output_grid=OctaHEALPixGrid(24))
output.field2D.grid
```

`healpy` and `cuHPX` implement only the standard HEALPix and **cannot** read an OctaHEALPix map.

By default the output grid **follows the model's grid** when the simulation already runs on
either supported grid, so both are written without interpolation.

### Skipping the interpolation

If the simulation already runs on the very HEALPix grid requested for output, no
interpolation is needed and none is set up — `output.interpolator` stays `nothing` and the
data is copied straight out of the model state:

```@example healpix
spectral_grid = SpectralGrid(truncation=31, nlayers=8, Grid=HEALPixGrid)
output = HEALPixOutput(spectral_grid, PrimitiveWet)
isnothing(output.interpolator)      # true, model and output share the grid
```

This saves both the per-output-step interpolation and the precomputation of the
interpolator's stencil indices and weights, and makes the written data bit-identical to the
model state (up to the number format and the `keepbits` mantissa rounding). Requesting a
different `nside` than the model's resolution brings the interpolator back.

````@docs; canonical=false
HEALPixOutput
````