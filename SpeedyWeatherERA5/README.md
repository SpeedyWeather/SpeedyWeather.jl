# SpeedyWeatherERA5 (WIP)

Initialise a SpeedyWeather.jl `PrimitiveWetModel` from ERA5, Julia only. Naive initialisation, simply takes ERA5 variables and regrids to model's own sigma levels and grid. Developed against SpeedyWeather.jl main (0.23.0+DEV), Julia 1.10.

```julia
using SpeedyWeather, SpeedyWeatherERA5, Dates
ENV["ERA5_CACHE"] = "/some/large/disk/era5" 
sg = SpectralGrid(truncation = 43, nlayers = 8, Grid = FullGaussianGrid)   # FullGaussianGrid is required
model = PrimitiveWetModel(sg; ocean = PrescribedOcean(sg), sea_ice = PrescribedSeaIce(sg))
sim, ic = era5_initialize!(model, DateTime(2021, 6, 26))
```

Or, to run with land initialised from ERA5 also: 

```julia
using SpeedyWeather, SpeedyWeatherERA5, Dates
ENV["ERA5_CACHE"] = "/path/to/large/disk/era5"
dt = DateTime(2021, 6, 26) #for convenience
sg    = SpectralGrid(truncation = 43, nlayers = 8, Grid = FullGaussianGrid)
model = PrimitiveWetModel(sg; ocean = PrescribedOcean(sg), sea_ice = PrescribedSeaIce(sg))  # fresh model for every initialisation`

sim, ic = era5_initialize!(model, dt)
land = prepare_era5_land(model, get_era5(dt), get_era5_land(dt))

set_era5_land!(sim, land; temperature = true, moisture = :relative, snow = true)

run!(sim, period = Day(7))
```

## Current code (src/)
| file | contents |
|---|---|
| `SpeedyWeatherERA5.jl` | module, includes, exports |
| `cache.jl` | `get_era5`, `get_era5_land`, `era5_path`, `era5_land_path`, `default_cache`. Knows nothing about SpeedyWeather |
| `surface_pressure.jl` | `surface_temperature`, `raw_pressure`, `orography_consistent_pressure` |
| `prepare.jl` | `pressure_to_sigma`, `native_to_model`, `masked_to_model`, `prepare_era5` |
| `land.jl` | `soil_layer_weights`, `prepare_era5_land`, `set_era5_land!` |
| `initialize.jl` | `set_era5!`, `era5_initialize!` |

Flow of `era5_initialize!(model, dt; cache, kwargs...)`: `get_era5(dt)` -> `initialize!(model; time = dt)` -> `prepare_era5(model, path; kwargs...)` -> `set_era5!(sim, ic)`;
returns `(sim, ic)`.  This leaves land on climatology.

`set_era5!` itself returns `sim`. 

Land is NOT part of this chain: call `get_era5_land(dt)`, `prepare_era5_land(model, get_era5(dt), land_path)`
and `set_era5_land!(sim, land; temperature, moisture, snow)` afterwards.

For more detailed notes on how data is assimilated, see NOTES.md