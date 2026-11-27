# SpeedyWeatherERA5: design notes and validation log

Developed against SpeedyWeather.jl main @ 14aa52b (2026-10-16), version 0.23.0+DEV, monorepo layout, Julia 1.10.12. This doc explains the different stages in more detail than the readme, and includes scientific decisions/equations used.


## Stages:

### Data Acquisition
- `get_era5`: cache, knows nothing about SpeedyWeather (`src/cache.jl`):

    - Fetches one hour of ARCO-ERA5 from Google Cloud (public) with Zarr.jl, at a 0.25deg grid, with 37 pressure levels. 

    - Writes one NetCDF into a cache dir, ENV["ERA5_CACHE"] (as above) sets this directory location

    - Variables: 3D `temperature, specific_humidity, u/v_component_of_wind`; 2D `surface_pressure, geopotential_at_surface, land_sea_mask, sea_surface_temperature, sea_ice_cover, snow_depth, soil_temperature_level_1/2, volumetric_soil_water_layer_1/2`. The land file has all four soil layers of temperature and water. 
- `get_era5_land`: same thing but for land, maps ERA5's four soil layers onto the model's two layers by depth overlap. 

### Preparing Initial Conditions

- `prepare_era5`: pressure to sigma, orography-consistent surface pressure spectral truncation, ocean-only averaging of SST and sea ice (`src/prepare.jl, src/surface_pressure.jl`)
    - **$p_s$ correction**: Something to note here is that ERA5's surface pressure is for ERA5's orography, but SpeedyWeather's orography is different. When converting ERA5's $p_s$ to SpeedyWeather's, you can either leave $p_s$ as it is (`:raw`) or shift it hydrostatically (default, `:consistent`): 
        - $ln(p_{s, speedy}) = ln(p_{s, ERA5}) + \frac{g \cdot (h_{ERA5} - h_{speedy})}{R \cdot T_s}$
        - Shifting $p_s$ hydrostatically resulted in far lower RMSE for windspeeds, so the latter was incorporated as default.

    - **Vertical interpolation**: the `pressure_to_sigma` is an analogue of the existing `interpolate_pressure_layers` function, but in the other direction. It samples each vertical column as $p = \sigma \cdot p_{s, speedy}$. 
        - if p lies between two above-ground ERA5 levels, then the value of each variable ($T, q, u, v$) is linearly interpolated in $ln(p)$. 
    SpeedyWeather's lowest sigma layers can sit below ERA5's ground, so here variables need to be extrapolated (where $p > \sigma \cdot p_{s, ERA5}$). 
        - T is extrapolated with a 6.5K/km lapse rate, where $T(p) = T_0 \cdot \left( \frac{p}{p_0} \right) ^{\frac{R \Gamma}{g}}$, combining the hydrostatic equ and ideal gas equ. 
        - $q, u, v$ are held constant at the lowest above-ground ERA5 level. 
    - **Horizontal truncation** once vertical steps are done at native resolution, `native_to_model` takes ERA5 Reduced Gaussian data, drops the pole rows (so it becomes a FullClenshawGrid), and then the data is truncated spectrally (to T*n+1* and transformed to a FullGaussianGrid. 
        - Only FullGaussianGrid is supported because SST/ice/and averaging assumes it in code (in `masked_to_model`)
    SST and sea ice are not spectrally truncated, `masked_to_model` does an ocean-only cos(lat)-weghted box average of native cells into each model cell.
If using an ERA5 land initialisation, `prepare_era5_land` comes into play:
    - **Land variables** Here, the ERA5 layers do not match SpeedyWeather's two-layer system.
        - Soil temperature and soil moisture at each SpeedyWeather layer are the overlap-weighted mean of the ERA5 layers.
        - Snow uses ERA5's snow_depth directly, capped at SpeedyWeather's `snow_depth_cap` (10m)
        - Horizontally, variables are averaged onto the model grid using ERA5 land cells, and cells without land are filled from neighbours and set to SpeedyWeather's default ocean values. 

### Setting initial conditions
 - `set_era5!`: sets the remapped variables as the initialisation inputs for SpeedyWeather (`src/initialize.jl`), or if including land variables, `set_era5_land!` (`src/land.jl`)
    - Uses SpeedyWeather's existing `set!` to return (sim, ic) (SpeedyWeather simulation, prepared initial conditions)
    - **Setting land from ERA5** set_era5_land! is separate and optional
    - **Unphysical values:** SST is capped at the freezing temperature and land points are masked
    - Works for `SlabOcean/ThermodynamicSeaIce` and for `PrescribedOcean/PrescribedSeaIce`
 - `set_era5_land!`: an additional step for the land initialisation 

## Additional Notes:
- Re-initialising a model changes its land-sea mask: always use a fresh PrimitiveWetModel per initialisation.
- Coarse resolutions (T10) seem to cause problems 
- Land soil moisure flagged as inconsistent: SpeedyWeather either uses field capacity (default, `:relative`) or volumetric convention (`:volumetric`, used in climatology code)
- The above was mostly verified by running `truncation=64, Grid=FullGaussianGrid, nlayers=8, Float32` (T63)

Stopped at Step 1 result
