# Stage 1: acquire + cache. Knows nothing about SpeedyWeather.
# ERA5 (ARCO, native 0.25°, 37 pressure levels, hourly) -> one NetCDF per datetime in the cache dir.

const ARCO_URL = "https://storage.googleapis.com/gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
"""Default cache directory: `ENV["ERA5_CACHE"]` if set, else `~/.cache/SpeedyWeatherERA5`. Never the working directory."""
default_cache() = get(ENV, "ERA5_CACHE", joinpath(homedir(), ".cache", "SpeedyWeatherERA5"))

const VARS_3D = ["temperature", "specific_humidity", "u_component_of_wind", "v_component_of_wind"]
const VARS_2D = ["surface_pressure", "geopotential_at_surface", "land_sea_mask", "sea_surface_temperature", "sea_ice_cover",
                 "snow_depth", "soil_temperature_level_1", "soil_temperature_level_2",
                 "volumetric_soil_water_layer_1", "volumetric_soil_water_layer_2"]

cached_variables(path) = NCDataset(ds -> split(ds.attrib["variables"], ","), path)

era5_path(dt::DateTime; cache = default_cache()) = joinpath(cache, "era5_native_" * Dates.format(dt, "yyyymmddTHH") * ".nc")

"""Fetch `variables3d` (37 pressure levels) and `variables2d` at `dt` from ARCO-ERA5 into the NetCDF `path` (atomic)."""
function _fetch_slice!(path, dt, variables3d, variables2d)
    @assert dt == floor(dt, Hour) "ERA5 is hourly, got $dt"
    mkpath(dirname(path))

    g = zopen(ARCO_URL; consolidated = true)
    i = Dates.value(Hour(dt - DateTime(1900, 1, 1))) + 1        # `time` counts hours since 1900-01-01 (1-based index)
    @assert g["time"][i] == i - 1 "unexpected time axis in ARCO store"
    @assert dt >= DateTime(Date(g.attrs["valid_time_start"][1:10])) "$dt is before the start of the ARCO ERA5 record"

    tmp = path * ".part"
    ds = NCDataset(tmp, "c")
    try
        ds.attrib["source"] = ARCO_URL
        ds.attrib["datetime"] = Dates.format(dt, "yyyy-mm-ddTHH:MM:SS")
        ds.attrib["retrieved"] = string(now(UTC)) * "Z"
        ds.attrib["store_last_updated"] = string(g.attrs["last_updated"])
        ds.attrib["variables"] = join([variables3d; variables2d], ",")
        defVar(ds, "longitude", Float32.(g["longitude"][:]), ("longitude",); attrib = Dict("units" => "degrees_east"))
        defVar(ds, "latitude", Float32.(g["latitude"][:]), ("latitude",); attrib = Dict("units" => "degrees_north"))
        defVar(ds, "level", Int32.(g["level"][:]), ("level",); attrib = Dict("units" => "hPa"))
        for v in variables3d
            x = defVar(ds, v, Float32, ("longitude", "latitude", "level"); deflatelevel = 3, attrib = Dict("units" => g[v].attrs["units"]))
            x[:, :, :] = g[v][:, :, :, i]
        end
        for v in variables2d
            x = defVar(ds, v, Float32, ("longitude", "latitude"); deflatelevel = 3, attrib = Dict("units" => g[v].attrs["units"]))
            x[:, :] = g[v][:, :, i]
        end
    finally
        GC.gc()                                                 # drop variable handles before closing (netCDF-C segfaulted in nc_close otherwise, rarely)
        close(ds)
    end
    mv(tmp, path; force = true)                                 # atomic: a half-written file never looks cached
    return path
end

"""Path to the cached native-grid ERA5 slice at `dt` (whole hours), fetched once from ARCO-ERA5 if not cached.
Knows nothing about SpeedyWeather. `variables3d` are on the 37 pressure levels, `variables2d` are single-level."""
function get_era5(dt::DateTime; cache = default_cache(), variables3d = VARS_3D, variables2d = VARS_2D)
    path = era5_path(dt; cache)
    isfile(path) && cached_variables(path) ⊇ [variables3d; variables2d] && return path   # re-fetch if a needed variable is missing
    return _fetch_slice!(path, dt, variables3d, variables2d)
end

const VARS_LAND = ["soil_temperature_level_1", "soil_temperature_level_2", "soil_temperature_level_3", "soil_temperature_level_4",
                   "volumetric_soil_water_layer_1", "volumetric_soil_water_layer_2", "volumetric_soil_water_layer_3", "volumetric_soil_water_layer_4"]

era5_land_path(dt::DateTime; cache = default_cache()) = joinpath(cache, "era5_land_native_" * Dates.format(dt, "yyyymmddTHH") * ".nc")

"""Path to the cached ERA5 soil temperature and soil moisture (all four ERA5 soil layers, 0-7, 7-28, 28-100, 100-289 cm) at `dt`,
a small separate file (a few tens of MB) so that the large atmospheric file is not refetched."""
function get_era5_land(dt::DateTime; cache = default_cache())
    path = era5_land_path(dt; cache)
    isfile(path) && cached_variables(path) ⊇ VARS_LAND && return path
    return _fetch_slice!(path, dt, String[], VARS_LAND)
end
