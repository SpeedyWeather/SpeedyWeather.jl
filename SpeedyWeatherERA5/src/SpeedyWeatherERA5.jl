"""
SpeedyWeatherERA5: initialise a SpeedyWeather.jl PrimitiveWetModel from ERA5.

    get_era5(dt)                 acquire + cache one hour of ARCO-ERA5 on its native 0.25° grid (no SpeedyWeather inside)
    prepare_era5(model, path)    pressure -> sigma, orography-consistent surface pressure, truncation, SST/ice averaging
    set_era5!(sim, ic)           set! the state of a simulation
    era5_initialize!(model, dt)  all of the above
"""
module SpeedyWeatherERA5

using Dates
using Zarr, NCDatasets
using SpeedyWeather, SpeedyTransforms, RingGrids
import SpeedyWeatherInternals

include("cache.jl")
include("surface_pressure.jl")
include("prepare.jl")
include("land.jl")
include("initialize.jl")

export get_era5, get_era5_land, era5_path, era5_land_path, default_cache, prepare_era5, native_to_model, masked_to_model, pressure_to_sigma
export raw_pressure, orography_consistent_pressure
export set_era5!, era5_initialize!, prepare_era5_land, set_era5_land!

end
