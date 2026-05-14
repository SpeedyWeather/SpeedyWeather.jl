module SpeedyWeatherBrowzarrExt

using SpeedyWeather
using Browzarr

function Browzarr.browzarr(simulation::SpeedyWeather.AbstractSimulation; kwargs...)
    # for netCDF it's a file, for zarr it's a directory call both "store"
    store_path = SpeedyWeather.get_output_path(simulation)
    return Browzarr.browzarr(; store = store_path, kwargs...)
end

end
