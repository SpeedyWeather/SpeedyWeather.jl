module SpeedyWeatherBrowzarrExt

using SpeedyWeather
using Browzarr

function Browzarr.browzarr(simulation::SpeedyWeather.AbstractSimulation; kwargs...)
    (; output) = simulation.model
    @assert output.active "Output has not been stored (yet) for this simulation. `run!(simulation)` with `output=true`."
    store_path = joinpath(output.run_path, output.filename)
    return Browzarr.browzarr(; store = store_path, kwargs...)
end

end
