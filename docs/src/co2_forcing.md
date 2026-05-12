# CO₂ forcing experiments 

With our [greenhouse gases](radiation.md#greenhouse-gases) and [radiation interface](@ref radiation-ssm), we can also run idealized climate forcing scenarios. We will perform an experiment in which we let the CO₂ concentration rise exponentially starting from 280 ppmv and compare the globally averaged surface temperature to a control run in which we keep the CO₂ concentration constant at 280 ppmv. 

To track the globally averaged surface temperature, we first set up a callback that computes and saves it. Then we set up the two simulations, both using the Simple Spectral Radiation Model of [`AnalyticBandRadiation.jl`](@ref radiation-ssm). This script takes significant compute, so that it is not built live in the documentation. 

```julia 
using Statistics
using Dates
using Makie

using SpeedyWeather
using AnalyticBandRadiation
const SpeedyAnalyticBand = Base.get_extension(AnalyticBandRadiation, :AnalyticBandRadiationSpeedyWeatherExt)

# Custom callback to compute and save global mean surface temperature
mutable struct GlobalMeanTemperatureCallback <: SpeedyWeather.AbstractCallback
    temperatures::Vector{Float32}
    timestep_counter::Int
    save_interval::Int  # save every N timesteps
end

GlobalMeanTemperatureCallback(; save_interval=1) =
    GlobalMeanTemperatureCallback(Float32[], 0, save_interval)

function SpeedyWeather.initialize!(
    callback::GlobalMeanTemperatureCallback,
    ::SpeedyWeather.Variables,
    ::SpeedyWeather.AbstractModel,
)
    # Pre-allocate storage for one sample per day
    # Assuming ~24 timesteps per day with typical parameters
    n_days = 365 * 10  # 10 years
    callback.temperatures = zeros(Float32, n_days)
    callback.timestep_counter = 0
    return nothing
end

function SpeedyWeather.callback!(
    callback::GlobalMeanTemperatureCallback,
    vars::SpeedyWeather.Variables,
    model::SpeedyWeather.AbstractModel,
)
    callback.timestep_counter += 1

    # Save every N timesteps
    if mod(callback.timestep_counter, callback.save_interval) == 0
        save_index = div(callback.timestep_counter, callback.save_interval)

        # Get surface layer temperature (last layer)
        nlayers = model.geometry.nlayers
        temp_grid = vars.grid.temperature[:, nlayers]

        # Compute global mean
        global_mean = mean(temp_grid)

        # Store result
        if save_index <= length(callback.temperatures)
            callback.temperatures[save_index] = global_mean
        else
            push!(callback.temperatures, global_mean)
        end
    end
    return nothing
end

# Setup
arch = SpeedyWeather.CPU()
spectral_grid = SpectralGrid(trunc=31, nlayers=8, architecture=arch)

# Control simulation (280 ppm CO2)
radiation = SpeedyAnalyticBand.SpeedyAnalyticBandLongwave(spectral_grid)
model_control = PrimitiveWetModel(spectral_grid, longwave_radiation=radiation)
# get the save_interval from the time stepping of the model
callback_control = GlobalMeanTemperatureCallback(save_interval=Int(round(Day(1) / Second(model_control.time_stepping.Δt_sec))))
add!(model_control, :surface_temp_callback => callback_control)

simulation_control = initialize!(model_control)

println("Running control simulation (280 ppm CO2) for 10 years...")
run!(simulation_control, period=Year(10))
println("Control simulation complete!")

# Simulation with exponential CO2 increase
radiation = SpeedyAnalyticBand.SpeedyAnalyticBandLongwave(spectral_grid)
model_exp = PrimitiveWetModel(spectral_grid, longwave_radiation=radiation)
radiation = SpeedyAnalyticBand.SpeedyAnalyticBandLongwave(spectral_grid)
simulation_exp = initialize!(model_exp)
callback_exp = GlobalMeanTemperatureCallback(save_interval=Int(round(Day(1) / Second(model_control.time_stepping.Δt_sec))))
add!(model_exp, :surface_temp_callback => callback_exp)

println("Running experimental simulation (elevated CO2) for 10 years...")
run!(simulation_exp, period=Year(10))
println("Experimental simulation complete!")

# Filter out unused slots
temps_control = callback_control.temperatures[callback_control.temperatures .!= 0]
temps_exp = callback_exp.temperatures[callback_exp.temperatures .!= 0]

# Create time axis (days)
n_samples = length(temps_control)
time_days = 1:n_samples

# Create Makie figure
fig = Figure(size=(1200, 600))
ax = Axis(
    fig[1, 1],
    xlabel="Time (days)",
    ylabel="Global Mean Surface Temperature (K)",
    title="Surface Temperature Evolution: Control vs Elevated CO2 (10-year simulations)",
)

# Plot both time series
lines!(ax, time_days, temps_control, label="Control (280 ppm CO2)", linewidth=2, color=:blue)
lines!(ax, time_days, temps_exp, label="Exp CO2 Rise", linewidth=2, color=:red)

# Add legend
axislegend(ax, position=:rb)

# Add grid
grid!(ax, alpha=0.3)

display(fig)

# Optionally save the figure
save("global_warming.png", fig)
```