module SpeedyWeatherInternalsMLDataDevicesExt

import MLDataDevices
import SpeedyWeatherInternals.Architectures: CPU, GPU, ReactantDevice

# Map a SpeedyWeather architecture to the corresponding MLDataDevices device.
# `gpu_device()`/`reactant_device()` select the actual device based on which
# backend package is loaded and functional, mirroring the architecture at hand.
MLDataDevices.get_device(::CPU) = MLDataDevices.cpu_device()
MLDataDevices.get_device(::GPU) = MLDataDevices.gpu_device()
MLDataDevices.get_device(::ReactantDevice) = MLDataDevices.reactant_device()

end
