# SpeedyWeatherInternals.jl

SpeedyWeatherInternals is a submodule of [SpeedyWeather.jl](https://github.com/SpeedyWeather/SpeedyWeather.jl) that contains internal utility functions, constants, and helper types used throughout the SpeedyWeather ecosystem. This module is primarily for internal use and provides the foundational components that support the atmospheric modeling capabilities of SpeedyWeather.

## Overview

SpeedyWeatherInternals contains two modules

- **Architectures** - Defines computing architectures (CPU/GPU) and how to move data between them
- **Utils** - Defines kernel launching infrastructure to use KernelAbstractions in SpeedyWeather, and other miscellaneous functionality

## Note for Users

SpeedyWeatherInternals is primarily intended for internal use within SpeedyWeather.jl. The API of this module may change without notice between versions as it supports the internal architecture. For stable, user-facing functionality, please refer to the main [SpeedyWeather.jl documentation](https://speedyweather.github.io/SpeedyWeatherDocumentation/dev/).

## Related Modules

- **RingGrids** - Spherical grid definitions and operations
- **LowerTriangularArrays** - Spectral coefficient storage
- **SpeedyTransforms** - Spherical harmonic transforms
- **SpeedyWeather** - The full atmospheric modeling library

All of these modules depend on SpeedyWeatherInternals.jl.

For more information about SpeedyWeather.jl, see the [SpeedyWeather documentation](https://speedyweather.github.io/SpeedyWeatherDocumentation/dev/).