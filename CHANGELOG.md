# Changelog

## Unreleased

- Optical depth introduced and array-agnostic ColumnVariables [#606](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/606)
- Include large-scale condensation tests [#615](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/615)
- bugfix: large-scale condensation also at <100% [#609](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/609)

## v0.12.1

- ConstantLandTemperature implemented [#612](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/612)
- set! for more boundary conditions [#611](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/611)
- SpectralFilter for horizontal diffusion [#601](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/601)
- GeoMakie weak dependency, globe function for 3D data visualisation [#600](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/600)
- Zonal mean for AbstractGridArray [#603](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/603)
- Rossby-Haurwitz wave with initial conditions for interface displacement for shallow water models[#604](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/604)

## v0.12.0

- OctaminimalGaussianArray/Grid to start with 4 points around the poles [#595](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/595)
- Output both accumulated and precipitation rate as netCDF [#596](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/596)
- Random processes for random pattern generation [#592](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/592)
- Also allow SpectralGrid as positional argument to model constructors [#593](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/593)
- De-interweave SpectralTransform [#587](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/587)
- Rossby-Haurwitz wave initial conditions [#591](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/591)
- Haversine formula and AbstractSphericalDistance [#588](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/588)
- Array-agnostic SpectralTransform [#583](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/583)
- Move CUDA dependency into extension [#586](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/586)
- Stop supporting Julia v1.9 [#585](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/585)
- `feedback.verbose` (de/activate the progressbar) is now set to `isinteractive()` to disable automatically for documentation [#582](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/582)

## v0.11.0

- Extend `set!` with `orography` keyword argument to set the orography with `set!` [#578](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/578)
- Added a new benchmark suite for the dynamics functions [#577](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/577)
- Introduced a new `set!` function that allows to set `PrognosticVariables` to new values with keyword arguments [#563](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/563)
- Restructured dynamical core with prognostic/diagnostic variables array-agnostic and 3-dimensional [#525](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/525)
- Modularised NetCDF output [#573](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/573)
- Fixed a bug in RingGrids, now broadcasts are defined even when the dimensions mismatch in some cases [#568](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/568)
- RingGrids: To wrap an Array with the horizontal dimension in matrix shape into a full grid, one has to use e.g. `FullGaussianGrid(map, input_as=Matrix)` now. [#572](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/572)
- CompatHelper: Allow for JLD2.jl v0.5
