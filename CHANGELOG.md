# Changelog

## Unreleased

- Random processes for random pattern generation [#592](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/592)
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