# Changelog

## Unreleased

- RingGrids: To wrap an Array with the horizontal dimension in matrix shape into a full grid, one has to use e.g. `FullGaussianGrid(map, input_as=Matrix)` now. [#572](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/572)
* RingGrids: Fixed a bug in `RingGrids`, so that now broadcasts are defined even when the dimensions mismatch in some cases 
  
## v0.11.0
