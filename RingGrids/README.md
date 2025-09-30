# RingGrids.jl

[![docs](https://img.shields.io/badge/documentation-latest_release-blue.svg)](https://speedyweather.github.io/SpeedyWeatherDocumentation/stable/ringgrids/)
[![docs](https://img.shields.io/badge/documentation-main-blue.svg)](https://speedyweather.github.io/SpeedyWeatherDocumentation/dev/ringgrids/)

RingGrids.jl is a package that has been developed for SpeedyWeather.jl but can also be used standalone. It's seperately registered as well, so you can install it via easily `using Pkg; Pkg.add("RingGrids")`.

RingGrids defines several iso-latitude grids that are discretizations of the space and the data format `Field` that stores data on these grids. 

## Example Use

In order to work with gridded data with RingGrids, we need to first define a grid and then create a field on that grid. All grids take in a resolution parameter `nlat_half` (the number of latitudes
on one hemisphere, Equator included) as an input argument. Then, we can create fields on that grid using the `zeros`, `ones`, `rand` functions as follows: 

```julia
using RingGrids

grid = HEALPixGrid(6)  # we define a HEALPix grid with nlat_half=6
field = zeros(grid)                 # 2D field
field = rand(grid, 10)              # 3D field with 10 vertical layers or time steps
field = randn(Float32, grid, 2, 2)  # 4D using Float32 as element type
```

RingGrids.jl also provised extensive functions to work with these fields and grids to e.g. visualize them or interpolate between them. The implementation is device agnostic so that it can be used on the CPU or GPU.

For more details see the [documentation](https://speedyweather.github.io/SpeedyWeatherDocumentation/dev/ringgrids/).

## Related Modules

RingGrids.jl is designed to be used in conjuction with the other modules of the SpeedyWeather.jl ecosystem: 

- **LowerTriangularArrays** - Spectral coefficient storage
- **SpeedyTransforms** - Spherical harmonic transforms
- **SpeedyWeather** - The full atmospheric modeling library
