# SpeedyTransforms.jl

SpeedyTransforms is a submodule of [SpeedyWeather.jl](https://github.com/SpeedyWeather/SpeedyWeather.jl) that provides fast spherical harmonic transforms for atmospheric modeling. This module implements efficient algorithms for transforming data between spectral (spherical harmonic) and gridpoint representations on the sphere.

SpeedyTransforms uses LowerTriangularArrays.jl to represent
n-dimensional (2D sphere + arbitrary others) data as coefficients
of the spherical harmonics (the spectral space) and RingGrids.jl
to represent gridded n-dimensional data on the sphere (the grid space).

## Usage

SpeedyTransforms is automatically available when using SpeedyWeather.jl:

```julia
using SpeedyWeather
# SpeedyTransforms functionality is integrated into the model
# or use `SpeedyTransforms.` to reach its scope
```

For direct use of transform functionality:

```julia
using SpeedyWeather: SpeedyTransforms
```

## Documentation

See 