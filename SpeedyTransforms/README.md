# SpeedyTransforms.jl

[![docs](https://img.shields.io/badge/documentation-latest_release-blue.svg)](https://speedyweather.github.io/SpeedyWeatherDocumentation/stable/speedytransforms/)
[![docs](https://img.shields.io/badge/documentation-main-blue.svg)](https://speedyweather.github.io/SpeedyWeatherDocumentation/dev/speedytransforms/)

SpeedyTransforms is a submodule of [SpeedyWeather.jl](https://github.com/SpeedyWeather/SpeedyWeather.jl) that provides fast spherical harmonic transforms for atmospheric modeling. This module implements efficient algorithms for transforming data between spectral (spherical harmonic) and gridpoint representations on the sphere. It also implements gradient operators ``\nabla, \nabla \cdot, \nabla \times, \nabla^2, \nabla^{-2}`` in spectral space.

SpeedyTransforms uses LowerTriangularArrays.jl to represent
n-dimensional (2D sphere + arbitrary others) data as coefficients
of the spherical harmonics (the spectral space) and RingGrids.jl
to represent gridded n-dimensional data on the sphere (the grid space).

SpeedyTransforms is written to be device agnostic and differentiable by Enzyme. 

## Example Usage

SpeedyTransforms is automatically available when using SpeedyWeather.jl, but can also be installed stand-alone:

```julia
using SpeedyWeather
# SpeedyTransforms functionality is integrated into the model
# or use `SpeedyTransforms.` to reach its scope
```

For direct use of transform functionality in which case we also need to load RingGrids and LowerTriangualarArrays for the spectral and grid data formats:

```julia
using SpeedyTransforms, RingGrids, LowerTriangularArrays
```

With the package loaded we can now use the transform functionality:

```julia
using SpeedyTransforms
```

With the package loaded we can now `transform` data between spectral and gridpoint space:

```julia
alms = rand(LowerTriangularArray{Float32}, 5, 5) # spectral coefficients
map = transform(alms) # grid space 
```

and back 

```julia
alms2 = transform(map)
```

For more information, e.g. how to pre-plan the transforms and how to use the gradient implementations, see the [documentation](https://speedyweather.github.io/SpeedyWeatherDocumentation/dev/speedytransforms/) for more information on its usage.

## Related Modules

SpeedyTransforms.jl is designed to be used in conjuction with the other modules of the SpeedyWeather.jl ecosystem: 
 
- **LowerTriangularArrays** - Spectral coefficient storage
- **RingGrids** - Spherical grid definitions and operations
- **SpeedyWeather** - The full atmospheric modeling library
