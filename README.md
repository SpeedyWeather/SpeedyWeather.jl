# SpeedyWeather.jl
[![CI](https://github.com/SpeedyWeather/SpeedyWeather.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/SpeedyWeather/SpeedyWeather.jl/actions/workflows/CI.yml)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://speedyweather.github.io/SpeedyWeather.jl/stable/)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://speedyweather.github.io/SpeedyWeather.jl/dev/)

SpeedyWeather.jl is a global spectral atmospheric model with simple physics which is developed as a research playground
with an everything-flexible attitude as long as it is speedy. With minimal code redundancies it supports

- different physical models (barotropic vorticity, shallow water, primitive equations dry and wet core)
- physics parameterizations for precipitation, boundary layer, etc; new ones are passed as argument to the model constructor
- different spatial grids (full and octahedral grids, Gaussian and Clenshaw-Curtis, HEALPix, OctaHEALPix)
- different horizontal resolutions (T31 to T1023 and higher, i.e. 400km to 10km using linear, quadratic or cubic truncation)
- 32-bit single and 64-bit double precision (16 bits in development)
- multi-threading, layer-wise for dynamics, grid point-wise for physics (multi-processing in development)
- different architectures (x86 and arm, GPUs in development)

and Julia will compile to these choices just-in-time.

For an overview of the functionality and explanation see the
[documentation](https://speedyweather.github.io/SpeedyWeather.jl/dev).

## Example use

SpeedyWeather.jl is currently developed. Some things work, some don't. Stay tuned or talk to us by raising an
[issue]([url](https://github.com/SpeedyWeather/SpeedyWeather.jl/issues))
to express interest as a user or developer. All contributions are always welcome.

With v0.6 the interface to SpeedyWeather.jl consist of 4 steps: define the grid, create the model, initialize, run

```julia
spectral_grid = SpectralGrid(trunc=31, Grid=OctahedralGaussianGrid, nlev=8)
model = PrimitiveDryModel(;spectral_grid, orography = EarthOrography(spectral_grid))
simulation = initialize!(model)
run!(simulation,period=Day(10),output=true)
```
and you will see

<img src="https://github.com/SpeedyWeather/SpeedyWeather.jl/assets/25530332/a04fbb10-1cc1-4f77-93f2-7bdf047f277d" width="450"><br>

Hurray🥳 In 5 seconds we just simulated 10 days of the Earth's atmosphere at a speed of 440 years per day.
This simulation used a T31 spectral resolution on an
[octahedral Gaussian grid](https://speedyweather.github.io/SpeedyWeather.jl/dev/grids/#Implemented-grids)
(~400km resolution) solving the dry primitive equations on 8 vertical levels.
The [UnicodePlot](https://github.com/JuliaPlots/UnicodePlots.jl) will give
you a snapshot of surface vorticity at the last time step. The plotted resolution is not representative,
but allows a quick check of what has been simulated. The [NetCDF output](https://speedyweather.github.io/SpeedyWeather.jl/dev/output/) is independent of the unicode plot.

More examples in the [How to run SpeedyWeather](https://speedyweather.github.io/SpeedyWeather.jl/dev/how_to_run_speedy/)
section of the [documentation](https://speedyweather.github.io/SpeedyWeather.jl/dev).

## Gallery

Here is video of some relative vorticity in the shallow water model, simulated at T1023 spectral resolution
(about 10km) on an
[octahedral Clenshaw-Curtis grid](https://speedyweather.github.io/SpeedyWeather.jl/dev/grids/#Implemented-grids)
with more than 4 million grid points

https://user-images.githubusercontent.com/25530332/190443050-d5b8d093-86c0-46c9-b515-8420059ac8dc.mp4

The primitive equation core (wet or dry) is in development, this is temperature at the surface
at T511 (~20km resolution) and 31 vertical levels. The simulation was multi-threaded in Float32 (single precision).

https://github.com/SpeedyWeather/SpeedyWeather.jl/assets/25530332/95897b82-9b81-4980-934b-cfdcf4d5a4b0

SpeedyWeather.jl can also solve the 2D barotropic vorticity equations on the sphere.
Here, we use single-threaded Float32 (single precision) at a resolution of T340 (40km) on
an [octahedral Gaussian grid](https://speedyweather.github.io/SpeedyWeather.jl/dev/grids/#Implemented-grids). 
Initial conditions are randomly distributed relative vorticity on a slowly rotating Earth ($\Omega = 10^{-6}\text{ s}^{-1}$) and no forcing is applied

https://github.com/SpeedyWeather/SpeedyWeather.jl/assets/25530332/8a7c6758-950f-424d-8ece-0480295386b3

## History

SpeedyWeather.jl is a reinvention of the atmospheric general circulation model
[SPEEDY](http://users.ictp.it/~kucharsk/speedy-net.html) in Julia. While conceptually a similar model,
it is entirely restructured, features have been added, changed and removed, such that only the numerical
schemes share similarities (but we start to diverge here too). Speedy's dynamical core has an obscure history:
Large parts were written by Isaac Held at GFDL in/before the 90ies with an unknown amount of
contributions/modifications from Steve Jewson (Oxford) in the 90ies.
The physical parametrizations were then added by Franco Molteni, Fred Kucharski, and Martin P. King
afterwards while the model was still written in Fortran77.
Around 2018-19, SPEEDY was then translated to Fortran90 by Sam Hatfield in [speedy.f90](https://github.com/samhatfield/speedy.f90).
SpeedyWeather.jl is then adopted from [first translations to Julia](https://github.com/samhatfield/speedy.jl) by Sam Hatfield.

## Submodules

SpeedyWeather.jl defines several submodules that are technically stand-alone (with dependencies) but aren't separated
out to their own packages for now

- [__RingGrids__](https://speedyweather.github.io/SpeedyWeather.jl/dev/ringgrids/),
a module that defines several iso-latitude ring-based spherical grids (like the FullGaussianGrid or the HEALPixGrid)
and interpolations between them
- [__LowerTriangularMatrices__](https://speedyweather.github.io/SpeedyWeather.jl/dev/lowertriangularmatrices/),
a module that defines `LowerTriangularMatrix` used for the spherical harmonic coefficients
- [__SpeedyTransforms__](https://speedyweather.github.io/SpeedyWeather.jl/dev/speedytransforms/), a module that defines
the spherical harmonic transform between spectral space (for which LowerTriangularMatrices is used) and grid-point space
(as defined by RingGrids).

These modules can also be used independently of SpeedyWeather like so
```julia
julia> using SpeedyWeather: LowerTriangularMatrices, RingGrids, SpeedyTransforms
```
check out their documentation: [RingGrids](https://speedyweather.github.io/SpeedyWeather.jl/dev/ringgrids/), 
[LowerTriangularMatrices](https://speedyweather.github.io/SpeedyWeather.jl/dev/lowertriangularmatrices/), 
[SpeedyTransforms](https://speedyweather.github.io/SpeedyWeather.jl/dev/speedytransforms/).

## Installation

SpeedyWeather.jl is registered in Julia's registry, so open the package manager with `]` and
```julia
(@v1.8) pkg> add SpeedyWeather
```
which will install the [latest release]([url](https://github.com/SpeedyWeather/SpeedyWeather.jl/releases))
and all dependencies automatically. For more information see the
[Installation](https://speedyweather.github.io/SpeedyWeather.jl/dev/installation/) in the documentation.
Please use the current minor version of Julia,
compatibilities with older versions are not guaranteed.

## Copyright and license

Copyright (c) 2020 Milan Klöwer for SpeedyWeather.jl  
Copyright (c) 2021 The SpeedyWeather.jl Contributors for SpeedyWeather.jl  
Copyright (c) 2022 Fred Kucharski and Franco Molteni for SPEEDY parametrization schemes  

Software licensed under the [MIT License](LICENSE.txt).
