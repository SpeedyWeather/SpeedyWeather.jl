# SpeedyWeather.jl
[![CI](https://github.com/milankl/SpeedyWeather.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/milankl/SpeedyWeather.jl/actions/workflows/CI.yml)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://milankl.github.io/SpeedyWeather.jl/dev)

SpeedyWeather.jl is a global atmospheric model with simple physics which is developed as a computational playground
with an everything-flexible attitude as long as it is speedy. With minimal code redundancies it supports

- any number format and precision (32 and 64 bits, 16 bits in development)
- different architectures (x86 and arm, GPUs in development)
- different physical models (barotropic vorticity and shallow water, primitive equations in development)
- different spatial grids (full and octahedral grids, Gaussian and Clenshaw-Curtis, HEALPix)
- different horizontal resolutions (T31 to T1023, i.e. 400km to 10km using linear, quadratic or cubic truncation)

and Julia will compile to these choices just-in-time. Parallelisation is designed (but not fully implemented yet)
towards efficiency within a compute node and at rather low resolutions: The dynamical core layer-by-layer in the vertical
and the physical parameterizations gridpoint-by-gridpoint in the horizontal.

For an overview of the functionality and explanation see the (always somehow incomplete)
[documentation](https://milankl.github.io/SpeedyWeather.jl/dev).

## Example use

SpeedyWeather.jl is currently developed. Some things work, some don't. Stay tuned.
Here is video of some relative vorticity, simulated at T682 (20km at the Equator) spectral resolution.

https://user-images.githubusercontent.com/25530332/185218396-620a0887-3860-496f-a265-aa59a2079768.mp4

The main interface to SpeedyWeather.jl is 

```julia
julia> using SpeedyWeather
julia> run_speedy(n_days=30, trunc=63, Grid=OctahedralGaussianGrid, model=:shallowwater, output=true)
Weather is speedy run 1: 100%|███████████████████████████████████| Time: 0:00:04 (1498.70 years/day)
```

Hurray! In 4 seconds we just simulated 30 days of the Earth's atmosphere at a speed of 1500 years per day.
This simulation used a T63 spectral resolution on an octahedral Gaussian grid (165km resolution) solving
the shallow water equations. The arguments for `run_speedy` are described in
[`src/default_parameters.jl`](https://github.com/milankl/SpeedyWeather.jl/blob/main/src/default_parameters.jl).

## History

SpeedyWeather.jl is a reinvention of the atmospheric general circulation model
[SPEEDY](http://users.ictp.it/~kucharsk/speedy-net.html) in Julia. While conceptually the same model,
it is entirely restructured, features have been added, changed and removed, such that only the core
algorithms share similarities. Speedy's dynamical core is originally written by Isaac Held
at GFDL and the physical parametrizations by Fred Kucharski, Franco Molteni and Martin P. King in Fortran77.
Speedy was then translated to Fortran90 by Sam Hatfield in [speedy.f90](https://github.com/samhatfield/speedy.f90).
SpeedyWeather.jl is then adopted from [first translations to Julia](https://github.com/samhatfield/speedy.jl) by Sam Hatfield.

## Installation

SpeedyWeather.jl is registered in Julia's registry, so open the package manager with `]` and
```julia
(@v1.7) pkg> add SpeedyWeather
```
which will install the latest release and all dependencies automatically. The very latest version is installed with
```julia
(@v1.7) pkg> add https://github.com/milankl/SpeedyWeather.jl#main
```
which pulls directly from the `#main` branch. Please use the current minor version of Julia,
compatibilities with older versions are not guaranteed.
