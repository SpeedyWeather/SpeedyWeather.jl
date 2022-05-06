# SpeedyWeather.jl
[![CI](https://github.com/milankl/SpeedyWeather.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/milankl/SpeedyWeather.jl/actions/workflows/CI.yml)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://milankl.github.io/SpeedyWeather.jl/dev)

The focus of SpeedyWeather.jl is to develop a global atmospheric model with simple physics, that can run at various levels of
precision (16, 32 and 64-bit) on different architectures (x86 and ARM are currently planned, GPUs maybe in the future).

Additionally, the model is written in an entirely number format-flexible way, such that any custom number format can
be used and Julia will compile to the format automatically. In contrast to the original speedy the resolution
(horizontal and vertical) is adjustable and some simple parallelism is included to run SpeedyWeather.jl efficiently
on small clusters.

For an overview of the functionality and explanation see the
[documentation](https://milankl.github.io/SpeedyWeather.jl/dev).

## Example use

SpeedyWeather.jl is currently developed. Some things work, some don't.
Stay tuned. Here is a teaser picture of some barotropic vorticity, simulated at T341 spectral resolution in single precision.

![vor](docs/img/barotropic_vorticity.png?raw=true "First simulation of barotropic vorticity with SpeedyWeather.jl")

The main interface to SpeedyWeather.jl is 

```julia
julia> using SpeedyWeather
julia> run_speedy(Float32,n_days=100,output=true);
Weather is speedy run 1:  49%|███████▍       |  ETA: 0:00:09 ( 3.71 ms/it)
```

and the arguments for `run_speedy` are described in [`src/default_parameters.jl`](https://github.com/milankl/SpeedyWeather.jl/blob/main/src/default_parameters.jl).

## History

SpeedyWeather.jl is the Julia version of the atmospheric general circulation model
[SPEEDY](http://users.ictp.it/~kucharsk/speedy-net.html). Speedy's dynamical core is originally written by Isaac Held
at GFDL and the physical parametrizations by Fred Kucharski, Franco Molteni and Martin P. King in Fortran77.
Speedy was then translated to Fortran90 by Sam Hatfield in [speedy.f90](https://github.com/samhatfield/speedy.f90).
SpeedyWeather.jl is then adopted from [first translations to Julia](https://github.com/samhatfield/speedy.jl) by Sam Hatfield,
but entirely restructured and only the parts of the conceptual algorithms are shared with the original Fortran versions.

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
