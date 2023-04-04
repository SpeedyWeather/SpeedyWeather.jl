# SpeedyWeather.jl
[![CI](https://github.com/SpeedyWeather/SpeedyWeather.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/SpeedyWeather/SpeedyWeather.jl/actions/workflows/CI.yml)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://speedyweather.github.io/SpeedyWeather.jl/dev/)

SpeedyWeather.jl is a global atmospheric model with simple physics which is developed as a computational playground
with an everything-flexible attitude as long as it is speedy. With minimal code redundancies it supports

- any number format and precision (32 and 64 bits, 16 bits in development)
- different architectures (x86 and arm, GPUs in development)
- different physical models (barotropic vorticity, shallow water, primitive equations dry core; wet core in development)
- different spatial grids (full and octahedral grids, Gaussian and Clenshaw-Curtis, HEALPix, OctaHEALPix)
- different horizontal resolutions (T31 to T1023 and higher, i.e. 400km to 10km using linear, quadratic or cubic truncation)
- multi-threading, layer-wise for dynamics, grid point-wise for physics (multi-processing in development)

and Julia will compile to these choices just-in-time.

For an overview of the functionality and explanation see the (always somehow incomplete)
[documentation](https://speedyweather.github.io/SpeedyWeather.jl/dev).

## Gallery

Here is video of some relative vorticity in the shallow water model, simulated at T1023 spectral resolution (about 10km) on an
[octahedral Clenshaw-Curtis grid](https://github.com/milankl/SpeedyWeather.jl/issues/112#issuecomment-1219644323)
with more than 4 million grid points

https://user-images.githubusercontent.com/25530332/190443050-d5b8d093-86c0-46c9-b515-8420059ac8dc.mp4

The primitive equation core (wet or dry) is in development, this is temperature at the surface and at the tropopause
at T511 (~20km resolution) and 31 vertical levels. The simulation was multi-threaded in Float32 (single precision).
The orography is visible at the tropopause level because we currently use sigma coordinates

https://user-images.githubusercontent.com/25530332/229872856-bdcab69a-2226-4e9b-9470-9c9f90aa31e7.mp4

## Example use

SpeedyWeather.jl is currently developed. Some things work, some don't. Stay tuned.
The main interface to SpeedyWeather.jl is 

```julia
julia> using SpeedyWeather
julia> run_speedy(ShallowWater, n_days=30, trunc=63, Grid=OctahedralGaussianGrid, output=true)
Weather is speedy run 1: 100%|███████████████████████████████████| Time: 0:00:04 (1498.70 years/day)
```

Hurray! In 4 seconds we just simulated 30 days of the Earth's atmosphere at a speed of 1500 years per day.
This simulation used a T63 spectral resolution on an octahedral Gaussian grid (~165km resolution) solving
the shallow water equations. The arguments for `run_speedy` are described in
[`src/default_parameters.jl`](https://github.com/milankl/SpeedyWeather.jl/blob/main/src/default_parameters.jl).

To run the primitive equation dry core you can for example do

```julia
julia> run_speedy(Float32,PrimitiveDryCore,physics=true,diffusion=HyperDiffusion(power=2))
Weather is speedy: 100%|███████████████████████████████████| Time: 0:00:03 (753.71 years/day)
```

The arguments here use single precision (Float32), enable physics (Held-Suarez forcing and boundary layer drag) and use biharmonic diffusion (default `power=4`).

## History

SpeedyWeather.jl is a reinvention of the atmospheric general circulation model
[SPEEDY](http://users.ictp.it/~kucharsk/speedy-net.html) in Julia. While conceptually the same model,
it is entirely restructured, features have been added, changed and removed, such that only the numerical
schemes share similarities. Speedy's dynamical core has an obscure history: Large parts were written by Isaac Held
at GFDL in/before the 90ies with an unknown amount of contributions/modifications from Steve Jewson (Oxford) in the 90ies.
The physical parametrizations were then added by Franco Molteni, Fred Kucharski, and Martin P. King in Fortran77.
Speedy was then translated to Fortran90 by Sam Hatfield in [speedy.f90](https://github.com/samhatfield/speedy.f90).
SpeedyWeather.jl is then adopted from [first translations to Julia](https://github.com/samhatfield/speedy.jl) by Sam Hatfield.

## Submodules

SpeedyWeather.jl defines several submodules that are technically stand-alone (with dependencies) but aren't separated
out to their own packages for now

- __LowerTriangularMatrices__, a module that defines `LowerTriangularMatrix` (among others) used for the spherical harmonic coefficients
- __RingGrids__, a module that defines several iso-latitude ring-based spherical grids (like the FullGaussianGrid or the HEALPixGrid) and interpolations between them
- __SpeedyTransforms__, a module that defines the spherical harmonic transform between spectral space (for which LowerTriangularMatrices is used) and grid-point space (as defined by RingGrids).

These modules can also be used independently of SpeedyWeather like so
```julia
julia> import SpeedyWeather: LowerTriangularMatrices, RingGrids, SpeedyTransforms
```

## Installation

SpeedyWeather.jl is registered in Julia's registry, so open the package manager with `]` and
```julia
(@v1.8) pkg> add SpeedyWeather
```
which will install the latest release and all dependencies automatically. The very latest version is installed with
```julia
(@v1.8) pkg> add https://github.com/SpeedyWeather/SpeedyWeather.jl#main
```
which pulls directly from the `#main` branch. Please use the current minor version of Julia,
compatibilities with older versions are not guaranteed.
