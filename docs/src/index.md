# SpeedyWeather.jl documentation

Welcome to the documentation for [SpeedyWeather.jl](https://github.com/milankl/SpeedyWeather.jl) a global
atmospheric circulation model with simple parametrizations to represent physical processes such as clouds,
precipitation and radiation.

## Overview

SpeedyWeather.jl is a global spectral model that uses a spherical harmonic transform to perform some calculations
in spectral space (time integration, gradients, linear terms) and some in grid-point space (advection, non-linear terms).
The prognostic variables used are vorticity, divergence, absolute temperature, logarithm of surface
pressure and specific humidity. The time stepping uses a leapfrog scheme with additional filters and a
semi-implicit formulation for gravity waves. The default resolution is T31 (96x48 grid points on a
regular Gaussian grid, about 400km at the Equator) and 8 vertical levels.

Simple parameterizations are used to represent the physical processes convection, large-scale condensation,
clouds, short-wave radiation, long-waves radiation, surface fluxes of momentum and energy, and vertical diffusion.

## Manual outline

See the following pages of the documentation for more details

- [How to run SpeedyWeather.jl](how_to_run_speedy.md)
- [Spherical harmonic transform](spectral_transform.md)
- [Grids](grids.md)
- [Dynamical core](dynamical_core.md)
- [Parametrizations](parametrizations.md)
- [New model setups](new_model_setups.md)
- [Function and type index](functions.md)

and the [original documentation](http://users.ictp.it/~kucharsk/speedy_description/km_ver41_appendixA.pdf)
by Molteni and Kucharski.

## Scope

The focus of SpeedyWeather.jl is to develop a global atmospheric model of intermediate complexity,
that can run at various levels of precision (16, 32 and 64-bit) on different architectures (x86 and ARM, 
GPUs in the future). Additionally, the model is written in an entirely number format-flexible way,
such that any custom number format can be used and Julia will compile to the format automatically.
Similarly, many model components are written in an abstract way to support modularity and extandability.

## History

SpeedyWeather.jl is a Julia implementation of [SPEEDY](http://users.ictp.it/~kucharsk/speedy-net.html),
which is written in Fortran 77. Sam Hatfield
[translated SPEEDY to Fortran 90](https://github.com/samhatfield/speedy.f90) and started the project to
[port it to Julia](https://github.com/samhatfield/speedy.jl). However, we are making an effort to
overhaul the implementation of the mathematical model behind speedy completely and it is unlikely
that a single line of code survived.

## Installation

SpeedyWeather.jl is registered in the Julia Registry. Open Julia's package manager from the REPL with `]`
and `add` the github repository to install SpeedyWeather.jl and all dependencies
```julia
(@v1.8) pkg> add SpeedyWeather
```
which will automatically install the latest release. However, you may want to install directly from the
main branch with
```julia
(@v1.78) pkg> add https://github.com/SpeedyWeather/SpeedyWeather.jl#main
```
other branches than `#main` can be installed by adding `#branch_name` instead.

## Developers

The development of  SpeedyWeather.jl is lead by [Milan Klöwer](https://github.com/milankl) and current and
past contributors include

- [Tom Kimpson](https://github.com/tomkimpson)
- [Alistair White](https://github.com/white-alistair)
- [Maximilian Gelbrecht](https://github.com/maximilian-gelbrecht)
- [David Meyer](https://github.com/dmey)
- [Daisuke Hotta](https://github.com/hottad)

Any contributions are always welcome!

## Funding

Contributors received funding by the European Research Council under Horizon 2020 within the ITHACA project,
grant agreement number 741112 from 2021-2022.