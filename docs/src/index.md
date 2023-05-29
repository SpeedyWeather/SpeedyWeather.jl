# SpeedyWeather.jl documentation

Welcome to the documentation for [SpeedyWeather.jl](https://github.com/milankl/SpeedyWeather.jl) a global
atmospheric circulation model with simple parametrizations to represent physical processes such as clouds,
precipitation and radiation.

## Overview

SpeedyWeather.jl is a global spectral model that uses a spherical harmonic transform to perform some calculations
in spectral space (time integration, gradients, linear terms) and some in grid-point space (advection, non-linear terms,
parameterizations).
The prognostic variables used are vorticity, divergence, absolute temperature, logarithm of surface
pressure and specific humidity. The time stepping uses a leapfrog scheme with additional filters and a
semi-implicit formulation for gravity waves. The default resolution is T31 (96x48 grid points on a
regular Gaussian grid, about 400km at the Equator) and 8 vertical levels.

Simple parameterizations are used to represent the physical processes convection, large-scale condensation,
clouds, short-wave radiation, long-waves radiation, surface fluxes of momentum and energy, and vertical diffusion.

## Manual outline

See the following pages of the documentation for more details

- [Installation](installation.md)
- [How to run SpeedyWeather.jl](how_to_run_speedy.md)
- [Spherical harmonic transform](spectral_transform.md)
- [Grids](grids.md)
- [Dynamical core](dynamical_core.md)
- [Parametrizations](parametrizations.md)
- [Extending SpeedyWeather](extending.md)

and the submodules

- [RingGrids](@ref) and their interpolation   
- [LowerTriangularMatrices](@ref)   
- [SpeedyTransforms](@ref)

and the [original documentation](http://users.ictp.it/~kucharsk/speedy_description/km_ver41_appendixA.pdf)
by Molteni and Kucharski.

## Developers

The development of  SpeedyWeather.jl is lead by [Milan Klöwer](https://github.com/milankl) and
[current and past contributors](https://github.com/SpeedyWeather/SpeedyWeather.jl/graphs/contributors) include

- [Tom Kimpson](https://github.com/tomkimpson)
- [Alistair White](https://github.com/white-alistair)
- [Maximilian Gelbrecht](https://github.com/maximilian-gelbrecht)
- [David Meyer](https://github.com/dmey)
- [Daisuke Hotta](https://github.com/hottad)
- [Navid Constantinou](https://github.com/navidcy)

Any contributions are always welcome!

## Funding

MK received funding by the European Research Council under Horizon 2020 within the ITHACA project,
grant agreement number 741112 from 2021-2022. Since 2023 this project is also funded by the
National Science Foundation NSF.