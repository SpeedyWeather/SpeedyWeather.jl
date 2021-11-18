# SpeedyWeather.jl documentation

Welcome to the documentation for [SpeedyWeather.jl](https://github.com/milankl/SpeedyWeather.jl) a global atmospheric
circulation model with simple parameterizations to represent physical processes such as clouds, precipitation and radiation.

## Overview

SpeedyWeather.jl is a spectral model that uses a Fourier and Legendre transform to calculcate tendencies of the prognostic variables
vorticity, divergence, absolute temperature, logarithm of surface pressure and specific humidity. The time stepping uses a leapfrog scheme
with additional filters and a semi-implicit formulation for gravity waves. The default resolution is T30 (96x48 grid points on a
Gaussian grid, about 400km at the Equator) and 8 vertical levels.

Simple parameterizations are used to represent the physical processes convection, large-scale condensation, clouds, short-wave radiation,
long-waves radiation, surface fluxes of momentum and energy, and vertical diffusion.

## Manual outline

See the following pages of the documentation for more details

- [How to run SpeedyWeather.jl](how_to_run_speedy.md)
- [Dynamical core](dynamical_core.md)
- [Parameterizations](parameterizations.md)
- [New model setups](new_model_setups.md)
- [Function and type index](functions.md)

and the [original documentation](http://users.ictp.it/~kucharsk/speedy_description/km_ver41_appendixA.pdf) by Molteni and Kucharski.

## Scope

The focus of SpeedyWeather.jl is to develop a global, yet simple, atmospheric model, that can run at various levels of precision
(16, 32 and 64-bit) on different architectures (x86 and ARM, currently planned, GPUs probably in the future). Additionally, the
model is written in an entirely number format-flexible way, such that any custom number format can be used and Julia will compile
to the format automatically.

## History

SpeedyWeather.jl is a Julia implementation of [SPEEDY](http://users.ictp.it/~kucharsk/speedy-net.html), which is written in Fortran 77.
Sam Hatfield [translated SPEEDY to Fortran 90](https://github.com/samhatfield/speedy.f90) and started the project to port it to Julia in
[first translations to Julia](https://github.com/samhatfield/speedy.jl).

## Installation

SpeedyWeather.jl is not yet registered in the Julia Registry. So at the moment, open Julia's package manager from the REPL with `]` and
`add` the github repository to install SpeedyWeather.jl and all dependencies
```julia
(@v1.6) pkg> add https://github.com/milankl/SpeedyWeather.jl
```
other branches can be installed by adding `#branch_name`, e.g. `add https://github.com/milankl/SpeedyWeather.jl#branch_name`.

## Developers

SpeedyWeather.jl is currently developed by [Milan Klöwer](https://github.com/milankl) and [Tom Kipson](https://github.com/tomkimpson), any contributions are always welcome.

## Funding

This project is funded by the European Research Council under Horizon 2020 within the ITHACA project, grant agreement number 741112.