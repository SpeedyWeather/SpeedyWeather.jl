# SpeedyWeather.jl
[![CI](https://github.com/milankl/SpeedyWeather.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/milankl/SpeedyWeather.jl/actions/workflows/CI.yml)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://milankl.github.io/SpeedyWeather.jl/dev)

The focus of SpeedyWeather.jl is to develop a global, yet simple, atmospheric model, that can run at various levels of
precision (16, 32 and 64-bit) on different architectures (x86 and ARM, currently planned, GPUs maybe in the future).
Additionally, the model is written in an entirely number format-flexible way, such that any custom number format can
be used and Julia will compile to the format automatically. In contrast to the original speedy, some adjustments to the
resolution (horizontal and vertical) are also planned with some simple parallelism to run SpeedyWeather.jl efficiently
on small clusters.

## History

SpeedyWeather.jl is the Julia version of the atmospheric general circulation model [SPEEDY](http://users.ictp.it/~kucharsk/speedy-net.html),
originally written by Fred Kucharski, Franco Molteni and Martin P. King in Fortran77. Then translated to Fortran90 by
Sam Hatfield in [speedy.f90](https://github.com/samhatfield/speedy.f90). SpeedyWeather.jl is then adopted from
[first translations to Julia](https://github.com/samhatfield/speedy.jl) by Sam Hatfield, but entirely restructured
and only the algorithms are shared with the original Fortran versions.

## Functionality

For an overview of the functionality and explanation see the
[documentation](https://milankl.github.io/SpeedyWeather.jl/dev).
For a mathematical description of the model see the
[original documentation](http://users.ictp.it/~kucharsk/speedy_description/km_ver41_appendixA.pdf) by Molteni and Kucharski.

## Installation

SpeedyWeather.jl is not yet registered, so open the package manager with `]` and
```julia
(@v1.7) pkg> add https://github.com/milankl/SpeedyWeather.jl
```
which will install the `main` branch and all dependencies automatically.
