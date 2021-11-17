# SpeedyWeather.jl
[![CI](https://github.com/milankl/SpeedyWeather.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/milankl/SpeedyWeather.jl/actions/workflows/CI.yml)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://milankl.github.io/SpeedyWeather.jl/dev)

The Julia-version of the atmospheric general circulation model [SPEEDY](http://users.ictp.it/~kucharsk/speedy-net.html),
originally written by Fred Kucharski, Franco Molteni and Martin P. King in Fortran77. Then translated to Fortran90 by
Sam Hatfield in [speedy.f90](https://github.com/samhatfield/speedy.f90). SpeedyWeather.jl is then adopted from
[first translations to Julia](https://github.com/samhatfield/speedy.jl) by Sam Hatfield.

Requires: Julia 1.4

## Functionality

For an overview of the functionality and explanation see the
[documentation](https://milankl.github.io/SpeedyWeather.jl/dev).
For a mathematical description of the model see the [original documentation](http://users.ictp.it/~kucharsk/speedy_description/km_ver41_appendixA.pdf) by Molteni and Kucharski.


## Installation

SpeedyWeather.jl is not yet registered, so open the package manager with `]` and
```julia
(@v1.6) pkg> add https://github.com/milankl/SpeedyWeather.jl
```
which will install the `main` branch and all dependencies automatically.
