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

SpeedyWeather.jl is currently developed. The spherical harmonic transform is already implemented
```julia
julia> using SpeedyWeather
julia> alms = zeros(ComplexF64,3,3)    # spectral coefficients
julia> alms[2,2] = 1                   # only l=1,m=1 Legendre polynomial
julia> map = gridded(alms)             # convert to grid space
8×4 Matrix{Float64}:
 -0.324541  -0.600363  -0.600363  -0.324541
 -0.134429  -0.248678  -0.248678  -0.134429
  0.134429   0.248678   0.248678   0.134429
  0.324541   0.600363   0.600363   0.324541
  0.324541   0.600363   0.600363   0.324541
  0.134429   0.248678   0.248678   0.134429
 -0.134429  -0.248678  -0.248678  -0.134429
 -0.324541  -0.600363  -0.600363  -0.324541
 
julia> spectral(map)                   # back to spectral space
3×3 Matrix{ComplexF64}:
 0.0+0.0im  0.0+0.0im          0.0+0.0im
 0.0+0.0im  1.0+3.60727e-17im  0.0+0.0im
 0.0+0.0im  0.0+0.0im          0.0+0.0im
```
And we have succesfully reobtained the `l=1,m=1` Legendre polynomial.

## History

SpeedyWeather.jl is the Julia version of the atmospheric general circulation model [SPEEDY](http://users.ictp.it/~kucharsk/speedy-net.html),
originally written by Fred Kucharski, Franco Molteni and Martin P. King in Fortran77. Then translated to Fortran90 by
Sam Hatfield in [speedy.f90](https://github.com/samhatfield/speedy.f90). SpeedyWeather.jl is then adopted from
[first translations to Julia](https://github.com/samhatfield/speedy.jl) by Sam Hatfield, but entirely restructured
and only the algorithms are shared with the original Fortran versions.

## Installation

SpeedyWeather.jl is not yet registered, so open the package manager with `]` and
```julia
(@v1.7) pkg> add https://github.com/milankl/SpeedyWeather.jl
```
which will install the `main` branch and all dependencies automatically.
