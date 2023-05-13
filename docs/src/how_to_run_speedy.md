# How to run SpeedyWeather.jl

The simplest way to run SpeedyWeather.jl with default parameters is

```julia
using SpeedyWeather
run_speedy()
```

Hooray, you have just simulated the Earth's atmosphere. Parameters, their meanings and
defaults are documented in [`Parameters`](@ref). For example, if you want to run
the primitive equation dry core (no humidity) simulation in double precision (`Float64`),
at higher resolution (`trunc`, the triangular spectral truncation), slow down the rotation
of the Earth (`rotation in ``s^{-1}``), and create some netCDF ouput, do

```julia
run_speedy(Float64,PrimitiveDryCore,trunc=42,planet=Earth(rotation=1e-5),output=true)
```

If provided, the number format has to be the first argument, the model (`Barotropic`, `ShallowWater`,
`PrimitiveDryCore`, `PrimitiveWetCore` are available) the second, and all other arguments are keyword
arguments.

## The `run_speedy` interface

```@docs
run_speedy
```

## The `initialize_speedy` interface

```@docs
initialize_speedy
```
