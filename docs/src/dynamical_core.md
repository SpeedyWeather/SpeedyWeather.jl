# Dynamical core

A mathematical and implementation-specific description of the dynamical core used in SpeedyWeather.jl

## Mathematical background

The [primitive equations](https://en.wikipedia.org/wiki/Primitive_equations) solved by SpeedyWeather.jl are

```math
\begin{aligned}
\partial_t u = ... \\
\partial_t v = ... \\
\partial_t T = ... \\ 
\partial_t Q = ... \\
\end{aligned}
```

more to come

### Implementation details

```julia
using SpeedyWeather

P = Params(T=Float64)
G = GeoSpectral{P.T}(P)
B = Boundaries{P.T}(P,G)

fourier(B.ϕ0trunc,G),G)
```

## Time integration

SpeedyWeather.jl uses a leapfrog time scheme with a Robert's and William's filter
to dampen the computational mode and achieve 3rd order accuracy.

### Oscillation equation

```math
\frac{dF}{dt} = i\omega F
```