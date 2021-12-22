# Time integration

SpeedyWeather.jl uses a leapfrog time scheme with a Robert's and William's filter
to dampen the computational mode and achieve 3rd order accuracy.

## Oscillation equation

```math
\frac{dF}{dt} = i\omega F
```


## Implementation details
