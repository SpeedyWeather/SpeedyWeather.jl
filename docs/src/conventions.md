# Style and convention guide

In SpeedyWeather.jl we've been following the several conventions that are documented here.

## Variable naming

The prognostic variables in spectral space are called

```julia
    vor         # Vorticity of horizontal wind field
    div         # Divergence of horizontal wind field
    temp        # Absolute temperature [K]
    pres_surf   # Logarithm of surface pressure [log(Pa)]
    humid       # Specific humidity [g/kg]
```

their transforms into grid-point space get a `_grid` suffix, their tendencies a `_tend` suffix. Further derived diagnostic dynamic variables are

```julia
    u
    v
    geopot
    ...
```

## Preallocation

All arrays representing variables are preallocated and grouped into two structs

```julia
    Prog::PrognosticVariables
    Diag::DiagnosticVariables
```

The `Diag` struct contains further structs which represent the grid-point transformations of the prognostic variables and their tendencies.

```julia
    gridvars::GridVariables
    tendencies::Tendencies
    ...
```

Constant arrays are grouped into several structs

```julia
Boundaries
```