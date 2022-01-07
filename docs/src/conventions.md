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

## Julian conventions

We follow Julia's [style guide](https://docs.julialang.org/en/v1/manual/style-guide/#Style-Guide) and highlight here some important aspects of it.

- __Bang (!) convention__. A function `func` does not change its input arguments, however, `func!` does.
Hence, `func!` is often the in-place version of `func`, avoiding as much memory allocation as possible
and often changing it's first argument, e.g. `func!(out,in)` so that argument `in` is used to calculate
`out` which has been preallocated before function call.
- __Number format flexibility__. Numeric literals such as `2.0` or `1/3` are only used in the model setup
but avoided throughout the code to obtain a fully number format-flexible package using the number format
`NF` as a compile-time variable throughout the code. This often leads to overly specific code whereas
a `Real` would generally suffice. However, this is done to avoid any implicit type conversions.