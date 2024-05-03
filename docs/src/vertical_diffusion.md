# Vertical Diffusion

Vertical diffusion in SpeedyWeather.jl is implemented as a Laplacian in
the vertical [Sigma coordinates](@ref) with a diffusion coefficient ``K``
that in general depends on space and time and is flow-aware, meaning
it is recalculated on every time step depending on the vertical stability
of the atmospheric column.

Vertical diffusion can be applied to velocities ``u, v``, temperature ``T``
and specific humidity ``q``.

## Implementations

The following schemes for vertical diffusion are currently implemented

```@example surface_fluxes
using InteractiveUtils # hide
using SpeedyWeather
subtypes(SpeedyWeather.AbstractVerticalDiffusion)
```

`NoVerticalDiffusion` disabled all vertical diffusion,
`BulkRichardsonDiffusion` is explained in the following.

## Laplacian in sigma coordinates

The vertical diffusion of a variable, say ``u``, takes the in
sigma ``\sigma`` coordinates the form

```math
\frac{\partial u}{\partial t} = \frac{\partial}{\partial \sigma} K
\frac{\partial u}{\partial \sigma}
```

as a tendency to ``u`` with no-flux boundary conditions
``\frac{\partial u}{\partial \sigma} = 0``
at ``\sigma = 0`` (top of the atmosphere) and ``\sigma = 1`` (the surface).
That way the diffusion preserves the integral of the variable ``u`` from
``\sigma = 0`` to ``\sigma = 1``.