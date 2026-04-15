# SpeedyWeather.jl documentation

Welcome to the documentation for [SpeedyWeather.jl](https://github.com/milankl/SpeedyWeather.jl)!
SpeedyWeather.jl is a global atmospheric model developed as a research playground
with an everything-flexible attitude as long as it is speedy. Technically it is a climate model with simple,
yet interactive representations of ocean, land and sea-ice. It is easy to use and easy to extend, making 
atmospheric modelling an interactive experience -- in the terminal, in a notebook or conventionally through scripts.

## Overview

SpeedyWeather.jl uses a spherical harmonic transform to simulate
the general circulation of the atmosphere using a vorticity-divergence formulation,
a semi-implicit time integration and simple parameterizations to represent various
climate processes: Convection, clouds, precipitation, radiation, surface fluxes, among others.

SpeedyWeather.jl defines 
- [`BarotropicModel`](@ref barotropic_vorticity_model) for the 2D barotropic vorticity equation
- [`ShallowWaterModel`](@ref shallow_water_model) for the 2D shallow water equations
- [`PrimitiveDryModel`](@ref primitive_equation_model) for the 3D primitive equations without humidity
- [`PrimitiveWetModel`](@ref primitive_equation_model) for the 3D primitive equations with humidity

and solves these equations in spherical coordinates as described in this documentation.

## Vision

Why another model? You may ask. We believe that most currently available are stiff, difficult to use
and extend, and therefore slow down research whereas a modern code in a modern language wouldn't have to.
We decided to use Julia because it combines the best of Fortran and Python: Within a single language
we can interactively run SpeedyWeather but also extend it, inspect its components, evaluate
individual terms of the equations, and analyse and visualise output on the fly.

We do not aim to make SpeedyWeather an atmospheric model similar to the production-ready models used
in weather forecasting, at least not at the cost of our current level of interactivity and ease of
use or extensibility. If someone wants to implement a cloud parameterization that is very complicated
and expensive to run then they are more than encouraged to do so, but it will probably live in
its own repository and we are happy to provide a general interface to do so. But SpeedyWeather's
defaults should be balanced: Physically accurate yet general; as independently as possible from other
components and parameter choices; not too complicated to implement and understand; and computationally cheap.
Finding a good balance is difficult but we try our best. 

## Developers and contributing

The development of  SpeedyWeather.jl is led by 

- [Milan Klöwer](https://github.com/milankl), Oxford
- [Maximilian Gelbrecht](https://github.com/maximilian-gelbrecht), PIK Potsdam

with many [current and past contributors](https://github.com/SpeedyWeather/SpeedyWeather.jl/graphs/contributors).

Open-source lives from large teams of (even occasional) contributors. If you are interested to
fix something, implement something, or just use it and provide feedback you are always welcome.
We are more than happy to guide you, especially when you don't know where to start.
We can point you to the respective code, highlight how everything is connected and tell you
about dos and don'ts. Just express your interest to contribute and we'll be happy to have you!

## Citing

If you use SpeedyWeather.jl in research, teaching, or other activities, we would be grateful 
if you could mention SpeedyWeather.jl and cite our paper in JOSS:

Klöwer et al., (2024). SpeedyWeather.jl: Reinventing atmospheric general circulation models towards interactivity and extensibility. _Journal of Open Source Software_, **9(98)**, 6323, doi:[10.21105/joss.06323](https://doi.org/10.21105/joss.06323).

The bibtex entry for the paper is:

```bibtex
@article{SpeedyWeatherJOSS,
    doi = {10.21105/joss.06323},
    url = {https://doi.org/10.21105/joss.06323},
    year = {2024},
    publisher = {The Open Journal},
    volume = {9},
    number = {98},
    pages = {6323},
    author = {Milan Klöwer and Maximilian Gelbrecht and Daisuke Hotta and Justin Willmert and Simone Silvestri and Gregory L. Wagner and Alistair White and Sam Hatfield and Tom Kimpson and Navid C. Constantinou and Chris Hill},
    title = {{SpeedyWeather.jl: Reinventing atmospheric general circulation models towards interactivity and extensibility}},
    journal = {Journal of Open Source Software}
}
```

## Funding

MK received funding by the European Research Council under Horizon 2020 within the ITHACA project,
grant agreement number 741112 from 2021-2022. From 2022-2024 this project was also funded by the
National Science Foundation NSF. From 2024-2025, the funding is from Schmidt Sciences LLC through
MK's Eric & Wendy Schmidt AI in Science Fellowship. Since 2025, MK's funding is provided through
a NERC Independent Research Fellowship under grant number UKRI191.

All other contributors bring in their own funding from various national and international
grants which is highly acknowledged.
