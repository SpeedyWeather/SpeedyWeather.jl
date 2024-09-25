# SpeedyWeather.jl documentation

Welcome to the documentation for [SpeedyWeather.jl](https://github.com/milankl/SpeedyWeather.jl) a global
atmospheric circulation model with simple parametrizations to represent physical processes such as clouds,
precipitation and radiation. SpeedyWeather in general is more a library than just a model as it exposes
most of its internal functions to the user such that simulations and analysis can be interactively
combined. Its user interface is built in a very modular way such that new components can be easily
defined and integrated into SpeedyWeather.

## Overview

SpeedyWeather.jl is uses a spherical harmonic transform to simulate
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

The development of  SpeedyWeather.jl is lead by [Milan Klöwer](https://github.com/milankl) and
[current and past contributors](https://github.com/SpeedyWeather/SpeedyWeather.jl/graphs/contributors) include

- [Tom Kimpson](https://github.com/tomkimpson)
- [Alistair White](https://github.com/white-alistair)
- [Maximilian Gelbrecht](https://github.com/maximilian-gelbrecht)
- [David Meyer](https://github.com/dmey)
- [Daisuke Hotta](https://github.com/hottad)
- [Navid Constantinou](https://github.com/navidcy)
- [Simone Silvestri](https://github.com/simone-silvestri)

(Apologies if you've recently started contributing but this isn't reflected here yet, create a pull request!)
Any contributions are always welcome!

Open-source lives from large teams of (even occasional) contributors. If you are interested to
fix something, implement something, or just use it and provide feedback you are always welcome.
We are more than happy to guide you, especially when you don't know where to start.
We can point you to the respective code, highlight how everything is connected and tell you
about dos and don'ts. Just express your interest to contribute and we'll be happy to have you.

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
grant agreement number 741112 from 2021-2022. From 2022-2024 this project is also funded by the
National Science Foundation NSF. Since 2024, the main funding is from Schmidt Sciences through
a Eric & Wendy Schmidt AI in Science Fellowship.
