# SpeedyWeather.jl <img src="https://github.com/user-attachments/assets/977f5f46-ccd3-49d8-950a-8b619df863c3" width="100" />


[![CI](https://github.com/SpeedyWeather/SpeedyWeather.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/SpeedyWeather/SpeedyWeather.jl/actions/workflows/CI.yml) 
[![status](https://joss.theoj.org/papers/515c81a4d6a69e31cc71ded65ac9c36a/status.svg)](https://joss.theoj.org/papers/515c81a4d6a69e31cc71ded65ac9c36a)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6510139.svg)](https://doi.org/10.5281/zenodo.6510139)  
[![docs](https://img.shields.io/badge/documentation-latest_release-blue.svg)](https://speedyweather.github.io/SpeedyWeatherDocumentation/stable/)
[![docs](https://img.shields.io/badge/documentation-main-blue.svg)](https://speedyweather.github.io/SpeedyWeatherDocumentation/dev/)

SpeedyWeather.jl is a global atmospheric model with simple physics developed as a research playground
with an everything-flexible attitude as long as it is speedy. It is easy to use and easy to extend, making 
atmospheric modelling an interactive experience -- in the terminal, in a notebook or conventionally through scripts.
With minimal code redundancies it supports

**Dynamics and physics**
- Different physical equations (barotropic vorticity, shallow water, primitive equations, with and without humidity)
- Particle advection in 2D for all equations
- Tracer advection in 2D/3D that can be added, deleted, (de)activated anytime
- Physics parameterizations for convection, precipitation, boundary layer, etc.

**Numerics and computing**
- Different spatial grids (full and octahedral grids, Gaussian and Clenshaw-Curtis, HEALPix, OctaHEALPix)
- Different resolutions (T31 to T1023 and higher, i.e. 400km to 10km using linear, quadratic or cubic truncation)
- Different arithmetics: Float32 (default), Float64, and (experimental) BFloat16, stochastic rounding
- a very fast and flexible spherical harmonics transform library SpeedyTransforms

**User interface**
- Data visualisation: 2D, 3D, interactive (you can zoom and rotate!) powered by Makie
- Extensibility: New model components (incl. parameterizations) can be externally defined
- Modularity: Models are constructed from its components, non-defaults are passed on as argument
- Interactivity: SpeedyWeather.jl runs in a notebook or in the REPL as well as from scripts
- Callbacks can be used to inject any piece of code after every time step, e.g. custom output, event handling, changing the model while it's running

and Julia will compile to these choices just-in-time.

For an overview of the functionality and explanation see the
[documentation](https://speedyweather.github.io/SpeedyWeatherDocumentation/dev).
You are always encouraged to [raise an issue](https://github.com/SpeedyWeather/SpeedyWeather.jl/issues)
(even it is not actually an issue but an idea, a suggestion or really anything)
describing what you'd like to use SpeedyWeather for. We're keen to help!

## Vision and roadmap

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

This means in practice, that while SpeedyWeather is currently developed, many more physical processes
and other features will be implemented. On our TODO is

- A (somewhat) realistic radiation scheme with a daily cycle, depending on clouds and humidity
- Longwave radiation that depends on (global) CO2 concentrations to represent climate change
- Slab ocean and a (seasonal cycle) sea ice interacting with radiation
- Exoplanet support
- 3D particle advection
- single GPU support to accelerate medium to high resolution simulations
- differentiability with Enzyme

## Contributing

Open-source lives from large teams of (even occasional) contributors. If you are interested to
fix something, implement something, or just use it and provide feedback you are always welcome.
We are more than happy to guide you, especially when you don't know where to start.
We can point you to the respective code, highlight how everything is connected and tell you
about dos and don'ts. Just express your interest to contribute and we'll be happy to have you.

## Example use

For a more comprehensive tutorial with several examples, see
[Examples](https://speedyweather.github.io/SpeedyWeatherDocumentation/dev/examples_2D/) in the documentation.
The basic interface to SpeedyWeather.jl consist of 4 steps: define the grid,
construct the model, initialize, run

```julia
spectral_grid = SpectralGrid(trunc=31, nlayers=8)   # define resolution
model = PrimitiveWetModel(spectral_grid)            # construct model
simulation = initialize!(model)                     # initialize all model components
run!(simulation, period=Day(10), output=true)       # aaaand action!
```
and you will see

<img src="https://github.com/SpeedyWeather/SpeedyWeather.jl/assets/25530332/a04fbb10-1cc1-4f77-93f2-7bdf047f277d" width="450"><br>

Hurray🥳 In a few seconds seconds we just simulated 10 days of the Earth's atmosphere at a speed of 440 years per day.
This simulation used a T31 spectral resolution on an
[octahedral Gaussian grid](https://speedyweather.github.io/SpeedyWeatherDocumentation/dev/grids/#Implemented-grids)
(~400km resolution) solving the primitive equations on 8 vertical levels.
The [UnicodePlot](https://github.com/JuliaPlots/UnicodePlots.jl) will give
you a snapshot of surface vorticity at the last time step. The plotted resolution is not representative,
but allows a quick check of what has been simulated.
The [NetCDF output](https://speedyweather.github.io/SpeedyWeatherDocumentation/dev/output/) is independent of the unicode plot.

More examples in the [How to run SpeedyWeather](https://speedyweather.github.io/SpeedyWeatherDocumentation/dev/how_to_run_speedy/)
section of the [documentation](https://speedyweather.github.io/SpeedyWeatherDocumentation/dev).

## Gallery

Specific humidity in the primitive equation model simulated at T340 spectral resolution (about 40km)
with 16 vertical levels (shown here is level 15, just above the surface) on the octahedral Gaussian grid
computed in single precision multi-threaded on 16 CPUs. With convection, large-scale condensation, surface fluxes
and some simplified radiation (the daily cycle is visible)

https://github.com/SpeedyWeather/SpeedyWeather.jl/assets/25530332/614f04cf-5080-4c89-9fd6-35efd54103a7

Relative vorticity in the shallow water model, simulated at T1023 spectral resolution (about 10km) on an
[octahedral Clenshaw-Curtis grid](https://speedyweather.github.io/SpeedyWeatherDocumentation/dev/grids/#Implemented-grids)
with more than 4 million grid points

https://user-images.githubusercontent.com/25530332/190443050-d5b8d093-86c0-46c9-b515-8420059ac8dc.mp4

Surface temperature in the primitive equation model without surface fluxes or radiation
at T511 (~20km resolution) and 31 vertical levels. The simulation was multi-threaded in Float32 (single precision).

https://github.com/SpeedyWeather/SpeedyWeather.jl/assets/25530332/95897b82-9b81-4980-934b-cfdcf4d5a4b0

SpeedyWeather.jl can also solve the 2D barotropic vorticity equations on the sphere.
Here, we use Float32 (single precision) at a resolution of T340 (40km) on
an [octahedral Gaussian grid](https://speedyweather.github.io/SpeedyWeatherDocumentation/dev/grids/#Implemented-grids).
Forcing is a stochastic stirring on northern hemisphere mid-latitudes following Barnes and Hartmann, 2011.
Map projection is orthographic centred on the north pole.

https://github.com/SpeedyWeather/SpeedyWeather.jl/assets/25530332/3d7fccd5-b66d-42e3-9f73-64dcf21d00ee

Here, SpeedyWeather.jl simulates atmospheric gravity waves, initialised randomly interacting with orography
over a period of 2 days. Each frame is one time step, solved with a centred semi-implicit scheme that
resolves gravity waves with a timestep of CFL=1.2-1.4 despite a single-stage RAW-filtered leapfrog integration.

https://github.com/SpeedyWeather/SpeedyWeather.jl/assets/25530332/510c38c7-12cb-42d5-b905-c66b4eaa514d

Advection of 5000 particles with wind in the lower-most layer of the primitive equation model at
T85 (150km) resolution and 8 vertical layers.

https://github.com/SpeedyWeather/SpeedyWeather.jl/assets/25530332/a6192374-24d9-4065-9fcc-8b719190472f

## Data visualisation

Difficult to plot spherical data? SpeedyWeather also includes extensions for Makie and GeoMakie
making it supereasy to create plots and interactively investigate a variables from a simulation.
Two examples (screen recording those makes it a bit laggy, it's pretty smooth otherwise): Humidity plotted
on a 50km HEALPix grid 

https://github.com/user-attachments/assets/b02b31eb-e139-4193-89d1-7e277a2af5cc

or the visualising cell centres and faces of the OctaminimalGaussianGrid 

https://github.com/user-attachments/assets/6dfa212a-c5dc-4c54-b274-7755d5baf15c

## History

SpeedyWeather.jl started off as a reinvention of the atmospheric general circulation model
[SPEEDY](http://users.ictp.it/~kucharsk/speedy-net.html) in Julia. While conceptually a similar model,
it is entirely restructured, features have been added, changed and removed, such that only some of the numerical
schemes share similarities. Fortran SPEEDY's dynamical core has an obscure history:
Large parts were written by Isaac Held at GFDL in/before the 90ies with an unknown amount of
contributions/modifications from Steve Jewson (Oxford) in the 90ies.
The physical parametrizations were then added by Franco Molteni, Fred Kucharski, and Martin P. King
afterwards while the model was still written in Fortran77.
Around 2018-19, SPEEDY was then translated to Fortran90 by Sam Hatfield in [speedy.f90](https://github.com/samhatfield/speedy.f90).
SpeedyWeather.jl is then adopted from [first translations to Julia](https://github.com/samhatfield/speedy.jl) by Sam Hatfield.

## Submodules

SpeedyWeather.jl defines several submodules that are technically stand-alone (with dependencies) but aren't separated
out to their own packages for now

- [__RingGrids__](https://speedyweather.github.io/SpeedyWeather.jl/dev/ringgrids/),
a module that defines several iso-latitude ring-based spherical grids (like the FullGaussianGrid or the HEALPixGrid)
and interpolations between them
- [__LowerTriangularMatrices__](https://speedyweather.github.io/SpeedyWeather.jl/dev/lowertriangularmatrices/),
a module that defines `LowerTriangularMatrix` used for the spherical harmonic coefficients
- [__SpeedyTransforms__](https://speedyweather.github.io/SpeedyWeather.jl/dev/speedytransforms/), a module that defines
the spherical harmonic transform between spectral space (for which LowerTriangularMatrices is used) and grid-point space
(as defined by RingGrids).

These modules can also be used independently of SpeedyWeather like so
```julia
julia> using SpeedyWeather: LowerTriangularMatrices, RingGrids, SpeedyTransforms
```
check out their documentation: [RingGrids](https://speedyweather.github.io/SpeedyWeatherDocumentation/dev/ringgrids/),
[LowerTriangularMatrices](https://speedyweather.github.io/SpeedyWeatherDocumentation/dev/lowertriangularmatrices/),
[SpeedyTransforms](https://speedyweather.github.io/SpeedyWeatherDocumentation/dev/speedytransforms/).

## Installation

SpeedyWeather.jl is registered in Julia's registry, so open the package manager with `]` and
```julia
(@v1.11) pkg> add SpeedyWeather
```
which will install the [latest release]([url](https://github.com/SpeedyWeather/SpeedyWeather.jl/releases))
and all dependencies automatically. For more information see the
[Installation](https://speedyweather.github.io/SpeedyWeatherDocumentation/dev/installation/) in the documentation.
Please use the current minor version of Julia,
compatibilities with older versions are not guaranteed.

## Benchmarks

The primitive equations at 400km resolution with 8 vertical layers can be simulated by
SpeedyWeather.jl at 1800 simulated years per day (SYPD) on a single core of newer CPUs with arm architecture
(M-series MacBooks for example). At that speed, simulating one year takes about 50 seconds
without output. The complex fused-multiply adds of the spectral transform compile efficiently to
the large vector extensions of the newer arm chips in single precision.
Another considerable speed-up comes from the reduced grids minimizing the number of columns for
which expensive parameterizations like convection have to be computed. The parameterizations
take up 40-60% of the total simulation time, depending on the grid. Particularly the
`OctaminimalGaussianGrid`, `OctaHEALPixGrid` and the `HEALPixGrid` are increasingly faster,
at a small accuracy sacrifice of the then inexact spectral transforms. 

On older CPUs, like the Intel CPU MacBooks, the 1800 SYPD drop to about 500-600 SYPD,
which is still 2x faster than Fortran SPEEDY which is reported to reach 240 SYPD.

For an overview of typical simulation speeds a user can expect under different model setups see
[Benchmarks](https://github.com/SpeedyWeather/SpeedyWeather.jl/blob/main/benchmark).

## Citing

If you use SpeedyWeather.jl in research, teaching, or other activities, we would be grateful 
if you could mention SpeedyWeather.jl and cite our paper in JOSS:

> Klöwer et al., (2024). SpeedyWeather.jl: Reinventing atmospheric general circulation models towards interactivity and extensibility. Journal of Open Source Software, 9(98), 6323, doi:[10.21105/joss.06323](https://doi.org/10.21105/joss.06323).

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

## Copyright and license

Copyright (c) 2020 Milan Klöwer for SpeedyWeather.jl  
Copyright (c) 2021 The SpeedyWeather.jl Contributors for SpeedyWeather.jl  
Copyright (c) 2022 Fred Kucharski and Franco Molteni for SPEEDY parametrization schemes  

Software licensed under the [MIT License](LICENSE.txt).
