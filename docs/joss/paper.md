---
title: 'SpeedyWeather.jl: Reinventing atmospheric general circulation models towards interactivity and extensibility'

tags:
  - Julia
  - weather
  - climate
  - general circulation model
  - spectral
  - spherical harmonic transform

authors:
  - name: Milan Klöwer
    orcid: 0000-0002-3920-4356
    email: milank@mit.edu
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: "1, 2" # (Multiple affiliations must be quoted)

  - name: Maximilian Gelbrecht
    orcid: 0000-0002-0729-6671
    affiliation: "3,4"

  - name: Daisuke Hotta
    orcid: 0000-0003-2287-0608
    affiliation: "5,6"

  - name: Justin Willmert
    orcid: 0000-0002-6452-4693
    affiliation: 7

  - name: Simone Silvestri
    orcid: 0000-0002-7156-946X
    affiliation: 1

  - name: Gregory L Wagner
    orcid: 0000-0001-5317-2445
    affiliation: 1

  - name: Alistair White
    orcid: 0000-0003-3377-6852
    affiliation: "3,4"

  - name: Sam Hatfield
    orcid: 0000-0001-7235-6450
    affiliation: 6

  - name: Tom Kimpson
    orcid: 0000-0002-6542-6032
    affiliation: "2,8"

  - name: Navid C Constantinou
    orcid: 0000-0002-8149-4094
    affiliation: "8, 9"

  - name: Chris Hill
    affiliation: 1
    
affiliations:
 - name: Massachusetts Institute of Technology, Cambridge, MA, USA
   index: 1
 - name: University of Oxford, UK
   index: 2
 - name: Technical University of Munich, Germany
   index: 3
 - name: Potsdam Institute for Climate Impact Research, Germany
   index: 4
 - name: Japan Meteorological Agency, Tsukuba, Japan
   index: 5
 - name: European Centre for Medium-Range Weather Forecasts, Reading, UK
   index: 6
 - name: University of Minnesota, Minneapolis, MN, USA
   index: 7
 - name: University of Melbourne, Parkville, VIC, Australia
   index: 8
 - name: ARC Centre of Excellence for the Weather of the 21st Century, University of Melbourne, Parkville, VIC, Australia
   index: 9
   
date: 7 June 2024
bibliography: paper.bib

---


# Summary

SpeedyWeather.jl is a library to simulate and analyze the global atmospheric
circulation on the sphere. It implements several 2D and 3D
models which solve different sets of equations:

- the primitive equations with and without humidity (\autoref{fig:primitive}),
- the shallow water equations (\autoref{fig:swm}), and
- the barotropic vorticity equation (\autoref{fig:particles}).

The primitive equation model in SpeedyWeather.jl is an
atmospheric general circulation model [@Kucharski2013] with simple parameterizations
for unresolved physical processes including precipitation or boundary layer mixing.
It can be thought of as a conceptual reinvention of the Fortran SPEEDY model [@Molteni2003]
in the Julia programming language [@Bezanson2017]. However, all models here are written in a modular
way to make its components easily extensible. For example, a new parameterization can
be externally defined and passed as an argument to the model constructor.
Operators used inside SpeedyWeather.jl are exposed to the user,
facilitating analysis of the simulation data.
SpeedyWeather.jl is therefore, beyond its main purpose of simulating 
atmospheric motion, also a library for the analysis of gridded data on the sphere.
Running and analyzing simulations can be interactively combined, enhancing user
experience and productivity.

The user interface of SpeedyWeather.jl is heavily influenced by
the Julia ocean model Oceananigans.jl [@Ramadhan2020].
A monolithic interface [@Mazlami2017], controlling most of the model's functionality
through arguments of a single function or through parameter files (often called namelists in Fortran),
is avoided in favor of a library-style interface.
A model is constructed bottom-up by first defining the discretization
and any non-default model components with their respective parameters.
All components are then collected into a single model object which, once
initialized, returns a simulation object. A simulation contains everything,
the model with all parameters as constructed before but also all prognostic and diagnostic variables.
Such a simulation can then be run, but also accessed before and after to analyze or
visualize the current variables, or individual terms of the equations.
One can also adjust some parameters before resuming the simulation.
While these steps can be written into a script for reproducibility,
the same steps can be executed and interacted with one-by-one in
Julia's read-evaluate-print loop (REPL) or in a Jupyter or Pluto notebook.
We thereby achieve an interactivity of a simulation and its various model components
far beyond the options provided in a monolithic interface.
At the same time, defaults, set to well-established test cases, 
enable even inexperienced users to run simulations in just a few lines of code. 

![Surface humidity, air temperature, wind speed and precipitation simulated
with the primitive equation model in SpeedyWeather.jl. Spectral resolution is
T340 (about 40km) on an octahedral Gaussian grid [@Malardel2016] with simple
physics to represent unresolved processes such as surface fluxes including
evaporation, and precipitation due to large-scale condensation and convection.
\label{fig:primitive}](primitive.jpg)

SpeedyWeather.jl relies on Julia's multiple dispatch programming paradigm [@Bezanson2017]
to be extensible with new components including parameterizations, forcing, drag,
or even the grid.
All such supported model components define an abstract type that can be
subtyped to introduce, for example, a new parameterization.
To define a new parameterization for convection in a given vertical column of the atmosphere,
one would define `MyConvection` as a new subtype of `AbstractConvection`.
One then only needs to extend the `initialize!` (executed once during model initialization)
and `convection!` (executed on every time step)
functions for this new type. Passing on `convection = MyConvection()`
to the model constructor then implements this new model component without
the need to branch off or overwrite existing model components.
Conceptually similar scientific modelling paradigms have been very successful
in the Python-based generic partial differential equation solver Dedalus [@Burns2020],
the process-oriented climate model CLIMLAB [@Rose2018],
and the Julia ocean model Oceananigans.jl [@Ramadhan2020].

The dynamical core of SpeedyWeather.jl uses established numerics
[@Bourke1972; @Hoskins1975; @Simmons1978; @Simmons1981],
widely adopted in numerical weather prediction. It is based on the spherical
harmonic transform [@reinecke2013; @stompor2011] with a leapfrog-based semi-implicit
time integration [@Hoskins1975] and a Robert-Asselin-Williams filter [@Williams2011; @Amezcua2011].
The spherical harmonic transform is grid-flexible [@Willmert2020]. Any iso-latitude ring-based
grid can be used and new grids can be externally defined and passed in
as an argument. Many grids are already implemented: the conventional
Gaussian grid, a regular longitude-latitude grid, 
the octahedral Gaussian grid [@Malardel2016], the octahedral
Clenshaw-Curtis grid [@Hotta2018], and the HEALPix grid [@Gorski2005].
Both SpeedyWeather.jl and its spherical harmonic transform are also
number format-flexible. Single-precision floating-point numbers
(Float32) are the default as adopted by other modelling efforts [@Vana2017; @Nakano2018],
but Float64 and other custom number formats can be used with a single
code basis [@Klower2022; @Klower2020].
Julia will compile to the choice of number format, the grid,
and and other model components just-in-time. A simple parallelization
(across vertical layers for the dynamical core, across horizontal grid points
for the parameterizations) is supported by Julia's multithreading.
No distributed-memory parallelization is currently supported,
GPU support is planned.

SpeedyWeather.jl internally uses three sub-modules `RingGrids`,
`LowerTriangularMatrices`, and `SpeedyTransforms`. `RingGrids` is a module that discretizes
the sphere on iso-latitude rings and implements interpolations between various such grids.
`LowerTriangularMatrices` facilitates the implementation of the spherical harmonics by organizing
their coefficients in a lower triangular matrix representation.
`SpeedyTransforms` implements the spectral transform between
the grid-point space as defined by `RingGrids` and the spectral space defined in
`LowerTriangularMatrices`. These three modules are independently usable
and therefore support SpeedyWeather's library-like user interface.
Output is stored as NetCDF files using NCDatasets.jl [@NCDatasets].

![Relative vorticity simulated with the shallow water model in SpeedyWeather.jl.
The simulation used a spectral resolution of T1023 (about 20 km) and Float32
arithmetic on an octahedral Clenshaw-Curtis grid [@Hotta2018]. Relative vorticity
is visualized with Matplotlib [@Hunter2007] and Cartopy [@Cartopy] using a
transparent-to-white colormap to mimic the appearance of clouds. Underlaid is
NASA's blue marble from June 2004. \label{fig:swm}](swm.png)

# Statement of need

SpeedyWeather.jl is a fresh approach to atmospheric models that have been
very influential in many areas of scientific  and high-performance computing
as well as climate change mitigation and adaptation.
Most weather, ocean and climate models are written in Fortran
(e.g. ICON [@ICON], CESM [@CESM], MITgcm [@MITgcm], NEMO [@NEMO]) and have been
developed over decades. From this tradition follows a specific programming
style and associated user interface.
SpeedyWeather.jl aims to overcome the constraints of traditional Fortran-based models.
The modern trend sees simulations in Fortran and data analysis in Python
(e.g. NumPy [@Numpy], Xarray [@Xarray], Dask [@Dask], Matplotlib [@Hunter2007]),
making it virtually impossible to interact with various model components directly.
In SpeedyWeather.jl, interfaces to the model components are exposed to the user.
Furthermore, data-driven climate modelling [@Rasp2018; @Schneider2023],
which replaces existing model components with machine learning,
is more difficult in Fortran due to the lack of established machine learning
frameworks [@Meyer2022a]. 
In Julia, Flux.jl [@Innes2019] is available for machine learning as well as automatic
differentiation with Enzyme [@Moses2020] for gradients-based optimization.

With SpeedyWeather.jl we hope to provide a platform for data-driven
atmospheric modelling and in general an interactive model that makes difficult
problems easy to simulate. Climate models that are user-friendly, trainable,
but also easily extensible will suddenly make many complex
research ideas possible.

![Particle trajectories advected in the barotropic vorticity model. The barotropic
vorticity equations were stochastically stirred at wave numbers 8 to 40 for homogeneous
turbulence on the sphere. The simulation was computed at T340 (about 40km global resolution).
Visualized are 20,000 particles (white dots) with trajectories (colored randomly)
for the previous 6 hours. \label{fig:particles}](particles.jpg)

# Acknowledgements

We acknowledge contributions from David Meyer, Mosè Giordano, Valentin Churavy, and Pietro Monticone
who have also committed to the SpeedyWeather.jl repository, and the wider Julia community
for help and support.
MK acknowledges funding from the National Science Foundation.
MK and TK acknowledge funding from the European Research Council
under the European Union's Horizon 2020 research and innovation programme for the ITHACA grant (no. 741112).
NCC acknowledges support by the Australian Research Council under DECRA Fellowship DE210100749 and the Center of Excellence for the Weather of the 21st Century CE230100012.

# References
