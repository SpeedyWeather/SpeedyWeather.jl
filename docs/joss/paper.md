---
title: 'SpeedyWeather.jl: Reinventing atmospheric general circulation models towards interactivity, extensibility and composability'

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
    affiliation: 7

  - name: Simone Silvestri
    affiliation: 1

  - name: Gregory L Wagner
    orcid: 0000-0001-5317-2445
    affiliation: 1

  - name: Alistair White
    affiliation: "3,4"

  - name: Sam Hatfield
    affiliation: 6

  - name: David Meyer
    affiliation: 6

  - name: Tom Kimpson
    affiliation: "2,8"

  - name: Navid C Constantinou
    orcid: 0000-0002-8149-4094
    affiliation: 9

  - name: Chris Hill
    affiliation: 1
    
affiliations:
 - name: Massachusetts Institute of Technology, Cambridge, MA, USA
   index: 1
 - name: University of Oxford, UK
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
 - name: University of Melbourne, Australia
   index: 8
 - name: Australian National University, Canberra, Australia
   index: 9
   
date: 14 September 2023
bibliography: paper.bib

---


# Summary

SpeedyWeather.jl is a library to simulate and analyse the global atmospheric
circulation on the sphere. It implements several 2D and 3D
models to solve the primitive equations with and without humidity (\autoref{fig:primitive}),
the shallow water equations (\autoref{fig:swm}), or the barotropic vorticity equations
with spherical harmonics. Several simple parameterizations for unresolved physical processes
such as precipitation or the boundary layer are implemented, and new ones can
be externally defined and passed as an argument to the model constructor.
SpeedyWeather.jl is an intermediate-complexity general circulation model [@Kucharski2013]
and research playground with an (almost) everything-flexible attitude.
It can be thought of as a conceptual reinvention of the Fortran SPEEDY model [@Molteni2003]
in the Julia programming language [@Bezanson2017].

SpeedyWeather.jl internally uses three sub-modules `SpeedyTransforms`, `RingGrids`, and
`LowerTriangularMatrices` to perform spherical harmonic transforms and interpolations
between various grids and the spectral space. `RingGrids` discretize the sphere
on iso-latitude rings and the spectral space is defined by the `LowerTriangularMatrices`
of the spherical harmonic coefficients. These three modules are independently usable
and therefore make SpeedyWeather.jl, beyond its main purpose of simulating the weather,
also a library for the analysis of gridded data on the sphere.
Running and analysing simulations can be interactively combined, enhancing user
experience and productivity.

The user interface of SpeedyWeather.jl is heavily influenced by
the Julia ocean model Oceananigans.jl [@Ramadhan2020].
A monolithic interface based on parameter files is avoided in favor of a
library-style interface in which users write code to run models rather than
merely supplying parameters and input arrays.
A model is created bottom-up by first defining the discretization
and any non-default model components with its respective parameters.
All components are then collected into a single model object, which, once
initialized, returns a simulation object that contains the entire model state,
work arrays and parameters, that can be run, analysed or changed.
While these steps can be written into a script for reproducibility,
the same steps can be executed and interacted with one-by-one in
Julia's read-evaluate-print loop (REPL). We thereby reach an interactivity
far beyond a monolithic interface that is limited to the options provided.

To be extensible and composable with new
model components, SpeedyWeather.jl relies on Julia's multiple dispatch
programming paradigm [@Bezanson2017]. Every model component
is defined as a new type. For example, to define a new way how to calculate
the precipitation due to the physical process of large-scale condensation,
one would define `MyCondensation` as a new subtype of `AbstractCondensation`.
One then only needs to extend the `initialize!` and `condensation!`
functions for this new type. Passing on `condensation = MyCondensation()`
to the model constructor then implements this new model component without
the need to branch off or overwrite existing model components.
Conceptually similar scientific modelling paradigms have been very successful
in the Python-based generic partial differential equation solver Dedalus [@Burns2020]
and the Julia ocean model Oceananigans.jl [@Ramadhan2020].

![Surface temperature simulated with the primitive equation model in SpeedyWeather.jl.
(Figure will be updated) \label{fig:primitive}](primitive.jpg)

The dynamical core of SpeedyWeather.jl uses established numerics
[@Bourke1972; @Hoskins1975; @Simmons1978; @Simmons1981],
widely adopted in numerical weather prediction. It is based on the spherical
harmonic transform with a leapfrog-based semi-implicit time integration [@Hoskins1975]
and a Robert-Asselin-Williams filter [@Williams2011; @Amezcua2011].
The spherical harmonic transform is grid-flexible. Any iso-latitude ring-based
grid can be used and new grids can be externally defined and passed in
as an argument. Many grids are already implemented: the conventional
Gaussian grid, a regular longitude-latitude grid, 
the octahedral Gaussian grid [@Malardel2016], the octahedral
Clenshaw-Curtis grid [@Hotta2018], and the HEALPix grid [@Gorski2005].
Both SpeedyWeather.jl and its spherical harmonic transform SpeedyTransforms are also
number format-flexible. 32-bit single-precision floating-point numbers
(Float32) are the default as adopted by other modelling efforts [@Vana2017, @Nakano2018],
but Float64 and other custom number formats can be used with a single
code basis. Julia will compile to the choice of the number format, the grid,
and and other model components just-in-time. A simple parallelisation
across vertical layers is supported with Julia's multi threading.

# Statement of need

SpeedyWeather.jl is a fresh approach to atmospheric
models that have been very influential in many areas of scientific 
and high-performance computing as well as climate change mitigation and adaptation.
Most weather, ocean and climate models are written in Fortran and have been developed over
decades. From this tradition follows a specific programming style and
associated user interface. Running a simulation in Fortran and analysing the
data in Python makes it virtually impossible to interact with various model
components interactively. Furthermore, data-driven climate modelling
[@Rasp2018,Schneider2023], which replaces existing model components with machine learning
is difficult due to the lack of established deep learning frameworks in Fortran [@Innes2019].
Let alone online learning, which trains a neural network-based component together
with the rest of the model, accounting for interactions between components.
Gradients, necessary to optimize training, can be computed 
with automatic differentiation [@Moses2020], but only if differentiable functions
in a coherent language framework are provided.

We hope to provide with SpeedyWeather.jl a first test platform for data-driven
atmospheric modelling and in general an interactive model that makes difficult
problems easy to simulate. Climate models that are user-friendly, trainable,
but also easily extensible will suddenly make many complex
research ideas possible.

![Relative vorticity simulated with the shallow water model in SpeedyWeather.jl.
The simulation used a spectral resolution of T1023 (about 20km) and Float32
arithmetic on an octahedral Clenshaw-Curtis grid [@Hotta2018]. Relative vorticity
is visualised with Matplotlib [@Hunter2007] and Cartopy [@Cartopy] using a
transparent-to-white colormap to mimic the appearance of clouds. Underlain is
NASA's blue marble from June 2004. \label{fig:swm}](swm.png)

# Acknowledgements

We acknowledge contributions from Mosè Giordano, Valentin Churavy, and Pietro Monticone
who have also committed to the SpeedyWeather.jl repository, and the wider Julia community
for help and support. MK acknowledges funding from the 
National Science Foundation (Chris please add). MK and TK acknowledge funding from the European Research Council
under the European Union’s Horizon 2020 research and innovation programme for the ITHACA grant (no. 741112).
NCC acknowledges support by the Australian Research Council DECRA Fellowship DE210100749.

# References
