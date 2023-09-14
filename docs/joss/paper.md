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
    affiliation: "3,4"

  - name: Daisuke Hotta
    affiliation: "5,6"

  - name: Justin Willmert
    affiliation: 7

  - name: Simone Silvestri
    affiliation: 1

  - name: Gregory L Wagner
    affiliation: 1

  - name: Alistair White
    affiliation: "3,4"

  - name: Sam Hatfield
    affiliation: 8

  - name: David Meyer
    affiliation: 8

  - name: Tom Kimpson
    affiliation: "2,9"

  - name: Navid C Constantinou
    affiliation: 10

  - name: Chris Hill
    affiliation: 1
    
affiliations:
 - name: Massachusetts Institute of Technology, Cambridge, MA, USA
   index: 1
 - name: University of Oxford, UK
   index: 2
 - name: Technical University Munich, Germany
   index: 3
 - name: Potsdam Institute for Climate Impact Research, Germany
   index: 4
 - name: Japan Meteorological Agency, Tsukuba, Japan
   index: 5
 - name: University of Reading, UK
   index: 6
 - name: University of Minnesota, Minneapolis, MN, USA
   index: 7
 - name: European Centre for Medium-Range Weather Forecasts, Reading, UK
   index: 8
 - name: University of Melbourne, Melbourne, Australia
   index: 9
 - name: Australian National University, Canberra, Australia
   index: 10
   
date: 14 September 2023
bibliography: paper.bib

---

# Summary

SpeedyWeather.jl is a library to simulate and analyse the global atmospheric
circulation on the sphere. It implements several 2D and 3D
models to solve the primitive equations with and without humidity,
the shallow water equations (\autoref{fig:swm}), or the barotropic vorticity equations
with spherical harmonics.

The user interface of SpeedyWeather.jl is heavily influenced by
Oceananigans.jl [@Ramadhan2020]. A monolithic interface is deliberately avoided

To be extendible and composable, SpeedyWeather.jl relies on multiple dispatch
from the Julia programming language [@Bezanson2017]. Every model component
is defined as a new type `SomeComponent <: AbstractComponent`, i.e. subtype of
an abstract super type `AbstractComponent`.
Extending SpeedyWeather.jl can therefore easily achieved by defining a new
`OtherComponent <: AbstractComponent` 

The dynamical core of SpeedyWeather.jl uses established numerics
[@Bourke1972; @Hoskins1975; @Simmons1978; @Simmons1981],
widely adapted in numerical weather prediction. It is based on the spherical
harmonic transform with a leapfrog-based semi-implicit time integration [@Hoskins1975]
and a Robert-Asselin-Williams filter [@Williams2011; @Amezcua2011].
The spherical harmonic transform is grid-flexible. Any iso-latitude ring-based
grid can be used and new grids can be externally defined and passed on
as argument. Many grids are already implemented: The conventional
Gaussian grid, a regular longitude-latitude grid, 
the octahedral Gaussian grid [@Malardel2016], the octahedral
Clenshaw-Curtis grid [@Hotta2018], and the HEALPix grid [@Gorski2005].
Both SpeedyWeather.jl and its spherical harmonic transform are also
number format-flexible. 32-bit single-precision floating-point numbers
(Float32) are the default as adapted by other modelling efforts [@Vana2017],
but Float64 and other custom number formats can be used with a single
code basis. Julia will compile to the choice of the number format
and grid (and other model components) just-in-time.

# Statement of need

Most weather and climate models are written in Fortran and have been developed over
decades. From this tradition follows a specific programming style and
associated user interface. SpeedyWeather.jl is a fresh approach to atmospheric
models that have been very influential in many areas of scientific 
and high-performance computing as well as climate change mitigation and adaptation.

![Shallow water simulation with SpeedyWeather.jl.
Relative vorticity with a spectral resolution of T1023 (about 20km) simulated
in Float32 on an octahedral Clenshaw-Curtis grid [@Hotta2018]. Relative vorticity
is visualised with Matplotlib [@Hunter2007] and Cartopy [@Cartopy] using a
transparent-to-white colormap to mimic the appearance of clouds. Underlain is
NASA's blue marble from June 2004. \label{fig:swm}](swm.png)

# Acknowledgements

We acknowledge contributions from Mosè Giordano, Valentin Churavy and Pietro Monticone
who have also committed to the SpeedyWeather.jl repository, and the wider Julia community
help and support. We gratefully acknowledge funding from the 
National Science Foundation (...) and the European Research Council
under the European Union’s Horizon 2020 research and innovation programme
for the ITHACA grant (no. 741112).

# References
