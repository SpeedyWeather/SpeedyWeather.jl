# Grids

The spectral transform (the [Spherical Harmonic Transform](@ref)) in SpeedyWeather.jl supports any ring-based
equi-longitude grid. Several grids are already implemented but other can be added. The following pages will
describe an overview of these grids and how they can be used.

#### Ring-based equi-longitude grids

SpeedyWeather.jl's spectral transform currently only supports ring-based equi-longitude grids.
These grids have their grid points located on rings with constant latitude and on rings the
points are equi-spaced in longitude. There is technically no constrain on the spacing of the
latitude rings, but the Legendre transform requires a quadrature to map those to spectral space
and back. Common choices for latitudes are the
[Gaussian latitudes](https://en.wikipedia.org/wiki/Gaussian_grid) which use the Gaussian
quadrature, or equi-angle latitudes (i.e. just regular latitudes but excluding the poles) that
use the [Clenshaw-Curtis quadrature](https://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature).
The longitudes have to be equi-spaced on every ring, which is necessary for the fast Fourier transform,
as one would otherwise need to use a
[non-uniform Fourier transform](https://en.wikipedia.org/wiki/Non-uniform_discrete_Fourier_transform#Nonuniform_fast_Fourier_transform).
In SpeedyWeather.jl the first grid point on any ring can have a longitudinal offset though, for example
by spacing 4 points around the globe at 45˚E, 135˚E, 225˚E, and 315˚E. In this case the offset is 45˚E
as the first point is not at 0˚E.

## Implemented grids

All grids in SpeedyWeather.jl are a subtype of `AbstractGrid`, i.e. `<: AbstractGrid`. We further distinguish
between _full_, and _reduced_ grids. Full grids have the same number of longitude points on every latitude
ring (i.e. points converge towards the poles) and reduced grids reduce the number of points towards the poles
to have them more evenly spread out across the globe. More evenly does not necessarily mean that a grid is
equal-area, meaning that every grid cell covers exactly the same area (although the shape changes).

Currently the following full grids `<: AbstractFullGrid` are implemented

- `FullGaussianGrid`, a full grid with Gaussian latitudes
- `FullClenshawGrid`, a full grid with equi-angle latitudes

and additionally we have `FullHEALPixGrid` and `FullOctaHEALPixGrid` which are the full grid equivalents to the
[HEALPix grid](@ref) and the [OctaHEALPix grid](@ref) discussed below. Full grid equivalent means that they have
the same latitude rings, but no reduction in the number of points per ring towards the poles and no longitude offset.
Other implemented reduced grids are

- `OctahedralGaussianGrid`, a reduced grid with Gaussian latitudes based on an [octahedron](https://en.wikipedia.org/wiki/Octahedron)
- `OctahedralClenshawGrid`, similar but based on equi-angle latitudes
- `HEALPixGrid`, an equal-area grid based on a [dodecahedron](https://en.wikipedia.org/wiki/Rhombic_dodecahedron) with 12 faces
- `OctaHEALPixGrid`, an equal-area grid from the class of HEALPix grids but based on an octahedron.

An overview of these grids is visualised here

![Overview of implemented grids in SpeedyWeather.jl](../img/grids_comparison.png)

Visualised are each grid's grid points (white dots) and grid faces (white lines).
All grids shown have 16 latitude rings on one hemisphere, Equator included.
The total number of grid points is denoted in the top left of every subplot.
The sphere is shaded with grey, orange and turquoise regions to denote the
hemispheres in **a** and **b**, the 8 octahedral faces **c**,
**d**,**f** and the 12 dodecahedral faces (or base pixels) in **e**.
Coastlines are added for orientation.

## Resolution

All grids use the same resolution parameter `nlat_half`, i.e. the number of rings on
one hemisphere, Equator included. The Gaussian grids (full and reduced) do not have
a ring on the equator, so their total number of rings `nlat` is always even and
twice `nlat_half`. Clenshaw-Curtis grids and the HEALPix grids have a ring on the
equator such their total number of rings is always odd and one less than the
Gaussian grids at the same `nlat_half`. 

!!! info "HEALPix grids do not use Nside as resolution parameter"
    The original formulation for HEALPix grids use ``N_{side}``, the number of grid
    points along the edges of each basepixel (8 in the figure above),
    SpeedyWeather.jl uses `nlat_half`, the number of rings on one hemisphere, Equator included,
    for all grids. This is done for consistency across grids.

## Truncation

A given spectral resolution can be matched to a variety of grid resolutions. A _cubic_ grid, for example,
combines a spectral truncation T with a grid resolution N (=`nlat_half`) such that `T + 1 = N`.
Using T31 and an O32 is therefore often abbreviated as Tco31 meaning that the spherical harmonics are
truncated at ``l_{max}=31`` in combination with `N=32`, i.e. 64 latitude rings in total on an octahedral
Gaussian grid.

Let `J` be the total number of rings. Then we have

- ``T \approx J`` for _linear_ truncation
- ``\frac{3}{2}T \approx J`` for _quadratic_ truncation
- ``2T \approx J`` for _cubic_ truncation

and in general ``\frac{m+1}{2}T \approx J`` for _m-th_ order truncation. So the higher the truncaction
order the more grid points are used in combination with the same spectral
resolution. A higher truncation order therefore makes all grid-point calculations more expensive,
but can represent products of terms on the grid (which will have higher wavenumber components) to
a higher accuracy as more grid points are available within a given wavelength.



## Full Gaussian grid

...

## Full Clenshaw-Curtis grid

...

## Octahedral Gaussian grid

...

## The HEALPix grid

Technically, HEALPix grids are a class of grids that tessalate the sphere into faces that are often
called basepixels. For each member of this class there are ``N_\varphi`` basepixels in zonal direction
and ``N_\theta`` basepixels in meridional direction. For ``N_\varphi = 4`` and ``N_\theta = 3`` we obtain
the classical HEALPix grid with ``N_\varphi N_\theta = 12`` basepixels shown above in [Implemented Grids](@ref).
Each basepixel has a quadratic number of grid points in them. There's an equatorial zone where the number
of zonal grid points is constant (always ``2N``, so 32 at ``N=16``) and there are polar caps above and below
the equatorial zone with the border at  ``\cos(\theta) = 2/3`` (``\theta`` in colatitudes).

Following Górski, 2004[^1], the ``z=cos(\theta)`` colatitude of the ``j``-th ring in the north polar cap,
``j=1,...,N_{side}`` with ``2N_{side} = N`` is 

```math
z = 1 - \frac{j^2}{3N_{side}^2}
```
and on that ring, the longitude ``\phi`` of the ``i``-th point is at
```math
\phi = \frac{\pi}{2j}(i-\tfrac{1}{2})
```
In the equatorial belt ``\theta > \theta^*`` for ``\cos(\theta^*) = 2/3`` this changes to
```math
z = \frac{4}{3} - \frac{2j}{3N_{side}}
```




## OctaHEALPix grid

The original

### References

[^1] Górski, Hivon, Banday, Wandelt, Hansen, Reinecke, Bartelmann, 2004. _HEALPix: A FRAMEWORK FOR HIGH-RESOLUTION DISCRETIZATION AND FAST ANALYSIS OF DATA DISTRIBUTED ON THE SPHERE_, The Astrophysical Journal. doi:[10.1086/427976](https://doi.org/10.1086/427976)

