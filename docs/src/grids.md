# Grids

The spectral transform (the [Spherical Harmonic Transform](@ref)) in SpeedyWeather.jl supports any ring-based
equi-longitude grid. Several grids are already implemented but other can be added. The following pages will
describe an overview of these grids and but let's start but how they can be used
```@example grids
using SpeedyWeather
spectral_grid = SpectralGrid(Grid = FullGaussianGrid)
```
The life of every SpeedyWeather.jl simulation starts with a `SpectralGrid` object which defines the
resolution in spectral and in grid-point space. The generator `SpectralGrid()` can take as a keyword
argument `Grid` which can be any of the grids described below. The resolution of the grid, however,
is not directly chosen, but determined from the spectral resolution `trunc` and the `dealiasing`
factor. More in [SpectralGrid](@ref) and [Matching spectral and grid resolution](@ref).

!!! info "RingGrids is a module too!"
    While RingGrids is the underlying module that SpeedyWeather.jl uses for data structs
    on the sphere, the module can also be used independently of SpeedyWeather, for example
    to interpolate between data on different grids. See [RingGrids](@ref)

#### Ring-based equi-longitude grids

SpeedyWeather.jl's spectral transform supports all ring-based equi-longitude grids.
These grids have their grid points located on rings with constant latitude and on these rings the
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

!!! info "Is the FullClenshawGrid a longitude-latitude grid?"
    Short answer: Yes. The `FullClenshawGrid` is *a* specific longitude-latitude grid with equi-angle spacing.
    The most common grids for geoscientific data use regular spacings for 0-360˚E in longitude and
    90˚N-90˚S. The `FullClenshawGrid` does that too, but it does not have a point on the North or South
    pole, and the central latitude ring sits exactly on the Equator. We name it *Clenshaw* following
    the [Clenshaw-Curtis quadrature](https://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature)
    that is used in the Legendre transfrom in the same way as *Gaussian*
    refers to the [Gaussian quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature).

## Implemented grids

All grids in SpeedyWeather.jl are a subtype of `AbstractGrid`, i.e. `<: AbstractGrid`. We further distinguish
between _full_, and _reduced_ grids. Full grids have the same number of longitude points on every latitude
ring (i.e. points converge towards the poles) and reduced grids reduce the number of points towards the poles
to have them more evenly spread out across the globe. More evenly does not necessarily mean that a grid is
equal-area, meaning that every grid cell covers exactly the same area (although the shape changes).

Currently the following full grids `<: AbstractFullGrid` are implemented

- [`FullGaussianGrid`](@ref FullGaussianGrid), a full grid with Gaussian latitudes
- [`FullClenshawGrid`](@ref FullClenshawGrid), a full grid with equi-angle latitudes

and additionally we have `FullHEALPixGrid` and `FullOctaHEALPixGrid` which are the full grid equivalents to the
[HEALPix grid](@ref HEALPixGrid) and the [OctaHEALPix grid](@ref OctaHEALPixGrid) discussed below. Full grid equivalent means that they have
the same latitude rings, but no reduction in the number of points per ring towards the poles and no longitude offset.
Other implemented reduced grids are

- [`OctahedralGaussianGrid`](@ref OctahedralGaussianGrid), a reduced grid with Gaussian latitudes based on an [octahedron](https://en.wikipedia.org/wiki/Octahedron)
- [`OctahedralClenshawGrid`](@ref OctahedralClenshawGrid), similar but based on equi-angle latitudes
- [`HEALPixGrid`](@ref HEALPixGrid), an equal-area grid based on a [dodecahedron](https://en.wikipedia.org/wiki/Rhombic_dodecahedron) with 12 faces
- [`OctaHEALPixGrid`](@ref OctaHEALPixGrid), an equal-area grid from the class of HEALPix grids but based on an octahedron.

An overview of these grids is visualised here, and a more detailed description follows below.

![Overview of implemented grids in SpeedyWeather.jl](https://raw.githubusercontent.com/SpeedyWeather/SpeedyWeather.jl/main/docs/img/grids_comparison.png)

Visualised are each grid's grid points (white dots) and grid faces (white lines).
All grids shown have 16 latitude rings on one hemisphere, Equator included.
The total number of grid points is denoted in the top left of every subplot.
The sphere is shaded with grey, orange and turquoise regions to denote the
hemispheres in **a** and **b**, the 8 octahedral faces **c**,
**d**, **f** and the 12 dodecahedral faces (or base pixels) in **e**.
Coastlines are added for orientation.

## Grid resolution

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
    for all grids. This is done for consistency across grids. We may use ``N_{side}``
    for the documentation or within functions though.

Related: [Effective grid resolution](@ref) and [Available horizontal resolutions](@ref).

## Matching spectral and grid resolution

A given spectral resolution can be matched to a variety of grid resolutions. A _cubic_ grid, for example,
combines a spectral truncation ``T`` with a grid resolution ``N`` (=`nlat_half`) such that ``T + 1 = N``.
Using T31 and an O32 is therefore often abbreviated as Tco31 meaning that the spherical harmonics are
truncated at ``l_{max}=31`` in combination with `N=32`, i.e. 64 latitude rings in total on an octahedral
Gaussian grid. In SpeedyWeather.jl the choice of the order of truncation is controlled with the
`dealiasing` parameter in the [SpectralGrid](@ref) construction.

Let `J` be the total number of rings. Then we have

- ``T \approx J`` for _linear_ truncation, i.e. `dealiasing = 1`
- ``\frac{3}{2}T \approx J`` for _quadratic_ truncation, i.e. `dealiasing = 2`
- ``2T \approx J`` for _cubic_ truncation, , i.e. `dealiasing = 3`

and in general ``\frac{m+1}{2}T \approx J`` for _m-th_ order truncation. So the higher the truncation
order the more grid points are used in combination with the same spectral
resolution. A higher truncation order therefore makes all grid-point calculations more expensive,
but can represent products of terms on the grid (which will have higher wavenumber components) to
a higher accuracy as more grid points are available within a given wavelength. Using a sufficiently
high truncation is therefore one way to avoid aliasing. A quick overview of how the grid resolution
changes when `dealiasing` is passed onto `SpectralGrid` on the `FullGaussianGrid`

| trunc | dealiasing | FullGaussianGrid size |
| ----- | ---------- | --------------------- |
| 31    | 1          | 64x32                 |
| 31    | 2          | 96x48                 |
| 31    | 3          | 128x64                |
| 42    | 1          | 96x48                 |
| 42    | 2          | 128x64                |
| 42    | 3          | 192x96                |
| ...   | ...        | ...                   |

You will obtain this information every time you create a `SpectralGrid(; Grid, trunc, dealiasing)`.

## Interactively exploring the grids

Based on [GeoMakie.jl](https://github.com/MakieOrg/GeoMakie.jl) SpeedyWeather.jl has an extension
(meaning only loaded when also `using GeoMakie`) that defines the `globe` function which
visualises the implemented grids at the desired resolution. With the `CairoMakie` backend
these visualisations are static, but `using GLMakie` they are interactive. For example

```julia
using SpeedyWeather
using GLMakie, GeoMakie

# grid type, resolution parameter nlat_half
globe(FullGaussianGrid, 24)
```

This will open a window with interactive zoom and rotation, visualising a full Gaussian grid
at a resolution of `nlat_half = 24` (i.e. 96x48 grid points). You can visualise all grids
at a wide range of resolutions for non-interactive plotting use `interactive=false`,
which is also what one should do when `using CairoMakie`. Additional keyword arguments are
`coastlines`, `background` among others, check `?globe`. You also can visualise data
on a grid directly this way which will draw polygons for the cell faces, e.g.

```julia
grid = rand(FullGaussianGrid, 24)
globe(grid)
```

## [Full Gaussian grid](@id FullGaussianGrid)

(called `FullGaussianGrid`)

```@example grids
using CairoMakie, GeoMakie    # when using GLMakie, use interactive=true (default) for zoom and rotation
globe(FullGaussianGrid, 24, interactive=false)
save("full_gaussian_grid.png", ans) # hide
nothing # hide
```
![FullGaussianGrid](full_gaussian_grid.png)

The full Gaussian grid is a grid that uses regularly spaced longitudes
with points that do not reduce in number towards the poles. That means for every latitude
``\theta`` the longitudes
``\phi`` are

```math
\phi_i = \frac{2\pi (i-1)}{N_\phi}
```

with ``i = 1, ..., N_\phi`` the in-ring index (1-based, counting from 0˚ eastward) and
``N_\phi`` the number of longitudinal points on the grid.
The first longitude is therefore 0˚, meaning that there is no longitudinal offset on this grid.
There are always twice as many points in zonal direction as there are in meridional,
i.e. ``N_\phi = 2N_\theta``. The latitudes, however, are not regular, but chosen from the
``j``-th zero crossing ``z_j(l)`` of the ``l``-th
[Legendre polynomial](https://en.wikipedia.org/wiki/Legendre_polynomials).
For ``\theta`` in latitudes
```math
\sin(\theta_j) = z_j(l) 
```
As it can be easy to mix up latitudes, colatitudes and as the Legendre polynomials are defined in ``[0, 1]``
an overview of the first Gaussian latitudes (approximated for ``l>2`` for brevity)

| ``l`` | Zero crossings ``z_j``      | Latitudes [˚N] |
| ----- | --------------------------- | -------------- |
| 2     | ``\pm \tfrac{1}{\sqrt{3}}`` | ``\pm 35.3...`` |
| 4     | ``\pm 0.34..., \pm 0.86...``| ``\pm 19.9..., \pm 59.44...`` |
| 6     | ``\pm 0.24..., \pm 0.66..., \pm 0.93...`` | ``\pm 13.8..., \pm 41.4..., \pm 68.8...``  |

Only even Legendre polynomials are used, such that there is always an even number of latitudes,
with no latitude on the Equator. As you can already see from this short table, the Gaussian
latitudes do not nest, i.e. different resolutions through different ``l`` do not share
latitudes. The latitudes are also only approximately evenly spaced. 
Due to the Gaussian latitudes, a spectral transform with a full Gaussian grid is *exact*
as long as the truncation is at least quadratic, see [Matching spectral and grid resolution](@ref).
Exactness here means that only rounding errors occur in the transform, meaning that the
transform error is very small compared to other errors in a simulation.
This property arises from that property of the
[Gauss-Legendre quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature),
which is used in the [Spherical Harmonic Transform](@ref) with a full Gaussian grid.

On the full Gaussian grid there are in total ``N_\phi N_\theta`` grid points, which are squeezed
towards the poles, making the grid area smaller and smaller following a cosine.
But no points are on the poles as ``z=-1`` or ``1`` is never a zero crossing of the Legendre polynomials.

## [Octahedral Gaussian grid](@id OctahedralGaussianGrid)

(called `OctahedralGaussianGrid`)

```@example grids
using CairoMakie, GeoMakie    # when using GLMakie, use interactive=true (default) for zoom and rotation
globe(OctahedralGaussianGrid, 24, interactive=false)
save("octahedral_gaussian_grid.png", ans) # hide
nothing # hide
```
![OctahedralGaussianGrid](octahedral_gaussian_grid.png)

The octahedral Gaussian grid is a reduced grid, i.e. the number of longitudinal points reduces
towards the poles. It still uses the Gaussian latitudes from the [full Gaussian grid](@ref FullGaussianGrid)
so the exactness property of the spherical harmonic transform also holds for this grid.
However, the longitudes ``\phi_i`` with ``i = 1, ..., 16+4j`` on the ``j``-th latitude ring
(starting with 1 around the north pole), ``j=1, ..., \tfrac{N_\theta}{2}``,
are now, on the northern hemisphere,

```math
\phi_i = \frac{2\pi (i-1)}{16 + 4j}.
```

We start with 20 points, evenly spaced, starting at 0˚E, around the first latitude ring below the north pole.
The next ring has 24 points, then 28, and so on till reaching the Equator (which is not a ring).
For the southern hemisphere all points are mirrored around the Equator.
For more details see Malardel, 2016[^M16].

Note that starting with 20 grid points on the first ring is a choice that ECMWF made with their grid
for accuracy reasons. An octahedral Gaussian grid can also be defined starting with fewer grid points
on the first ring. However, in SpeedyWeather.jl we follow ECMWF's definition. 

The grid cells of an octahedral Gaussian grid are not exactly equal area, but are usually within
a factor of two. This largely solves the efficiency problem of having too many grid points near
the poles for computational, memory and data storage reasons.

## [Octaminimal Gaussian grid](@id OctaminimalGaussianGrid)

(called `OctaminimalGaussianGrid`)

```@example grids
using CairoMakie, GeoMakie    # when using GLMakie, use interactive=true (default) for zoom and rotation
globe(OctaminimalGaussianGrid, 24, interactive=false)
save("octaminimal_gaussian_grid.png", ans) # hide
nothing # hide
```
![OctaminimalGaussianGrid](octaminimal_gaussian_grid.png)

The `OctaminimalGaussianGrid` is similar to the `OctahedralGaussianGrid` but starts with
4 points around the poles. It therefore minimizes the number of grid points for the
Gaussian grids at the cost of a (somewhat) inexact spectral transform. But for `nlat_half = 24`
(i.e. 48 latitude rings) this grid reduces the number of horizontal grid points from 3168
to 2400, i.e. -24% which speeds up the computation on the grid (dynamics and physics) as well
as the spectral transform. The octaminimal Gaussian grid is intended to be used for low
resolutions as for higher resolutions the relative increase in the number of grid points
with the octahedral Gaussian grid becomes negligible. 

Furthermore the longitudes are rotated: Instead of no offset where the first point on
every ring starts at 0˚E, an offset of ``360/2n`` degrees is chosen similar to how
the longitudinal points in the HEALPix grids
([HEALPixGrid](@ref HEALPixGrid) and [OctaHEALPixGrid](@ref OctaHEALPixGrid))
are chosen. This allows for a more even distribution of grid points near the poles.

## [Full Clenshaw-Curtis grid](@id FullClenshawGrid)

(called `FullClenshawGrid`)

```@example grids
using CairoMakie, GeoMakie    # when using GLMakie, use interactive=true (default) for zoom and rotation
globe(FullClenshawGrid, 24, interactive=false)
save("full_clenshaw_grid.png", ans) # hide
nothing # hide
```
![FullClenshawGrid](full_clenshaw_grid.png)

The full Clenshaw-Curtis grid is a regular longitude-latitude grid, but a specific one:
The colatitudes ``\theta_j``, and the longitudes ``\phi_i`` are
```math
\theta_j = \frac{j}{N_\theta + 1}\pi, \quad \phi_i = \frac{2\pi (i-1)}{N_\phi}
```
with ``i`` the in-ring zonal index ``i = 1, ..., N_\phi`` and ``j = 1, ... , N_\theta`` the
ring index starting with 1 around the north pole. There is no grid point on the poles,
but in contrast to the Gaussian grids there is a ring on the Equator.
The longitudes are shared with the full Gaussian grid.
Being a full grid, also the full Clenshaw-Curtis grid suffers from too many grid points
around the poles, this is addressed with the [octahedral Clenshaw-Curtis grid](@ref OctahedralClenshawGrid).

The full Clenshaw-Curtis grid gets its name from the
[Clenshaw-Curtis quadrature](https://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature)
that is used in the Legendre transform (see [Spherical Harmonic Transform](@ref)).
This quadrature relies on evenly spaced latitudes, which also means that this grid
nests, see Hotta and Ujiie[^HU18]. More importantly for our application,
the Clenshaw-Curtis grids (including the octahedral described below) allow for an exact
transform with cubic truncation (see [Matching spectral and grid resolution](@ref)).
Recall that the Gaussian latitudes allow for an exact transform with quadratic truncation,
so the Clenshaw-Curtis grids require more grid points for the same spectral resolution
to be exact. But compared to other errors during a simulation this error may be
masked anyway.

## [Octahedral Clenshaw-Curtis grid](@id OctahedralClenshawGrid)

(called `OctahedralClenshawGrid`)

```@example grids
using CairoMakie, GeoMakie    # when using GLMakie, use interactive=true (default) for zoom and rotation
globe(OctahedralClenshawGrid, 24, interactive=false)
save("octahedral_clenshaw_grid.png", ans) # hide
nothing # hide
```
![OctahedralClenshawGrid](octahedral_clenshaw_grid.png)

In the same as we constructed the octahedral Gaussian grid from the full Gaussian grid,
the octahedral Clenshaw-Curtis grid can be constructed from the full Clenshaw-Curtis grid.
It therefore shares the latitudes with the full grid, but the longitudes with the
octahedral Gaussian grid.

```math
\theta_j = \frac{j}{N_\theta + 1}\pi, \quad \phi_i = \frac{2\pi (i-1)}{16 + 4j}.
```

Notation as before, but note that the definition for ``\phi_i`` only holds for the
northern hemisphere, Equator included. The southern hemisphere is mirrored.
The octahedral Clenshaw-Curtis grid inherits the exactness properties from the full Clenshaw-Curtis grid,
but as it is a reduced grid, it is more efficient in terms of computational aspects
and memory than the full grid. Hotta and Ujiie[^HU18] describe this grid in more detail.

## [HEALPix grid](@id HEALPixGrid)

(called `HEALPixGrid`)

```@example grids
using CairoMakie, GeoMakie    # when using GLMakie, use interactive=true (default) for zoom and rotation
globe(HEALPixGrid, 24, interactive=false)
save("healpix_grid.png", ans) # hide
nothing # hide
```
![HEALPixGrid](healpix_grid.png)

Technically, HEALPix grids are a class of grids that tessalate the sphere into faces that are often
called basepixels. For each member of this class there are ``N_\varphi`` basepixels in zonal direction
and ``N_\theta`` basepixels in meridional direction. For ``N_\varphi = 4`` and ``N_\theta = 3`` we obtain
the classical HEALPix grid with ``N_\varphi N_\theta = 12`` basepixels shown above in [Implemented grids](@ref).
Each basepixel has a quadratic number of grid points in them. There's an equatorial zone where the number
of zonal grid points is constant (always ``2N``, so 32 at ``N=16``) and there are polar caps above and below
the equatorial zone with the border at  ``\cos(\theta) = 2/3`` (``\theta`` in colatitudes).

Following Górski, 2004[^G04], the ``z=cos(\theta)`` colatitude of the ``j``-th ring in the north polar cap,
``j=1, ..., N_{side}`` with ``2N_{side} = N`` is 

```math
z = 1 - \frac{j^2}{3N_{side}^2}
```
and on that ring, the longitude ``\phi`` of the ``i``-th point (``i`` is the in-ring-index) is at
```math
\phi = \frac{\pi}{2j}(i-\tfrac{1}{2})
```
The in-ring index ``i`` goes from ``i=1, ..., 4`` for the first (i.e. northern-most) ring,
``i=1, ..., 8`` for the second ring and ``i = 1, ..., 4j`` for the ``j``-th ring in the northern polar cap.

In the north equatorial belt ``j=N_{side}, ..., 2N_{side}`` this changes to
```math
z = \frac{4}{3} - \frac{2j}{3N_{side}}
```
and the longitudes change to (``i`` is always ``i = 1, ..., 4N_{side}`` in the equatorial belt meaning the
number of longitude points is constant here)
```math
\phi = \frac{\pi}{2N_{side}}(i - \frac{s}{2}), \quad s = (j - N_{side} + 1) \mod 2
```
The modulo function comes in as there is an alternating longitudinal offset from the prime meridian (see [Implemented grids](@ref)). For the southern hemisphere the grid point locations can be obtained by mirror symmetry.

### Grid cell boundaries

The cell boundaries are obtained by setting ``i = k + 1/2`` or ``i = k + 1/2 + j`` (half indices) into the
equations above, such that ``z(\phi, k)``, a function for the cosine of colatitude ``z`` of index ``k`` and
the longitude ``\phi`` is obtained. These are then (northern polar cap)

```math
z = 1 - \frac{k^2}{3N_{side}^2}\left(\frac{\pi}{2\phi_t}\right)^2, \quad z = 1 - \frac{k^2}{3N_{side}^2}\left(\frac{\pi}{2\phi_t - \pi}\right)^2
```
with ``\phi_t = \phi \mod \tfrac{\pi}{2}`` and in the equatorial belt
```math
z = \frac{2}{3}-\frac{4k}{3N_{side}} \pm \frac{8\phi}{3\pi}
```

## [OctaHEALPix grid](@id OctaHEALPixGrid)

(called `OctaHEALPixGrid`)

```@example grids
using CairoMakie, GeoMakie    # when using GLMakie, use interactive=true (default) for zoom and rotation
globe(OctaHEALPixGrid, 24, interactive=false)
save("octahealpix_grid.png", ans) # hide
nothing # hide
```
![OctaHEALPixGrid](octahealpix_grid.png)

While the classic HEALPix grid is based on a [dodecahedron](https://en.wikipedia.org/wiki/Rhombic_dodecahedron),
other choices for ``N_\varphi`` and ``N_\theta`` in the class of HEALPix grids will change the number of faces
there are in zonal/meridional direction. With ``N_\varphi = 4`` and ``N_\theta = 1`` we obtain a HEALPix grid that
is based on an octahedron, which has the convenient property that there are twice as many longitude points around
the equator than there are latitude rings between the poles. This is a desirable for truncation as this matches
the distances too, ``2\pi`` around the Equator versus ``\pi`` between the poles. ``N_\varphi = 6, N_\theta = 2``
or ``N_\varphi = 8, N_\theta = 3`` are other possible choices for this, but also more complicated. See 
Górski, 2004[^G04] for further examples and visualizations of these grids.

We call the ``N_\varphi = 4, N_\theta = 1`` HEALPix grid the OctaHEALPix grid, which combines the equal-area
property of the HEALPix grids with the octahedron that's also used in the `OctahedralGaussianGrid` or the
`OctahedralClenshawGrid`. As ``N_\theta = 1`` there is no equatorial belt which simplifies the grid.
The latitude of the ``j``-th isolatitude ring on the `OctaHEALPixGrid` is defined by
```math
z = 1 - \frac{j^2}{N^2},
```
with ``j=1, ..., N``, and similarly for the southern hemisphere by symmetry. On this grid ``N_{side} = N``
where ``N``= `nlat_half`, the number of latitude rings on one hemisphere, Equator included,
because each of the 4 basepixels spans from pole to pole and covers a quarter of the sphere.
The longitudes with in-ring- index ``i = 1, ..., 4j`` are
```math
\phi = \frac{\pi}{2j}(i - \tfrac{1}{2})
```
and again, the southern hemisphere grid points are obtained by symmetry.

### Grid cell boundaries

Similar to the grid cell boundaries for the HEALPix grid, the OctaHEALPix grid's boundaries are
```math
z = 1 - \frac{k^2}{N^2}\left(\frac{\pi}{2\phi_t}\right)^2, \quad z = 1 - \frac{k^2}{N^2}\left(\frac{\pi}{2\phi_t - \pi}\right)^2
```
The ``3N_{side}^2`` in the denominator of the HEALPix grid came simply ``N^2`` for the OctaHEALPix
grid and there's no separate equation for the equatorial belt (which doesn't exist in the OctaHEALPix grid).


### References

[^G04]: Górski, Hivon, Banday, Wandelt, Hansen, Reinecke, Bartelmann, 2004. _HEALPix: A FRAMEWORK FOR HIGH-RESOLUTION DISCRETIZATION AND FAST ANALYSIS OF DATA DISTRIBUTED ON THE SPHERE_, The Astrophysical Journal. doi:[10.1086/427976](https://doi.org/10.1086/427976)

[^M16]: S Malardel, et al., 2016: _A new grid for the IFS_, ECMWF Newsletter 146. [https://www.ecmwf.int/sites/default/files/elibrary/2016/17262-new-grid-ifs.pdf](https://www.ecmwf.int/sites/default/files/elibrary/2016/17262-new-grid-ifs.pdf)

[^HU18]: Daisuke Hotta and Masashi Ujiie, 2018: _A nestable, multigrid-friendly grid on a sphere for global spectralmodels based on Clenshaw–Curtis quadrature_, Quarterly Journal of the Royal Meteorological Society, DOI: [10.1002/qj.3282](https://doi.org/10.1002/qj.3282)