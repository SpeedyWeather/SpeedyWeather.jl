# RingGrids

RingGrids is a submodule that has been developed for SpeedyWeather.jl which is technically
independent (SpeedyWeather.jl however imports it and so does SpeedyTransforms) and can also
be used without running simulations. It is just not put into its own respective repository.

RingGrids defines several iso-latitude grids, which are mathematically described in the
section on [Grids](@ref). In brief, they include the regular latitude-longitude grids
(here called `FullClenshawGrid`) as well as grids which latitudes are shifted to
the Gaussian latitudes and _reduced_ grids, meaning that they have a decreasing number
of longitudinal points towards the poles to be more equal-area than _full_ grids.

RingGrids defines and exports the following grids:

- full grids: `FullClenshawGrid`, `FullGaussianGrid`, `FullHEALPix`, and `FullOctaHEALPix`
- reduced grids: `OctahedralGaussianGrid`, `OctahedralClenshawGrid`, `OctaHEALPixGrid` and `HEALPixGrid`

The following explanation of how to use these can be mostly applied to any of them, however,
there are certain functions that are not defined, e.g. the full grids can be trivially converted
to a `Matrix` but not the `OctahedralGaussianGrid`.

!!! note "What is a ring?"
    We use the term *ring*, short for *iso-latitude ring*, to refer to a sequence of grid points
    that all share the same latitude. A latitude-longitude grid is a ring grid, as it organises
    its grid-points into rings. However, other grids, like the
    [cubed-sphere](https://en.wikipedia.org/wiki/Quadrilateralized_spherical_cube)
    are not based on iso-latitude rings. SpeedyWeather.jl only works with ring grids
    because its a requirement for the [Spherical Harmonic Transform](@ref) to be efficient.
    See [Grids](@ref).

## Creating data on a RingGrid

Every `grid` in RingGrids has a `grid.data` field, which is a vector containing the data on the grid.
The grid points are unravelled west to east then north to south, meaning that it starts at
90˚N and 0˚E then walks eastward for 360˚ before jumping on the next latitude ring further south,
this way circling around the sphere till reaching the south pole. This may also be called _ring order_.

Data in a `Matrix` which follows this ring order can be put on a `FullGaussianGrid` like so
```@example ringgrids
using SpeedyWeather.RingGrids
map = randn(Float32,8,4)
```

```@example ringgrids
grid = FullGaussianGrid(map)
```
A full Gaussian grid has always ``2N`` x ``N`` grid points, but a `FullClenshawGrid` has ``2N`` x ``N-1``,
if those dimensions don't match, the creation will throw an error. To reobtain the data from a grid,
you can access its `data` field which returns a normal `Vector`
```@example ringgrids
grid.data
```
Which can be reshaped to reobtain `map` from above. Alternatively you can `Matrix(grid)` to do this in one step
```@example ringgrids
map == Matrix(FullGaussianGrid(map))
```
You can also use `zeros`,`ones`,`rand`,`randn` to create a grid, whereby `nlat_half`, i.e. the number of latitude
rings on one hemisphere, Equator included, is used as a resolution parameter and here as a second argument.
```@example ringgrids
nlat_half = 4
grid = randn(OctahedralGaussianGrid{Float16},nlat_half)
```
and any element type `T` can be used for `OctahedralGaussianGrid{T}` and similar for other grid types.

## Visualising RingGrid data

As only the full grids can be reshaped into a matrix, the underyling data structure of any `AbstractGrid`
is a vector. As shown in the examples above, one can therefore inspect the data as if it was a vector.
But as that data has, through its `<:AbstractGrid` type, all the geometric information available to plot
it on a map, RingGrids also exports `plot` function,
based on [UnicodePlots](https://github.com/JuliaPlots/UnicodePlots.jl)' `heatmap`.
```@example ringgrids
nlat_half = 24
grid = randn(OctahedralGaussianGrid,nlat_half)
plot(grid)
```

## Indexing RingGrids

All RingGrids have a single index `ij` which follows the ring order. While this is obviously not super
exciting here are some examples how to make better use of the information that the data sits on a grid.

We obtain the latitudes of the rings of a grid by calling `get_latd` (`get_lond` is only defined for full
grids, or use `get_latdlonds` for latitudes, longitudes per grid point not per ring)
```@example ringgrids
grid = randn(OctahedralClenshawGrid,5)
latd = get_latd(grid)
```
Now we could calculate Coriolis and add it on the grid as follows
```@example ringgrids
rotation = 7.29e-5                  # angular frequency of Earth's rotation [rad/s]
coriolis = 2rotation*sind.(latd)    # vector of coriolis parameters per latitude ring

rings = eachring(grid)
for (j,ring) in enumerate(rings)
    f = coriolis[j]
    for ij in ring
        grid[ij] += f
    end
end
```
`eachring` creates a vector of `UnitRange` indices, such that we can loop over the ring index `j`
(`j=1` being closest to the North pole) pull the coriolis parameter at that latitude and then
loop over all in-ring indices `i` (changing longitudes) to do something on the grid.
Something similar can be done to scale/unscale with the cosine of latitude for example.
We can always loop over all grid-points like so
```@example ringgrids
for ij in eachgridpoint(grid)
    grid[ij]
end
```
or use `eachindex` instead.

## Interpolation on RingGrids

In most cases we will want to use RingGrids so that our data directly comes with the geometric information
of where the grid-point is one the sphere. We have seen how to use `get_latd`, `get_lond`, ... for that
above. This information generally can also be used to interpolate our data from grid to another or to
request an interpolated value on some coordinates. Using our data on `grid` which is an `OctahedralGaussianGrid`
from above we can use the `interpolate` function to get it onto a `FullGaussianGrid` (or any other grid for
purpose)
```@example ringgrids
grid = randn(OctahedralGaussianGrid{Float32},4)
```
```@example ringgrids
interpolate(FullGaussianGrid,grid)
```
By default this will linearly interpolate (it's an [Anvil interpolator](@ref), see below) onto a grid with the same
`nlat_half`, but we can also coarse-grain or fine-grain by specifying `nlat_half` directly as 2nd argument
```@example ringgrids
interpolate(FullGaussianGrid,6,grid)
```
So we got from an `8-ring OctahedralGaussianGrid{Float16}` to a `12-ring FullGaussianGrid{Float64}`, so it did
a conversion from `Float16` to `Float64` on the fly too, because the default precision is `Float64` unless
specified. `interpolate(FullGaussianGrid{Float16},6,grid)` would have interpolated onto a grid with element type
`Float16`.

One can also interpolate onto a give cordinate ˚N, ˚E like so
```@example ringgrids
interpolate(30.0,10.0,grid)
```
we interpolated the data from `grid` onto 30˚N, 10˚E. To do this simultaneously for many coordinates they can
be packed into a vector too
```@example ringgrids
interpolate([30.0,40.0,50.0],[10.0,10.0,10.0],grid)
```
which returns the data on `grid` at 30˚N, 40˚N, 50˚N, and 10˚E respectively. Note how the interpolation here
retains the element type of `grid`.

## Performance for RingGrid interpolation

Every time an interpolation like `interpolate(30.0,10.0,grid)` is called, several things happen, which
are important to understand to know how to get the fastest interpolation out of this module in a given situation.
Under the hood an interpolation takes three arguments

- output vector
- input grid
- interpolator

The output vector is just an array into which the interpolated data is written, providing this prevents
unnecessary allocation of memory in case the destination array of the interpolation already exists.
The input grid contains the data which is subject to interpolation, it must come on a ring grid, however,
its coordinate information is actually already in the interpolator. The interpolator knows about the
geometry of the grid the data is coming on and the coordinates it is supposed to interpolate onto.
It has therefore precalculated the indices that are needed to access the right data on the input grid
and the weights it needs to apply in the actual interplation operation. The only thing it does not
know is the actual data values of that grid. So in the case you want to interpolate from grid A
to grid B many times, you can just reuse the same interpolator. If you want to change the coordinates
of the output grid but its total number of points remain constants then you can update the locator
inside the interpolator and only else you will need to create a new interpolator. Let's look at this
in practice. Say we have two grids an want to interpolate between them
```@example ringgrids
grid_in = rand(HEALPixGrid,4)
grid_out = zeros(FullClenshawGrid,6)
interp = RingGrids.interpolator(grid_out,grid_in)
```
Now we have created an interpolator `interp` which knows about the geometry where to interpolate
*from* and the coordinates there to interpolate *to*. It is also initialized, meaning it has
precomputed the indices to of `grid_in` that are supposed to be used. It just does not know about
the data of `grid_in` (and neither of `grid_out` which will be overwritten anyway). We can now do
```@example ringgrids
interpolate!(grid_out,grid_in,interp)
grid_out
```
which is identical to `interpolate(grid_out,grid_in)` but you can reuse `interp` for other data.
The interpolation can also handle various element types (the interpolator `interp` does not have
to be updated for this either)
```@example ringgrids
grid_out = zeros(FullClenshawGrid{Float16},6);
interpolate!(grid_out,grid_in,interp)
grid_out
```
and we have converted data from a `HEALPixGrid{Float64}` (`Float64` is always default if not specified)
to a `FullClenshawGrid{Float16}` including the type conversion Float64-Float16 on the fly.
Technically there are three data types and their combinations possible: The input data will come
with a type, the output array has an element type and the interpolator has precomputed weights
with a given type. Say we want to go from Float16 data on an `OctahedralGaussianGrid` to Float16
on a `FullClenshawGrid` but using Float32 precision for the interpolation itself, we would do
this by
```@example ringgrids
grid_in = randn(OctahedralGaussianGrid{Float16},24)
grid_out = zeros(FullClenshawGrid{Float16},24)
interp = RingGrids.interpolator(Float32,grid_out,grid_in)
interpolate!(grid_out,grid_in,interp)
grid_out
```
As a last example we want to illustrate a situation where we would always want to interplate onto
10 coordinates, but their locations may change. In order to avoid recreating an interpolator object
we would do (`AnvilInterpolator` is described in [Anvil interpolator](@ref))
```@example ringgrids
npoints = 10    # number of coordinates to interpolate onto
interp = AnvilInterpolator(Float32,HEALPixGrid,24,npoints)
```
with the first argument being the number format used during interpolation, then the input grid type,
its resolution in terms of `nlat_half` and then the number of points to interpolate onto. However,
`interp` is not yet initialized as it does not know about the destination coordinates yet. Let's define
them, but note that we already decided there's only 10 of them above.
```@example ringgrids
latds = collect(0.0:5.0:45.0)
londs = collect(-10.0:2.0:8.0)
nothing # hide
```
now we can update the locator inside our interpolator as follows
```@example ringgrids
RingGrids.update_locator!(interp,latds,londs)
```
With data matching the input from above, a `nlat_half=24` HEALPixGrid, and allocate 10-element output vector
```@example ringgrids
output_vec = zeros(10)
grid_input = rand(HEALPixGrid,24)
nothing # hide
```
we can use the interpolator as follows
```@example ringgrids
interpolate!(output_vec,grid_input,interp)
```
which is the approximately the same as doing it directly without creating an interpolator first and updating its locator
```@example ringgrids
interpolate(latds,londs,grid_input)
```
but allows for a reuse of the interpolator. Note that the two output arrays are not exactly identical because we manually
set our interpolator `interp` to use `Float32` for the interplation whereas the default is `Float64`.

## Anvil interpolator

Currently the only interpolator implemented is a 4-point bilinear interpolator, which schematically works as follows.
Anvil interpolation is the bilinear average of a,b,c,d which are values at grid points in an anvil-shaped configuration
at location x, which is denoted by Δab,Δcd,Δy, the fraction of distances between a-b,c-d, and ab-cd, respectively.
Note that a,c and b,d do not necessarily share the same longitude/x-coordinate.
```
        0..............1    # fraction of distance Δab between a,b
        |<  Δab   >|

0^      a -------- o - b    # anvil-shaped average of a,b,c,d at location x
.Δy                |
.                  |
.v                 x 
.                  |
1         c ------ o ---- d

          |<  Δcd >|
          0...............1 # fraction of distance Δcd between c,d

^ fraction of distance Δy between a-b and c-d.
```
This interpolation is chosen as by definiton of the ring grids, a and b share the same latitude, so do c and d,
but the longitudes can be different for all four, a,b,c,d.

## Function index

```@docs
RingGrids.each_index_in_ring
RingGrids.eachgridpoint
RingGrids.eachring
RingGrids.whichring
RingGrids.get_nlons
```