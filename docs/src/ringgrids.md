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

## Creating data on a RingGrid

Every grid in RingGrids has a `data` field, which is a vector containing the data on the grid.
The grid points are unravelled west to east then north to south, meaning that it starts at
90˚N and 0˚E then walks eastward for 360˚ before jumping on the next latitude ring further south,
this way circling around the sphere till reaching the south pole. This may also be called _ring order_.

Data on a 96x48 `Matrix` which follows this ring order can be put on a `FullGaussianGrid` like so
```julia
julia> map
96×48 Matrix{Float32}:
 -0.219801   -0.519598   1.2066    …   0.81304    -1.16023    1.0353
 -0.55615    -1.05712   -0.227948     -2.06369     1.10353    1.60918
  0.446913   -0.856431   1.58896      -1.0894     -0.894315   0.632353
  0.445915   -0.107201  -0.577785      0.574784   -0.825049   1.29677
  1.194      -0.353374   1.30581       0.465554    0.358457  -0.726567
  1.28693     1.43997    0.691283  …  -0.330544   -0.267588   0.181308
  ⋮                                ⋱   ⋮
 -0.432703    0.17233    0.89222       0.888913    1.32787   -0.248779
 -0.404498    0.127172  -0.64237       0.127979   -1.55253   -2.00749
 -0.857746   -0.433251  -0.468293      1.09825    -0.291169   1.07452
  0.375367   -0.218278   0.492855     -0.287976    0.878996  -1.19745
 -0.0619525  -0.129129  -1.35502   …   0.0824819   0.481736   0.845638

julia> grid = FullGaussianGrid(map)
4608-element, 48-ring FullGaussianGrid{Float32}:
 -0.21980117
 -0.5561496
  0.44691312
  0.4459149
  1.1940043
  1.2869298
  ⋮
 -0.24877885
 -2.007495
  1.0745221
 -1.197454
  0.84563845
```
A full Gaussian grid has always ``2N`` x ``N`` grid points, but a `FullClenshawGrid` has ``2N`` x ``N-1``, if those dimensions don't match, the creation will throw an error. To reobtain the data from a grid, you can access its `data` field
which returns a normal `Vector`
```julia
julia> grid.data
4608-element Vector{Float32}:
 -0.21980117
 -0.5561496
  0.44691312
  0.4459149
  1.1940043
  1.2869298
  ⋮
 -0.24877885
 -2.007495
  1.0745221
 -1.197454
  0.84563845
```
Which can be reshaped to reobtain `map` from above. Alternatively you can `Matrix(grid)` to do this in one step
```julia
julia> map == Matrix(FullGaussianGrid(map))
true
```
You can also use `zeros`,`ones`,`rand`,`randn` to create a grid, whereby `nlat_half`, i.e. the number of latitude
rings on one hemisphere, Equator included, is used as a resolution parameter and here as a second argument.
```julia
julia> nlat_half = 4
julia> grid = randn(OctahedralGaussianGrid{Float16},nlat_half)
208-element, 8-ring OctahedralGaussianGrid{Float16}:
 -1.868
  0.493
  0.3142
  1.871
  1.349
  0.623
  ⋮
  1.064
  0.4346
 -0.641
  0.1445
  0.3909
```
and any element type `T` can be used for `OctahedralGaussianGrid{T}` and similar for other grid types.

## Indexing RingGrids

All RingGrids have a single index `ij` which follows the ring order. While this is obviously not super
exciting here are some examples how to make better use of the information that the data sits on a grid.

We obtain the latitudes of the rings of a grid by calling `get_latd` (`get_lond` is only defined for full
grids, or use `get_latdlonds` for latitudes, longitudes per grid point not per ring)
```julia
julia> latd = get_latd(grid)
8-element Vector{Float64}:
  73.79921362856324
  52.81294318999426
  31.704091745007943
  10.569882312576093
 -10.569882312576093
 -31.704091745007943
 -52.81294318999426
 -73.79921362856324
```
Now we could calculate Coriolis and add it on the grid as follows
```julia
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
```julia
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
```julia
julia> grid
208-element, 8-ring OctahedralGaussianGrid{Float16}:
 -1.868
  0.493
  0.3142
  1.871
  1.349
  0.623
  ⋮
  1.064
  0.4346
 -0.641
  0.1445
  0.3909

julia> interpolate(FullGaussianGrid,grid)
128-element, 8-ring FullGaussianGrid{Float64}:
 -1.8681640625
  0.4482421875
  1.0927734375
  1.4794921875
  0.623046875
 -0.6435546875
  ⋮
 -0.57763671875
  1.064453125
  0.16552734375
 -0.248291015625
  0.329345703125
```
By default this will linearly interpolate (it's an anvil interpolator, see below) onto a grid with the same
`nlat_half`, but we can also coarse-grain or fine-grain by specifying `nlat_half` directly as 2nd argument
```julia
julia> interpolate(FullGaussianGrid,6,grid)
288-element, 12-ring FullGaussianGrid{Float64}:
 -1.248046875
  0.08984375
  0.2763671875
  0.76513671875
  1.1767578125
  ⋮
  0.26416015625
 -0.295166015625
 -0.271728515625
  0.0511474609375
  0.0814208984375
```
So we got from an `8-ring OctahedralGaussianGrid{Float16}` to a `12-ring FullGaussianGrid{Float64}`, so it did
a conversion from `Float16` to `Float64` on the fly too, because the default precision is `Float64` unless
specified. `interpolate(FullGaussianGrid{Float16},6,grid)` would have interpolated onto a grid with element type
`Float16`.

One can also interpolate onto a give cordinate ˚N, ˚E like so
```julia
julia> interpolate(30.0,10.0,grid)
1-element Vector{Float16}:
 0.9297
```
we interpolated the data from `grid` onto 30˚N, 10˚E. To do this simultaneously for many coordinates they can
be packed into a vector too
```julia
julia> interpolate([30.0,40.0,50.0],[10.0,10.0,10.0],grid)
3-element Vector{Float16}:
  0.9297
  0.08887
 -0.929
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
```julia
julia> grid_in = rand(HEALPixGrid,4)
julia> grid_out = zeros(FullClenshawGrid,6)
julia> interp = RingGrids.interpolator(grid_out,grid_in)
```
Now we have created an interpolator `interp` which knows about the geometry where to interpolate
*from* and the coordinates there to interpolate *to*. It is also initialized, meaning it has
precomputed the indices to of `grid_in` that are supposed to be used. It just does not know about
the data of `grid_in` (and neither of `grid_out` which will be overwritten anyway). We can now do
```julia
julia> interpolate!(grid_out,grid_in,interp);
julia> grid_out
264-element, 11-ring FullClenshawGrid{Float64}:
 0.47698810225785837
 0.49923033302273034
 0.5214725637876022
 0.5437147945524742
 ⋮
 0.6277318221906577
 0.5934538182075797
 0.6009226488782581
 0.6083914795489366
```
which is identical to `interpolate(grid_out,grid_in)` but you can reuse `interp` with more data.
The interpolation can also handle various element types (the interpolator `interp` does not have
to be updated for this either)
```julia
julia> grid_out = zeros(FullClenshawGrid{Float16},6);
julia> interpolate!(grid_out,grid_in,interp)
julia> grid_out
264-element, 11-ring FullClenshawGrid{Float16}:
 0.477
 0.4993
 0.5215
 0.544
 0.5493
 0.555
 ⋮
 0.662
 0.628
 0.5933
 0.601
 0.6084
```
and we have converted data from a `HEALPixGrid{Float64}` (which is always default if not specified)
to a `FullClenshawGrid{Float16}` including the type conversion Float64-Float16 on the fly.
Technically there are three data types and their combinations possible: The input data will come
with a type, the output array has an element type and the interpolator has precomputed weights
with a given type. Say we want to go from Float16 data on an `OctahedralGaussianGrid` to Float16
on a `FullClenshawGrid` but using Float32 precision for the interpolation itself, we would do
this by
```julia
julia> grid_in = randn(OctahedralGaussianGrid{Float16},24)
julia> grid_out = zeros(FullClenshawGrid{Float16},24)
julia> interp = RingGrids.interpolator(Float32,grid_out,grid_in)
julia> interpolate!(grid_out,grid_in,interp)
julia> grid_out
4512-element, 47-ring FullClenshawGrid{Float16}:
 -0.954
 -0.724
 -0.4941
 -0.264
 -0.03433
  0.1796
  ⋮
 -0.5703
 -0.3481
 -0.07666
  0.1958
  0.467
```
As a last example we want to illustrate a situation where we would always want to interplate onto
10 coordinates, but their locations may change. In order to avoid recreating an interpolator object
we would do (`AnvilInterpolator` is described in [Anvil interpolator](@ref))
```julia
julia> npoints = 10    # number of coordinates to interpolate onto
julia> interp = AnvilInterpolator(Float32,HEALPixGrid,24,npoints)
```
with the first argument being the number format used during interpolation, then the input grid type,
its resolution in terms of `nlat_half` and then the number of points to interpolate onto. However,
`interp` is not yet initialized as it does not know about the destination coordinates yet. Let's define
them, but note that we already decided there's only 10 of them above.
```julia
julia> latds = collect(0.0:5.0:45.0)
10-element Vector{Float64}:
  0.0
  5.0
  ⋮
 40.0
 45.0

julia> londs = collect(-10.0:2.0:8.0)
10-element Vector{Float64}:
 -10.0
  -8.0
  -6.0
  ⋮
   6.0
   8.0
```
now we can update the locator inside our interpolator as follows
```julia
julia> RingGrids.update_locator!(interp,latds,londs)
```
With data matching the input from above, a `nlat_half=24` HEALPixGrid, and allocate 10-element output vector
```julia
julia> output_vec = zeros(10)
julia> grid_input = rand(HEALPixGrid,24)
```
we can use the interpolator as follows
```julia
julia> interpolate!(output_vec,grid_input,interp)
10-element Vector{Float64}:
 0.3182548251299291
 0.7499448926757676
 0.25733825675836064
  ⋮
 0.2949249541923441
 0.6690698461409016
 0.6159433564856793
```
which is the approximately the same as doing it directly without creating an interpolator first and updating its locator
```julia
julia> interpolate(latds,londs,grid_input)
10-element Vector{Float64}:
 0.31825482404891603
 0.7499448923165136
 0.25733824520344434
  ⋮
 0.294924962125593
 0.6690698486360254
 0.6159433558779497
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