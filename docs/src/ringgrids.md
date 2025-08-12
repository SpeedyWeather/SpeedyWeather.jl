# RingGrids

RingGrids is a submodule that has been developed for SpeedyWeather.jl which is technically
independent (SpeedyWeather.jl however imports it and so does SpeedyTransforms) and can also
be used without running simulations. It is just not put into its own respective repository.

RingGrids defines several iso-latitude grids, which are mathematically described in the
section on [Grids](@ref). In brief, they include the regular latitude-longitude grids
(here called `FullClenshawGrid`) as well as grids which latitudes are shifted to
the Gaussian latitudes and _reduced_ grids, meaning that they have a decreasing number
of longitudinal points towards the poles to be more equal-area than _full_ grids.

# Defined grids

RingGrids defines and exports the following grids, see also [Implemented grids](@ref) for a
more mathematical description.

Full grids
```@example ringgrids
using SpeedyWeather.RingGrids
using InteractiveUtils # hide
subtypes(RingGrids.AbstractFullGrid)
```

and reduced grids
```@example ringgrids
subtypes(RingGrids.AbstractReducedGrid)
```

The following explanation of how to use these can be mostly applied to any of them, however,
there are certain functions that are not defined, e.g. the full grids can be trivially converted
to a `Matrix` (i.e. they are *rectangular* grids) but not the `OctahedralGaussianGrid`.

!!! note "What is a ring?"
    We use the term *ring*, short for *iso-latitude ring*, to refer to a sequence of grid points
    that all share the same latitude. A latitude-longitude grid is a ring grid, as it organises
    its grid-points into rings. However, other grids, like the
    [cubed-sphere](https://en.wikipedia.org/wiki/Quadrilateralized_spherical_cube)
    are not based on iso-latitude rings. SpeedyWeather.jl only works with ring grids
    because its a requirement for the [Spherical Harmonic Transform](@ref) to be efficient.
    See [Grids](@ref).

## Grid versus Field

With "grid" we mean the discretization of space. Also called tesselation given that we are tiling
a space with polygons, we subdivide the sphere into grid cells, with vertices, faces and centres.
In that sense, a grid does not contain any data it purely describes the location of grid cells.
Grids in RingGrids are identified by their name, e.g.
FullGaussianGrid, and a resolution parameter where we use `nlat_half` (the number of latitudes
on one hemisphere, Equator included) for all grids. This is because some grids have an even number
some an odd number number of latitudes so not all `nlat` are valid, but all `nlat_half` are.
While an instance of a grid stores some precomputed arrays to facilitate faster indexing
it does not store coordinates and similar grid information, these can be recomputed on the fly
whenever needed given a grid. All grids are considered to be two-dimensional, so they do not
contain information about the vertical or time, for example. The horizontal grid points are unravelled into
a vector, starting at 0˚E at the north pole, going first east, then ring by ring to the south pole.

Data on a grid is called a `Field`, many variables, like temperature are a field. Surface temperature
would be a 2D field (though represented as a vector), temperature on several vertical layers of
the atmosphere would be 3D (data represented as a matrix, horizontal x vertical), including
time would make it 4D. Several fields can share the same grid. Given that the grid is always
two-dimensional, a 2D and 3D field can also share the same grid, leaving the 3rd dimension
not further specified for flexibility. 

## Creating a grid

All grids are specified by name and the resolution parameter `nlat_half::Integer` (number of latitude rings
on one hemisphere, Equator included). An instance
of a grid is simply created by

```@example ringgrids
grid = FullGaussianGrid(24)
```

As a second argument `architecture` can be provided which helps to share information on the
computing architecture (CPU/GPU) but this will not be further explained here. 

## Accessing coordinates

With a `grid` you can get the coordinates through `get_lat`, `get_latd`, `get_colat`, `get_lond`
but note that the latter is only defined for full grids as the reduced grids do not share the same
longitudes across rings. To obtain the longitudes and latitudes for all grid points use
`get_londlatds`, `get_lonlats`, `get_loncolats`, e.g.

```@example ringgrids
grid = OctaminimalGaussianGrid(2)   # a tiny grid with 4 latitudes only
get_latd(grid)
```

## Creating a Field

Creating `Field`, that means data on a grid can be done in many ways, for example using
`zeros`, `ones, `rand`, or `randn`

```@example ringgrids
grid = HEALPixGrid(2)               # smallest HEALPix
field = zeros(grid)                 # 2D field
field = rand(grid, 10)              # 3D field with 10 vertical layers or time steps
field = randn(Float32, grid, 2, 2)  # 4D using Float32 as element type
```

Note that in this case all fields share the same `grid`.
Or the grid can be created on the fly when the grid type is specified, followed by `nlat_half`.
Every `?Grid` also has a corresponding `?Field` type, e.g.

```@example ringgrids
field = zeros(OctaminimalGaussianGrid, 2)   # nlat_half=2
field = HEALPixField(undef, 2)              # using undef initializor
field = HEALPixField{Float16}(undef, 2, 3)  # using Float16 as eltype
```

## Creating a Field from data

A `field` has `field.data` (some `AbstractArray{T, N}`) and `field.grid` (some `AbstractGrid` as described above).
The first dimension of `data` describes the horizontal as the grid points on every grid (full and reduced)
are unravelled west to east then north to south, meaning that it starts at
90˚N and 0˚E then walks eastward for 360˚ before jumping on the next latitude ring further south,
this way circling around the sphere till reaching the south pole. This may also be called _ring order_.

Data in a `Matrix` which follows this ring order can be put on a `FullGaussianGrid` like so
```@example ringgrids
data = randn(Float32, 8, 4)
field = FullGaussianGrid(data, input_as=Matrix)
field = FullGaussianField(data, input_as=Matrix)    # equivalent
```
Both return a field (there is no data in the grid itself) and create an instance of
`FullGaussianGrid` on the fly. The `input_as=Matrix` is required to denote that the horizontal
dimension is not unravelled into a vector, triggering a `reshape` internally, `input_as=Vector`
is the default. A full Gaussian grid has always ``2N`` x ``N`` grid points, but a `FullClenshawGrid`
has ``2N`` x ``N-1``, if those dimensions don't match, the creation will throw an error.

If you have the data and know which grid it comes one you can also create
a `Field` by providing both
```@example ringgrids
data = randn(Float32, 8, 4)         # data of some shape
grid = OctaminimalGaussianGrid(1)   # you need to know the nlat_half (here 1) of that grid!
field = Field(data, grid)
```
But you can also automatically let `nlat_half` be calculated from the shape of the data.
Note that you have to provide the name of the field though, `FullGaussianField` here to
create a `Field` on a `FullGaussianGrid`.
```@example ringgrids
field = FullGaussianField(data, input_as=Matrix)
```
To return to the original data array you can reshape the data of a full grid (which is representable as a matrix) as follows

```@example ringgrids
data == Matrix(FullGaussianField(data, input_as=Matrix))
```
which in general is `Array(field, as=Vector)` for no reshaping (equivalent to `field.data`
including possible conversion to `Array`) and `Array(field, as=Matrix)` with reshaping
(full grids only).

## Visualising Fields

As only the full fields can be reshaped into a matrix, the underlying data structure of any `AbstractField`
is a vector for consistency. As shown in the examples above, one can therefore inspect the data as if it was a vector.
But a field has through `field.grid` all the geometric information available to plot
it on a map, SpeedyWeather also implements extensions for [Makie's](https://github.com/MakieOrg/Makie.jl)
and [UnicodePlots's](https://github.com/JuliaPlots/UnicodePlots.jl)' `heatmap`, also see
[Visualisation via Makie](@ref) and [Visualisation via UnicodePlots](@ref).

```@example ringgrids
using CairoMakie    # triggers loading of Makie extension, or do using UnicodePlots instead!
grid = OctahedralGaussianGrid(24)
field = randn(grid)
heatmap(field)
```

Reduced fields are automatically interpolated to the corresponding full fields so that they can be visualised
as a matrix.

## Indexing Fields

All fields have a single index `ij` which follows the ring order. While this is obviously not super
exciting here are some examples how to make better use of the information that the data sits on a
two-dimensional grid covering the sphere. We obtain the latitudes of the rings of a grid by calling
`get_latd` (`get_lond` is only defined for full
grids, or use `get_londlatds` for latitudes, longitudes per grid point not per ring).
Now we could calculate Coriolis and add it on the grid as follows

```@example ringgrids
grid = OctahedralClenshawGrid(5)    # define a grid
field = randn(grid)                 # random data on that grid
latd = get_latd(grid)               # get a vector of latitudes [˚N] per ring, north to south

rotation = 7.29e-5                  # angular frequency of Earth's rotation [rad/s]
coriolis = 2rotation*sind.(latd)    # vector of coriolis parameters per latitude ring

for (j, ring) in enumerate(eachring(field))     # loop over every ring j
    f = coriolis[j]
    for ij in ring                              # loop over every longitude point on that ring
        field[ij] += f
    end
end
```

`eachring(field)` accesses `field.grid.rings` a precomputed vector of `UnitRange` indices, such that
we can loop over the ring index `j` (`j=1` being closest to the North pole) pull the coriolis parameter
at that latitude and then loop over all in-ring indices `i` (changing longitudes) to do something on the grid.
Something similar can be done to scale/unscale with the cosine of latitude for example.

We can always loop over all horizontal grid points with `eachgridpoint` and over every other dimensions
with `eachlayer`, e.g. for 2D fields you can do
```@example ringgrids
for ij in eachgridpoint(grid)
    field[ij]
end
```
or use `eachindex` instead. For 3D fields `eachindex` loops over all elements, including the 3rd dimension
but `eachgridpoint` would only loop over the horizontal. To loop with an index `k` over all
additional dimensions (vertical, time, ...) do
```@example ringgrids
field = zeros(grid, 2, 3)           # 4D, 2D defined by grid x 2 x 3
for k in eachlayer(field)           # loop over 2 x 3
    for ij in eachgridpoint(field)  # loop over 2D grid points
        field[ij, k]
    end
end
```

## Interpolation between grids

In most cases we will want to use RingGrids so that our data directly comes with the geometric information
of where the grid-point is one the sphere. We have seen how to use `get_latd`, `get_lond`, ... for that
above. This information generally can also be used to interpolate our data from one grid to another or to
request an interpolated value on some coordinates. Using a `OctahedralGaussianField`
(data on an `OctahedralGaussianGrid`) from above we can use the `interpolate` function to get it onto a
`FullGaussianGrid` (or any other grid for purpose)
```@example ringgrids
grid = OctahedralGaussianGrid(4)
field = randn(Float32, grid)
interpolate(FullGaussianGrid, field)
nothing # hide
```
By default this will linearly interpolate (it's an [Anvil interpolator](@ref), see below) onto a grid with the same
`nlat_half`, but we can also coarse-grain or fine-grain by specifying any other `grid` as an instance not a type
(or provide `nlat_half` as the second argument), e.g.
```@example ringgrids
output_grid = HEALPixGrid(6)
interpolate(output_grid, field)
interpolate(HEALPixGrid, 6, field)  # equivalent
nothing # hide
```
To change the number format during interpolation you can preallocate the output field
```@example ringgrids
output_field = Field(Float16, output_grid)
interpolate!(output_field, field)
nothing # hide
```
Which will convert from `Float64` to `Float16` on the fly too. Note that we use `interpolate!` here not
`interpolate` without `!` as the former writes into `output_field` and therefore changes it in place.

## Interpolation on coordinates

One can also interpolate a 2D field onto a given coordinate ˚N, ˚E like so
```@example ringgrids
interpolate(30.0, 10.0, field)
```
we interpolated the data from `field` onto 30˚N, 10˚E. To do this simultaneously for many coordinates they can
be packed into a vector too
```@example ringgrids
interpolate([30.0, 40.0, 50.0], [10.0, 10.0, 10.0], field)
```
which returns the data on `field` at 30˚N, 40˚N, 50˚N, and 10˚E respectively. Note how the interpolation here
retains the element type of `field`.

## Performance for RingGrid interpolation

Every time an interpolation like `interpolate(30.0, 10.0, field)` is called, several things happen, which
are important to understand to know how to get the fastest interpolation out of this module in a given situation.
Under the hood an interpolation takes three arguments

- output array
- input field
- interpolator

The output array is just an array into which the interpolated data is written, providing this prevents
unnecessary allocation of memory in case the destination array of the interpolation already exists.
Such a destination can be a field but it does not have to be (see [Interpolation on coordinates](@ref)).
The input field contains the data which is subject to interpolation, it must come on a ring grid, however,
its coordinate information is actually already in the interpolator. The interpolator knows about the
geometry of the grid of the field and the coordinates it is supposed to interpolate onto.
It has therefore precalculated the indices that are needed to access the right data on the input field
and the weights it needs to apply in the actual interpolation operation. The only thing it does not
know is the actual data values of that field. So in the case you want to interpolate from field A
to field B many times, you can just reuse the same interpolator. If you want to change the coordinates
of the output but its total number of points remain constant then you can update the locator
inside the interpolator and only else you will need to create a new interpolator. Let's look at this
in practice. Say we have two grids and want to interpolate between them
```@example ringgrids
grid_in = HEALPixGrid(4)
grid_out = FullClenshawGrid(6)
interp = RingGrids.interpolator(grid_out, grid_in)
```
Now we have created an interpolator `interp` which knows about the geometry where to interpolate
*from* and the coordinates there to interpolate *to*. It is also initialized, meaning it has
precomputed the indices to of `grid_in` that are supposed to be used. It just does not know about
the data (it has only seen grids, no fields). We can now do
```@example ringgrids
field_in = rand(grid_in)
field_out = zeros(grid_out)
interpolate!(field_out, field_in, interp)
nothing # hide
```
which is identical to `interpolate(field_out, field_in)` but you can reuse `interp` for other data.
The interpolation can also handle various element types (the interpolator `interp` does not have
to be updated for this either)
```@example ringgrids
field_out = zeros(Float16, grid_out)
interpolate!(field_out, field_in, interp)
nothing # hide
```
and we have converted data from a `HEALPixField{Float64}` (`Float64` is always default if not specified)
to a `FullClenshawField{Float16}` including the type conversion Float64-Float16 on the fly.
Technically there are three data types and their combinations possible: The input data will come
with a type, the output array has an element type and the interpolator has precomputed weights
with a given type. Say we want to go from Float16 data on an `OctahedralGaussianGrid` to Float16
on a `FullClenshawGrid` but using Float32 precision for the interpolation itself, we would do
this by
```@example ringgrids
field_in = randn(OctahedralGaussianField{Float16}, 24)
field_out = zeros(FullClenshawField{Float16}, 24)
interp = RingGrids.interpolator(field_out, field_in, NF=Float32)
interpolate!(field_out, field_in, interp)
nothing # hide
```

As a last example we want to illustrate a situation where we would always want to interpolate onto
10 coordinates, but their locations may change. In order to avoid recreating an interpolator object
we would do (`AnvilInterpolator` is described in [Anvil interpolator](@ref))
```@example ringgrids
npoints = 10    # number of coordinates to interpolate onto
grid = HEALPixGrid(24)
interp = AnvilInterpolator(grid, npoints, NF=Float32)
```
with the first argument being the input grid and then the number of points to interpolate onto.
The number format used for the interpolator is provided as `NF`.
However, `interp` is not yet initialized as it does not know about the destination coordinates yet.
Let's define them, but note that we already decided there's only 10 of them above.
```@example ringgrids
londs = collect(-10.0:2.0:8.0)
latds = collect(0.0:5.0:45.0)
nothing # hide
```
now we can update the locator inside our interpolator as follows
```@example ringgrids
RingGrids.update_locator!(interp, londs, latds)
```
With data matching the input from above, a `nlat_half=24` HEALPixGrid, and allocate 10-element output vector
```@example ringgrids
output_vec = zeros(10)
field_in = rand(grid)
nothing # hide
```
we can use the interpolator as follows
```@example ringgrids
interpolate!(output_vec, field_in, interp)
```
which is the approximately the same as doing it directly without creating an interpolator first and updating its locator
```@example ringgrids
interpolate(londs, latds, field_in)
```
but allows for a reuse of the interpolator. Note that the two output arrays are not exactly identical because we manually
set our interpolator `interp` to use `Float32` for the interpolation whereas the default is `Float64`.

## Anvil interpolator

Currently the only interpolator implemented is a 4-point bilinear interpolator, which schematically works as follows.
Anvil interpolation is the bilinear average of a, b, c, d which are values at grid points in an anvil-shaped configuration
at location x, which is denoted by Δab, Δcd, Δy, the fraction of distances between a-b, c-d, and ab-cd, respectively.
Note that a, c and b, d do not necessarily share the same longitude/x-coordinate.
```
        0..............1    # fraction of distance Δab between a, b
        |<  Δab   >|

0^      a -------- o - b    # anvil-shaped average of a, b, c, d at location x
.Δy                |
.                  |
.v                 x 
.                  |
1         c ------ o ---- d

          |<  Δcd >|
          0...............1 # fraction of distance Δcd between c, d

^ fraction of distance Δy between a-b and c-d.
```
This interpolation is chosen as by definition of the ring grids, a and b share the same latitude, so do c and d,
but the longitudes can be different for all four, a, b, c, d.

## ColumnField 

Additionally to `Field` there is also a `ColumnField` type. `ColumnField` store the data in a column-major format, which is more efficient for column-based computations. As such, when indexing `ColumnField` the first dimension is the vertical dimension, while the second dimension is the horizontal dimension. Otherwise it behaves just like `Field`. To create a `ColumnField` from a `Field` one can use the `transpose` function, which will transpose the data in place and return a `ColumnField`:

```@example ringgrids
grid = OctahedralClenshawGrid(5)    # define a grid
field = randn(grid, 5)  
column_field = transpose(field)
field == transpose(column_field)    # transposing again returns the original Field
```

## Function index

```@autodocs
Modules = [SpeedyWeather.RingGrids]
```
