# SpeedyTransforms

SpeedyTransforms is a submodule that has been developed for SpeedyWeather.jl which is technically
independent (SpeedyWeather.jl however imports it) and can also be used without running simulations.
It is just not put into its own respective repository for now.

The SpeedyTransforms are based on [RingGrids](@ref) and
[LowerTriangularMatrices](@ref lowertriangularmatrices) to hold
data in either grid-point space or in spectral space. So you want to read these sections first
for clarifications how to work with these. We will also not discuss mathematical details
of the [Spherical Harmonic Transform](@ref) here, but will focus on the usage of the
SpeedyTransforms module.

The SpeedyTransforms module also implements the gradient operators
``\nabla, \nabla \cdot, \nabla \times, \nabla^2, \nabla^{-2}`` in spectral space.
Combined with the spectral transform, you could for example start with a velocity
field in grid-point space, transform to spectral, compute its divergence and transform
back to obtain the divergence in grid-point space. Examples are outlined in [Gradient operators](@ref).

## Notation: Spectral resolution

There are different ways to describe the spectral resolution, the truncation wavenumber (e.g. T31),
the maximum degree ``l`` and order ``m`` of the spherical harmonics (e.g. ``l_{max}=31``, ``m_{max} = 31``),
or the size of the lower triangular matrix, e.g. 32x32. In this example, they are all equivalent.
We often use the truncation, i.e. T31, for brevity but sometimes it is important to describe
degree and order independently (see for example [One more degree for spectral fields](@ref)).
Note also how truncation, degree and order are 0-based, but matrix sizes are 1-based. 


## Example transform

Lets start with a simple transform. We could be `using SpeedyWeather` but to be more verbose
these are the modules required to load

```@example speedytransforms
using SpeedyWeather.RingGrids
using SpeedyWeather.LowerTriangularMatrices
using SpeedyWeather.SpeedyTransforms
```

As an example, we want to transform the ``l=m=1`` spherical harmonic from spectral space in `alms`
to grid-point space. Note, the ``+1`` on both degree (first index) and order (second index) for
0-based harmonics versus 1-based matrix indexing, see [Size of `LowerTriangularArray`](@ref).
Create a `LowerTriangularMatrix` for T5 resolution, i.e. 6x6 matrix size

```@example speedytransforms
alms = zeros(LowerTriangularMatrix{ComplexF64}, 6, 6)     # spectral coefficients T5
alms[2, 2] = 1                                            # only l=1, m=1 harmonic
alms
```
Now `transform` is the function that takes spectral coefficients `alms` and converts
them a grid-point space `map` (or vice versa)

```@example speedytransforms
map = transform(alms)
```
By default, the `transforms` transforms onto a [`FullGaussianGrid`](@ref FullGaussianGrid) unravelled here
into a vector west to east, starting at the prime meridian, then north to south, see [RingGrids](@ref).
We can visualize `map` quickly with a UnicodePlot via `heatmap` (see [Visualising RingGrid data](@ref)),
or alternatively in higher quality after `using CairoMakie` or `usin GLMakie`, see
[Visualisation via Makie](@ref) too

```@example speedytransforms
using UnicodePlots
heatmap(map)
```
Yay! This is the what the ``l=m=1`` spherical harmonic is supposed to look like!
Now let's go back to spectral space with `transform`
```@example speedytransforms
alms2 = transform(map)
```
Comparing with `alms` from above you can see that the transform is exact up to a typical rounding error
from `Float64`. 
```@example speedytransforms
alms ≈ alms2
```
YAY! The transform is typically idempotent, meaning that either space may hold information
that is not exactly representable in the other but the first two-way transform will remove that
so that subsequent transforms do not change this any further. However, also note here that
the default [`FullGaussianGrid`](@ref FullGaussianGrid) is an exact grid, inexact grids usually have
a transform error that is larger than the rounding error from floating-point arithmetic.

## Transform onto another grid

While the default grid for [SpeedyTransforms](@ref) is the [`FullGaussianGrid`](@ref FullGaussianGrid)
we can transform onto other grids by specifying `Grid` too
```@example speedytransforms
map = transform(alms, Grid=HEALPixGrid)
heatmap(map)
```
which, if transformed back, however, can yield a larger transform error as discussed above
```@example speedytransforms
transform(map)
```
On such a coarse grid the transform error (absolute and relative) is about ``10^{-2}``, this decreases
for higher resolution. The `transform` function will choose a corresponding
grid-spectral resolution (see [Matching spectral and grid resolution](@ref)) following quadratic
truncation, but you can always truncate/interpolate in spectral space with `spectral_truncation`,
`spectral_interpolation` which takes `trunc` = ``l_{max} = m_{max}`` as second argument
```@example speedytransforms
spectral_truncation(alms, 2)
```
Yay, we just chopped off ``l > 2`` from `alms` which contained the harmonics up to degree and
order 5 before.
If the second argument in `spectral_truncation` is larger than `alms` then it will automatically
call `spectral_interpolation` and vice versa. Also see [Interpolation on RingGrids](@ref)
to interpolate directly between grids. If you want to control directly the resolution of the
grid you want to `transform` onto, use the keyword `dealiasing` (default: 2 for quadratic,
see [Matching spectral and grid resolution](@ref)).
But you can also provide a `SpectralTransform` instance to reuse a precomputed spectral transform.
More on that now.

## [The `SpectralTransform` struct](@id SpectralTransform)

The function `transform` only with arguments as shown above,
will create an instance of `SpectralTransform` under the hood.
This object contains all precomputed information that is required for the transform, either way: 
The Legendre polynomials, pre-planned Fourier transforms, precomputed gradient, divergence and
curl operators, the spherical harmonic eigenvalues among others. Maybe the most intuitive way to
create a `SpectralTransform` is to start with a `SpectralGrid`, which already defines
which spectral resolution is supposed to be combined with a given grid.
```@example speedytransforms
using SpeedyWeather
spectral_grid = SpectralGrid(NF=Float32, trunc=5, Grid=OctahedralGaussianGrid, dealiasing=3)
```
(We `using SpeedyWeather` here as `SpectralGrid` is exported therein).
We also specify the number format `Float32` here to be used for the transform although this
is the default anyway. From `spectral_grid` we now construct a `SpectralTransform` as follows
```@example speedytransforms
S = SpectralTransform(spectral_grid)
```
Note that because we chose `dealiasing=3` (cubic truncation) we now match a T5 spectral field
with a 12-ring octahedral Gaussian grid, instead of the 8 rings as above. So going from
`dealiasing=2` (default) to `dealiasing=3` increased our resolution on the grid while the
spectral resolution remains the same.

Passing on `S` the `SpectralTransform` now allows us to transform directly on the grid
defined therein. Note that we recreate `alms` to be of size 7x6 instead of 6x6 for T5
spectral resolution because SpeedyWeather uses internally [One more degree for spectral fields](@ref)
meaning also that's the default when creating a `SpectralTransform` from a `SpectralGrid`.
But results don't change if the last degree (row) contains only zeros.

```@example speedytransforms
alms = zeros(LowerTriangularMatrix{ComplexF64}, 7, 6)     # spectral coefficients
alms[2, 2] = 1                                            # only l=1, m=1 harmonic

map = transform(alms, S)
heatmap(map)
```

Yay, this is again the ``l=m=1`` harmonic, but this time on a slightly higher resolution
`OctahedralGaussianGrid` as specified in the `SpectralTransform` `S`.
Note that also the number format was converted on the fly to Float32 because that is the number
format we specified in `S`! And from grid to spectral
```@example speedytransforms
alms2 = transform(map, S)
```
As you can see the rounding error is now more like ``10^{-8}`` as we are using Float32
(the [OctahedralGaussianGrid](@ref OctahedralGaussianGrid) is another _exact_ grid).
While for this interface to SpeedyTransforms this means that on a grid-to-spectral transform you will
get one more degree than orders of the spherical harmonics by default. You can, however, always truncate
this additional degree, say to T5 (hence matrix size is 6x6)
```@example speedytransforms
spectral_truncation(alms2, 5, 5)
```
`spectral_truncation(alms2, 5)` would have returned the same, a single argument is then assumed
equal for both degrees and orders. Alternatively, you can also pass on the `one_more_degree=false`
argument to the `SpectralTransform` constructor

```@example speedytransforms
S = SpectralTransform(spectral_grid, one_more_degree=false)
```
As you can see the `7x6 LowerTriangularMatrix` in the description above dropped down to
`6x6 LowerTriangularMatrix`, this is the size of the input that is expected (otherwise
you will get a `BoundsError`).

## `SpectralTransform` generators

While you can always create a `SpectralTransform` from a `SpectralGrid` (which defines
both spectral and grid space) there are other constructors/generators available:

```@example speedytransforms
SpectralTransform(alms)
```

Now we have defined the resolution of the spectral space through `alms` but create
a `SpectralTransform` by making assumption about the grid space. E.g. `Grid=FullGaussianGrid`
by default, `dealiasing=2` and `nlat_half` correspondingly. But you can also pass them
on as keyword arguments, for example

```@example speedytransforms
SpectralTransform(alms, Grid=OctahedralClenshawGrid, nlat_half=24)
```

Only note that you don't want to specify both `nlat_half` and `dealiasing` as you would
otherwise overspecify the grid resolution (`dealiasing` will be ignored in that case).
This also works starting from the grid space

```@example speedytransforms
grid = rand(FullClenshawGrid, 12)
SpectralTransform(grid)
```

where you can also provide spectral resolution `trunc` or `dealiasing`. You can also
provide both a grid and a lower triangular matrix to describe both spaces

```@example speedytransforms
SpectralTransform(grid, alms)
```

and you will precompute the transform between them. For more generators see the
docstrings at `?SpectralTransform`.

## Power spectrum

How to take some data and compute a power spectrum with SpeedyTransforms you may ask.
Say you have some global data in a matrix `m` that looks, for example, like
```@example speedytransforms2
using SpeedyWeather.RingGrids # hide
using SpeedyWeather.LowerTriangularMatrices # hide
using SpeedyWeather.SpeedyTransforms # hide
alms = randn(LowerTriangularMatrix{Complex{Float32}}, 32, 32) # hide
spectral_truncation!(alms, 10) # hide
map = transform(alms, Grid=FullClenshawGrid) # hide
m = Matrix(map) # hide
m
```
You hopefully know which grid this data comes on, let us assume it is a regular
latitude-longitude grid, which we call the `FullClenshawGrid` (in analogy to the Gaussian grid based
on the Gaussian quadrature). Note that for the spectral transform this should not include the poles,
so the 96x47 matrix size here corresponds to 23 latitudes north and south of the Equator respectively
plus the equator (=47).

We now wrap this matrix into a `FullClenshawGrid` (`input_as=Matrix` is required because all
grids organise their data as vectors, see [Creating data on a RingGrid](@ref))
therefore to associate it with the necessary grid information like its coordinates
```@example speedytransforms2
map = FullClenshawGrid(m, input_as=Matrix)

using CairoMakie
heatmap(map)
save("random_pattern.png", ans) # hide
nothing # hide
```
![Random pattern](random_pattern.png)

Now we transform into spectral space and call `power_spectrum(::LowerTriangularMatrix)`
```@example speedytransforms2
alms = transform(map)
power = power_spectrum(alms)
nothing # hide
```

Which returns a vector of power at every wavenumber. By default this is normalized
as average power per degree, you can change that with the keyword argument `normalize=false`.
Plotting this yields
```@example speedytransforms2
using UnicodePlots
k = 0:length(power)-1
lineplot(k, power, yscale=:log10, ylim=(1e-15, 10), xlim=extrema(k),
    ylabel="power", xlabel="wavenumber", height=10, width=60)
```

The power spectrum of our data is about 1 up to wavenumber 10 and then close to zero for
higher wavenumbers (which is in fact how we constructed this fake data). Let us
turn this around and use SpeedyTransforms to create random noise in spectral space
to be used in grid-point space!

## Example: Creating k^n-distributed noise

How would we construct random noise in spectral space that follows a certain
power law and transform it back into grid-point space? Define the wavenumber ``k``
for T31, the spectral resolution we are interested in.
(We start from 1 instead of 0 to avoid zero to the power of something negative).
Now create some normally distributed spectral coefficients but scale them down
for higher wavenumbers with ``k^{-2}``

```@example speedytransforms2
k = 1:32
A = randn(Complex{Float32}, 32, 32)
A .*= k.^-2
alms = LowerTriangularArray(A)
```
We first create a Julia `Matrix` so that the matrix-vector broadcasting `.*= k`
is correctly applied across dimensions of `A` and then convert to a
`LowerTriangularMatrix`.

Awesome. For higher degrees and orders the amplitude clearly decreases!
Now to grid-point space and let us visualize the result
```@example speedytransforms2
map = transform(alms)

using CairoMakie
heatmap(map, title="k⁻²-distributed noise")
save("random_noise.png", ans) # hide
nothing # hide
```
![Random noise](random_noise.png)

You can always access the underlying data in `map` via `map.data` in case you
need to get rid of the wrapping into a grid again!

## Precomputed polynomials and allocated memory

!!! info "Reuse `SpectralTransform`s wherever possible"
    Depending on horizontal and vertical resolution of spectral and grid space,
    a `SpectralTransform` can be become very large in memory. Also the recomputation
    of the polynomials and the planning of the FFTs are expensive compared to the
    actual transform itself. Therefore reuse a `SpectralTransform` wherever possible.

The spectral transform uses a Legendre transform in meridional direction. For this the
Legendre polynomials are required, at each latitude ring this is a ``l_{max} \times m_{max}``
lower triangular matrix. Storing precomputed Legendre polynomials therefore quickly
increase in size with resolution. It is therefore advised to reuse a precomputed
`SpectralTransform` object wherever possible. This prevents transforms to allocate
large memory which would otherwise be garbage collected again after the transform.

You get information about the size of that memory (both polynomials and required scratch memory)
in the terminal "show" of a `SpectralTransform` object, e.g. at T127 resolution
with 8 layers these are

```@example speedytransforms2
using SpeedyWeather
spectral_grid = SpectralGrid(trunc=127, nlayers=8)
SpectralTransform(spectral_grid)
```

## Batched Transforms 

SpeedyTransforms also supports batched transforms. With batched input data the `transform`
is performed along the leading dimension, and all further dimensions are interpreted as
batch dimensions. Take for example 

```@example speedytransforms2
alms = randn(LowerTriangularMatrix{Complex{Float32}}, 32, 32, 5) 
grids = transform(alms)
```

In this case we first randomly generated five (32x32) `LowerTriangularMatrix` that hold the
coefficients and then transformed all five matrices batched to the grid space with the 
transform command, yielding 5 `RingGrids` with each 48-rings. 

## Batched power spectra

SpeedyTransforms also supports power spectra calculated over any additional dimensions
to the leading spherical harmonic dimension (it is unravelled as a vector so the first
only, not the first two...). But the power spectrum is always calculated along that
first spherical-harmonic dimension. For example

```@example speedytransforms2 
alms = randn(LowerTriangularMatrix{Complex{Float32}}, 5, 5, 2) 
power_spectrum(alms)
```
returns the power spectrum for `[..., 1]` in the first column and `[..., 2]` in the second.
This avoids to loop over these additional dimensions, but the result would be the same:

```@example speedytransforms2
power_spectrum(alms[:, 1])
```


## Functions and type index

```@autodocs
Modules = [SpeedyWeather.SpeedyTransforms]
```
