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

## Example transform

Lets start with a simple transform. We could be `using SpeedyWeather` but to be more verbose
these are the modules required to load

```@example speedytransforms
using SpeedyWeather.RingGrids
using SpeedyWeather.LowerTriangularMatrices
using SpeedyWeather.SpeedyTransforms
```

As an example, we want to transform the ``l=m=1`` spherical harmonic from spectral space in `alms`
to grid-point space.

```@example speedytransforms
alms = zeros(LowerTriangularMatrix{ComplexF64},6,6)     # spectral coefficients
alms[2,2] = 1                                           # only l=1,m=1 harmonic
alms
```
Now `gridded` is the function that takes spectral coefficients `alms` and converts
them a grid-point space `map`

```@example speedytransforms
map = gridded(alms)
```
By default, the `gridded` transforms onto a [`FullGaussianGrid`](@ref FullGaussianGrid) unravelled here
into a vector west to east, starting at the prime meridian, then north to south, see [RingGrids](@ref).
We can visualize `map` quickly with a unicodeplot via `plot` (see [Visualising RingGrid data](@ref))
```@example speedytransforms
plot(map)
```
Yay! This is the what the ``l=m=1`` spherical harmonic is supposed to look like!
Now let's go back to spectral space with `spectral`
```@example speedytransforms
alms2 = spectral(map)
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
map = gridded(alms,Grid=HEALPixGrid)
plot(map)
```
which, if transformed back, however, can yield a larger transform error as discussed above
```@example speedytransforms
spectral(map)
```
On such a coarse grid the transform error (absolute and relative) is about ``10^{-2}``, this decreases
for higher resolution. The `gridded` and `spectral` functions will choose a corresponding
grid-spectral resolution (see [Matching spectral and grid resolution](@ref)) following quadratic
truncation, but you can always truncate/interpolate in spectral space with `spectral_truncation`,
`spectral_interpolation` which takes `trunc` = ``l_{max} = m_{max}`` as second argument
```@example speedytransforms
spectral_truncation(alms,2)
```
Yay, we just chopped off ``l > 2`` from `alms` which contained the harmonics up to degree and
order 5 before.
If the second argument in `spectral_truncation` is larger than `alms` then it will automatically
call `spectral_interpolation` and vice versa. Also see [Interpolation on RingGrids](@ref)
to interpolate directly between grids. If you want to control directly the resolution of the
grid `gridded` is supposed to transform onto you have to provide a `SpectralTransform` instance.
More on that now.

## The `SpectralTransform` struct

Both `spectral` and `gridded` create an instance of `SpectralTransform` under the hood. This object
contains all precomputed information that is required for the transform, either way: 
The Legendre polynomials, pre-planned Fourier transforms, precomputed gradient, divergence and
curl operators, the spherical harmonic eigenvalues among others. Maybe the most intuitive way to
create a `SpectralTransform` is to start with a `SpectralGrid`, which already defines
which spectral resolution is supposed to be combined with a given grid.
```@example speedytransforms
using SpeedyWeather
spectral_grid = SpectralGrid(Float32,trunc=5,Grid=OctahedralGaussianGrid,dealiasing=3)
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
spectral resolution remains the same. The `SpectralTransform` also has options for the
recomputation or precomputation of the Legendre polynomials, fore more information see
[(P)recompute Legendre polynomials](@ref).

Passing on `S` the `SpectralTransform` now allows us to transform directly on the grid
defined therein.
```@example speedytransforms
map = gridded(alms,S)
plot(map)
```
Yay, this is again the ``l=m=1`` harmonic, but this time on a slightly higher resolution
`OctahedralGaussianGrid` as specified in the `SpectralTransform` `S`.
Note that also the number format was converted on the fly to Float32 because that is the number
format we specified in `S`! And from grid to spectral
```@example speedytransforms
alms2 = spectral(map,S)
```
As you can see the rounding error is now more like ``10^{-8}`` as we are using Float32
(the [OctahedralGaussianGrid](@ref OctahedralGaussianGrid) is another _exact_ grid).
Note, however, that the returned `LowerTriangularMatrix` is of size 7x6, not 6x6
what we started from. We may change this in the future, but the underlying reason is
that internally SpeedyWeather uses `LowerTriangularMatrix`s of size ``l_{max} + 2 \times m_{max} + 1``.
One ``+1`` on both degree and order for 0-based harmonics versus 1-based matrix sizes,
but an additional ``+1`` for the degrees which is required by the meridional derivative.
For consistency, all `LowerTriangularMatrix`s in SpeedyWeather.jl carry this additional degree
but only the vector quantities explicitly make use of it.
See [Meridional derivative](@ref) for details.

You can, however, always truncate this additional degree, say to T5 (hence matrix size is 6x6)
```@example speedytransforms
spectral_truncation(alms2,5,5)
```
`spectral_truncation(alms2,5)` would have returned the same, a single argument is then assumed
equal for both degrees and orders. Alternatively, you can also pass on the `one_more_degree=false`
argument to the `SpectralTransform` constructor

```@example speedytransforms
S = SpectralTransform(spectral_grid,one_more_degree=false)
```
As you can see the `7x6 LowerTriangularMatrix` in the description above dropped down to
`6x6 LowerTriangularMatrix`, this is the size of the input that is expected (otherwise
you will get a `BoundsError`).

## Power spectrum

How to take some data and compute a power spectrum with SpeedyTransforms you may ask.
Say you have some global data in a matrix `m` that looks, for example, like
```@example speedytransforms
alms = randn(LowerTriangularMatrix{Complex{Float32}},32,32) # hide
spectral_truncation!(alms,10) # hide
map = gridded(alms,Grid=FullClenshawGrid) # hide
m = Matrix(map) # hide
m
```
You hopefully know which grid this data comes on, let us assume it is a regular
latitde-longitude grid, which we call the `FullClenshawGrid`. We now wrap this matrix
therefore to associate it with the necessary grid information
```@example speedytransforms
map = FullClenshawGrid(m)
plot(map)
```
Now we transform into spectral space and call `power_spectrum(::LowerTriangularMatrix)`
```@example speedytransforms
alms = spectral(map)
p = SpeedyTransforms.power_spectrum(alms)
nothing # hide
```

Which returns a vector of power at every wavenumber. By default this is normalized
as average power per degree, you can change that with the keyword argument `normalize=false`.
Plotting this yields
```@example speedytransforms
import PyPlot: PyPlot, savefig, tight_layout # hide
PyPlot.ioff() # hide
import PyPlot: semilogy, xlabel, ylabel
semilogy(p)
xlabel("wavenumber")
ylabel("power")
tight_layout() # hide
savefig("spectrum.svg"); close() # hide
```
![](spectrum.svg)

The power spectrum of our data is about 1 up to wavenumber 10 and then close to zero for
higher wavenumbers (which is in fact how we constructed this fake data). Let us
turn this around and use SpeedyTransforms to create random noise in spectral space
to be used in grid-point space!

## Example: Creating k^n-distributed noise


How would be construct random noise in spectral space that follows a certain
power law and transform it back into grid-point space? Define the wavenumber ``k``
for T31, the spectral resolution we are interested in.
(We start from 1 instead of 0 to avoid zero to the power of something negative).
Now create some normally distributed spectral coefficients but scale them down
for higher wavenumbers with ``k^{-2}``

```@example speedytransforms
k = 1:32
alms = randn(LowerTriangularMatrix{Complex{Float32}},32,32)
alms .*= k.^-2
```
Awesome. For higher degrees and orders the amplitude clearly decreases! Now
to grid-point space and let us visualize the result
```@example speedytransforms
map = gridded(alms)
plot(map)
```
![](power_noise.svg)

You can always access the underlying data in `map` via `map.data` in case you
need to get rid of the wrapping into a grid again!

## (P)recompute Legendre polynomials

The spectral transform uses a Legendre transform in meridional direction. For this the
Legendre polynomials are required, at each latitude ring this is a ``l_{max} \times m_{max}``
lower triangular matrix. Storing precomputed Legendre polynomials therefore quickly
increase in size with resolution. One can recompute them to save memory, but that 
uses more arithmetic operations. There is therefore a memory-compute tradeoff.

For a single transform, there is no need to precompute the polynomials as the
`SpectralTransform` object will be garbage collected again anyway. For low resolution
simulations with many repeated small transforms it makes sense to precompute the
polynomials and SpeedyWeather.jl does that automatically anyway. At very high resolution
the polynomials may, however, become prohibitively large. An example at T127,
about 100km resolution

```@example speedytransforms
spectral_grid = SpectralGrid(trunc=127)
SpectralTransform(spectral_grid,recompute_legendre=false)
```
the polynomials are about 3MB in size. Easy that is not much. But at T1023 on the
O1536 octahedral Gaussian grid, this is already 1.5GB, cubically increasing with the
spectral truncation T. `recompute_legendre=true` (default `false` when
constructing a `SpectralTransform` object which may be reused) would lower this
to kilobytes
```@example speedytransforms
SpectralTransform(spectral_grid,recompute_legendre=true)
```


## Functions and type index

```@autodocs
Modules = [SpeedyWeather.SpeedyTransforms]
```