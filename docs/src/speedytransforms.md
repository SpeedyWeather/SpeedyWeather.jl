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

!!! info "SpeedyTransforms assumes a unit sphere"
    The operators in SpeedyTransforms generally assume a sphere of radius ``R=1``.
    For the transforms themselves that does not make a difference, but the gradient operators
    `div`,`curl`,`∇`,`∇²`,`∇⁻²` omit the radius scaling. You will have to do this manually.
    Also note that the meridional derivate expects a ``\cos^{-1}(\theta)`` scaling.
    More in [Gradient operators](@ref).

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
We can visualize `map` quickly with a UnicodePlot via `plot` (see [Visualising RingGrid data](@ref))
```@example speedytransforms
import SpeedyWeather.RingGrids: plot    # not necessary when `using SpeedyWeather`
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
recomputation or pre-computation of the Legendre polynomials, fore more information see
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
what we started from. The underlying reason is that internally SpeedyWeather uses
`LowerTriangularMatrix`s of size ``l_{max} + 2 \times m_{max} + 1``.
One ``+1`` on both degree and order for 0-based harmonics versus 1-based matrix sizes,
but an additional ``+1`` for the degrees which is required by the meridional derivative.
For consistency, all `LowerTriangularMatrix`s in SpeedyWeather.jl carry this additional degree
but only the vector quantities explicitly make use of it. 
See [Meridional derivative](@ref) for details.

For this interface to SpeedyTransforms this means that on a grid-to-spectral transform you will get
one more degree than orders of the spherical harmonics by default. You can, however, always truncate
this additional degree, say to T5 (hence matrix size is 6x6)
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
latitude-longitude grid, which we call the `FullClenshawGrid`. Note that for the spectral
transform this should not include the poles, so the 96x47 matrix size here corresponds
to 

We now wrap this matrix
therefore to associate it with the necessary grid information
```@example speedytransforms
map = FullClenshawGrid(m)
plot(map)
```
Now we transform into spectral space and call `power_spectrum(::LowerTriangularMatrix)`
```@example speedytransforms
alms = spectral(map)
power = SpeedyTransforms.power_spectrum(alms)
nothing # hide
```

Which returns a vector of power at every wavenumber. By default this is normalized
as average power per degree, you can change that with the keyword argument `normalize=false`.
Plotting this yields
```@example speedytransforms
using UnicodePlots
k = 0:length(power)-1
lineplot(k,power,yscale=:log10,ylim=(1e-15,10),xlim=extrema(k),
    ylabel="power",xlabel="wavenumber",height=10,width=60)
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

## Gradient operators

`SpeedyTransforms` also includes many gradient operators to take derivatives in
spherical harmonics. These are in particular ``\nabla, \nabla \cdot, \nabla \times,
\nabla^2, \nabla^{-2}``. However, the actually implemented operators are,
in contrast to the mathematical [Derivatives in spherical coordinates](@ref)
due to reasons of scaling as follows. Let the implemented operators be
``\hat{\nabla}`` etc.

```math
\hat{\nabla} A = \left(\frac{\partial A}{\partial \lambda}, \cos(\theta)\frac{\partial A}{\partial \theta} \right) =
R\cos(\theta)\nabla A
```
So the zonal derivative omits the radius and the ``\cos^{-1}(\theta)`` scaling.
The meridional derivative adds a ``\cos(\theta)`` due to a recursion relation
being defined that way, which, however, is actually convenient because the whole
operator is therefore scaled by ``R\cos(\theta)``. The curl and divergence operators
expect the input velocity fields to be scaled by ``\cos^{-1}(\theta)``, i.e.

```math
\begin{aligned}
\hat{\nabla} \cdot (\cos^{-1}(\theta)\mathbf{u}) &= \frac{\partial u}{\partial \lambda} +
\cos\theta\frac{\partial v}{\partial \theta} = R\nabla \cdot \mathbf{u}, \\
\hat{\nabla} \times (\cos^{-1}(\theta)\mathbf{u}) &= \frac{\partial v}{\partial \lambda} -
\cos\theta\frac{\partial u}{\partial \theta} = R\nabla \times \mathbf{u}.
\end{aligned}
```

And the Laplace operators omit a ``R^2`` (radius ``R``) scaling, i.e.

```math
\hat{\nabla}^{-2}A = \frac{1}{R^2}\nabla^{-2}A , \quad \hat{\nabla}^{2}A = R^2\nabla^{2}A
```

## Example: Geostrophy

Now, we want to use the following example to illustrate
their use: We have ``u,v`` and want to calculate ``\eta`` in the shallow water
system from it following geostrophy. Analytically we have
```math
-fv = -g\partial_\lambda \eta, \quad fu = -g\partial_\theta \eta
```
which becomes, if you take the divergence of these two equations
```math
\zeta = \frac{g}{f}\nabla^2 \eta
```
Meaning that if we start with ``u,v`` we can obtain the relative vorticity
``\zeta`` and, using Coriolis parameter ``f`` and gravity ``g``, invert
the Laplace operator to obtain displacement ``\eta``. How to do this with
SpeedyTransforms? 

Let us start by generating some data
```@example speedytransforms
using SpeedyWeather

spectral_grid = SpectralGrid(trunc=31,nlev=1)
forcing = SpeedyWeather.JetStreamForcing(spectral_grid)
drag = QuadraticDrag(spectral_grid)
model = ShallowWaterModel(;spectral_grid,forcing,drag)
model.feedback.verbose = false # hide
simulation = initialize!(model);
run!(simulation,period=Day(30))
nothing # hide
```

Now pretend you only have `u,v` to get vorticity (which is actually the prognostic variable in the model,
so calculated anyway...).
```@example speedytransforms
u = simulation.diagnostic_variables.layers[1].grid_variables.u_grid
v = simulation.diagnostic_variables.layers[1].grid_variables.v_grid
vor = SpeedyTransforms.curl(u,v) / spectral_grid.radius
nothing # hide
```
Here, `u,v` are the grid-point velocity fields, and the function `curl` takes in either
`LowerTriangularMatrix`s (no transform needed as all gradient operators act in spectral space),
or, as shown here, arrays of the same grid and size. In this case, the function actually
runs through the following steps
```@example speedytransforms
RingGrids.scale_coslat⁻¹!(u)
RingGrids.scale_coslat⁻¹!(v)

S = SpectralTransform(u,one_more_degree=true)
us = spectral(u,S)
vs = spectral(v,S)

vor = curl(us,vs)
```
(Copies of) the velocity fields are unscaled by the cosine of latitude (see above),
then transformed into spectral space, and the returned `vor` requires a manual division
by the radius. We always unscale vector fields by the cosine of latitude before any
`curl`, or `div` operation, as these omit those.

## One more degree for spectral fields

The `SpectralTransform` in general takes a `one_more_degree` keyword argument,
if otherwise the returned `LowerTriangularMatrix` would be of size 32x32, setting this
to true would return 33x32. The reason is that while most people would expect
square lower triangular matrices for a triangular spectral truncation, all vector quantities
always need one more degree (= one more row) because of a recursion relation in the
meridional gradient. So as we want to take the curl of `us,vs` here, they need this
additional degree, but in the returned lower triangular matrix this row is set to zero.

!!! info "One more degree for vector quantities"
    All gradient operators expect the input lower triangular matrices of shape ``(N+1) \times N``.
    This one more degree of the spherical harmonics is required for the meridional derivative.
    Scalar quantities contain this degree too for size compatibility but they should not
    make use of it. Use `spectral_truncation` to add or remove this degree manually.

## Example: Geostrophy (continued)

Now we transfer `vor` into grid-point space, but specify that we want it on the grid
that we also used in `spectral_grid`. The Coriolis parameter for a grid like `vor_grid`
is obtained, and we do the following for ``f\zeta/g``.

```@example speedytransforms
vor_grid = gridded(vor,Grid=spectral_grid.Grid)
f = SpeedyWeather.coriolis(vor_grid)
fζ_g = spectral_grid.Grid(vor_grid .* f ./ model.planet.gravity)
nothing # hide
```
The last line is a bit awkward for now, as the element-wise multiplication between
two grids escapes the grid and returns a vector that we wrap again into a grid.
We will fix that in future releases. Now we need to apply the inverse
Laplace operator to ``f\zeta/g`` which we do as follows

```@example speedytransforms
fζ_g_spectral = spectral(fζ_g,one_more_degree=true);
η = SpeedyTransforms.∇⁻²(fζ_g_spectral) * spectral_grid.radius^2
η_grid = gridded(η,Grid=spectral_grid.Grid)
nothing # hide
```

Note the manual scaling with the radius ``R^2`` here. We now compare the results
```@example speedytransforms
plot(η_grid)
```
Which is the interface displacement assuming geostrophy. The actual interface
displacement contains also ageostrophy
```@example speedytransforms
plot(simulation.diagnostic_variables.surface.pres_grid)
```
Strikingly similar! The remaining differences are the ageostrophic motions but
also note that the mean is off. This is because geostrophy only use/defines the gradient
of ``\eta`` not the absolute values itself. Our geostrophic ``\eta_g`` has by construction
a mean of zero (that is how we define the inverse Laplace operator) but the actual ``\eta``
is some 1400m higher.

## Functions and type index

```@autodocs
Modules = [SpeedyWeather.SpeedyTransforms]
```