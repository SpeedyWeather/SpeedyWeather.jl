# Spherical Harmonic Transform

The following sections outline the implementation of the spherical harmonic transform (in short _spectral_ transform)
between the coefficients of the spherical harmonics (the _spectral_ space) and the grid space which can be any of the
[Implemented grids](@ref) as defined by [RingGrids](@ref). This includes the classical full
[Gaussian grid](https://confluence.ecmwf.int/display/FCST/Gaussian+grids), a regular longitude-latitude grid called
the full Clenshaw grid ([FullClenshawGrid](@ref FullClenshawGrid)), ECMWF's octahedral Gaussian grid[^Malardel2016],
and HEALPix grids[^Gorski2004].
SpeedyWeather.jl's spectral transform module SpeedyTransforms is grid-flexible and can be used with any of these,
see [Grids](@ref).

!!! info "SpeedyTransforms is a module too!"
    SpeedyTransform is the underlying module that SpeedyWeather imports to transform between spectral and grid-point
    space, which also implements [Derivatives in spherical coordinates](@ref). You can use this module independently of
    SpeedyWeather for spectral transforms, see [SpeedyTransforms](@ref).

## Inspiration

The spectral transform implemented by SpeedyWeather.jl follows largely Justin Willmert's
[CMB.jl](https://github.com/jmert/CMB.jl) and
[SphericalHarmonicTransforms.jl](https://github.com/jmert/SphericalHarmonicTransforms.jl) package and makes use of
[AssociatedLegendrePolynomials.jl](https://github.com/jmert/AssociatedLegendrePolynomials.jl) and
[FFTW.jl](https://github.com/JuliaMath/FFTW.jl) for the Fourier transform.
Justin described his work in a Blog series [^Willmert2020].

## Spherical harmonics

The [spherical harmonics](https://en.wikipedia.org/wiki/Spherical_harmonics) ``Y_{lm}`` of degree ``l`` and order ``m``
over the longitude ``\phi = (0,2\pi)`` and colatitudes ``\theta = (-\pi/2,\pi/2)``, are

```math
Y_{lm}(\phi, \theta) = \lambda_l^m(\sin\theta) e^{im\phi}
```

with ``\lambda_l^m`` being the pre-normalized associated Legendre polynomials, and ``e^{im\phi}`` are the
complex exponentials (the Fourier modes). Together they form a set of orthogonal basis functions on the sphere.
For an interactive visualisation of the spherical harmonics, see
[here](https://justinwillmert.com/posts/2020/plots-of-the-spherical-harmonics-eigenmodes/).

!!! info "Latitudes versus colatitudes"
    The implementation of the spectral transforms in SpeedyWeather.jl uses colatitudes ``\theta = (0,\pi)``
    (0 at the north pole) but the dynamical core uses latitudes ``\theta = (-\pi/2,\pi/2)`` (``\pi/2`` at the north pole).
    Note: We may also use latitudes in the spherical harmonic transform in the future for consistency. 

## [Synthesis (spectral to grid)](@id synthesis)

The synthesis (or inverse transform) takes the spectral coefficients ``a_{lm}`` and transforms them to grid-point values
``f(\phi,\theta)`` (for the sake of simplicity first regarded as continuous). The synthesis is a linear combination of
the spherical harmonics ``Y_{lm}`` with non-zero coefficients.

```math
f(\phi,\theta) = \sum_{l=0}^{\infty} \sum_{m=-l}^l a_{lm} Y_{lm}(\phi,\theta)
```

We obtain an approximation with a finite set of ``a_{l,m}`` by truncating the series in both degree ``l``
and order ``m`` somehow. Most commonly, a triangular truncation is applied, such that all degrees
after ``l = l_{max}`` are discarded. Triangular because the retained array of the coefficients ``a_{l,m}``
looks like a triangle. Other truncations like rhomboidal have been studied[^Daley78] but are rarely used
since. Choosing ``l_{max}`` also constrains ``m_{max}`` and determines the (horizontal) spectral resolution.
In SpeedyWeather.jl this resolution as chosen as `trunc` when creating the [SpectralGrid](@ref).

For ``f`` being a real-valued there is a symmetry
```math
a_{l,-m} = (-1)^m a^*_{l,+m},
```
meaning that the coefficients at ``-m`` and ``m`` are the same, but the sign of the real and imaginary component
can be flipped, as denoted with the ``(-1)^m`` and the complex conjugate ``a_{l,m}^*``. As we are only dealing with
real-valued fields anyway, we therefore never have to store the negative orders ``-m`` and end up with a lower
triangular matrix of size ``(l_{max}+1) \times (m_{max}+1)`` or technically ``(T+1)^2`` where ``T`` is
the truncation `trunc`. One is added here because the degree ``l`` and order ``m`` use 0-based indexing
but sizes (and so is Julia's indexing) are 1-based.

For correctness we want to mention here that vector quantities require one more degree ``l`` due to the recurrence
relation in the [Meridional derivative](@ref). Hence for practical reasons *all* spectral fields are represented
as a lower triangular matrix of size ``(m_{max} + 2) \times (m_{max} +1)``. And the scalar quantities would just
not make use of that last degree, and its entries would be simply zero. We will, however, for the following
sections ignore this and only discuss it again in [Meridional derivative](@ref).

Another consequence of the symmetry mentioned above is that the zonal harmonics, meaning ``a_{l,m=0}`` have
no imaginary component. Because these harmonics are zonally constant, a non-zero imaginary component would
rotate them around the Earth's axis, which, well, doesn't actually change a real-valued field. 

Following the notation of [^Willmert2020] we can therefore write the truncated synthesis as
```math
f(\phi,\theta) = \sum_{l=0}^{l_{max}} \sum_{m=0}^l (2-\delta_{m0}) a_{lm} Y_{lm}(\phi,\theta)
```
The ``(2-\delta_{m0})`` factor using the Kronecker ``\delta`` is used here because of the symmetry we have
to count both the ``m,-m`` order pairs (hence the ``2``) except for the zonal harmonics which do not have
a pair.

Another symmetry arises from the fact that the spherical harmonics are either symmetric or anti-symmetric
around the Equator. There is an even/odd combination of degrees and orders so that the sign flips like a
checkerboard
```math
Y_{l,m}(\phi,\pi-\theta) = (-1)^{l+m}Y_{lm}(\phi,\phi)
```
This means that one only has to compute the Legendre polynomials for one hemisphere and the other one
follows with this equality.

## [Analysis (grid to spectral)](@id analysis)

Starting in grid-point space we can transform a field ``f(\lambda,\theta)`` into the spectral space of the spherical harmonics by

```math
a_{l,m} = \int_0^{2\pi} \int_{0}^\pi f(\phi,\theta) Y_{l,m}(\phi,\theta) \sin \theta d\theta d\phi
```
Note that this notation again uses colatitudes ``\theta``, for latitudes the ``\sin\theta`` becomes a ``\cos\theta`` and the
bounds have to be changed accordingly to ``(-\frac{\pi}{2},\frac{\pi}{2})``. A discretization with ``N``
grid points at location ``(\phi_i,\theta_i)``, indexed by ``i`` can be written as [^Willmert2020]
```math
\hat{a}_{l,m} = \sum_i f(\phi_i,\theta_i) Y_{l,m}(\phi_i,\theta_i) \sin \theta_i \Delta\theta \Delta\phi
```
The hat on ``a`` just means that it is an approximation, or an estimate of the true ``a_{lm} \approx \hat{a}_{lm}``.
We can essentially make use of the same symmetries as already discussed in [Synthesis](@ref synthesis).
Splitting into the Fourier modes ``e^{im\phi}`` and the Legendre polynomials ``\lambda_l^m(\cos\theta)``
(which are defined over ``[-1,1]`` so the ``\cos\theta`` argument maps them to colatitudes) we have
```math
\hat{a}_{l,m} = \sum_j \left[ \sum_i f(\phi_i,\theta_j) e^{-im\phi_i} \right] \lambda_{l,m}(\theta_j) \sin \theta_j \Delta\theta \Delta\phi
```
So the term in brackets can be separated out as long as the latitude ``\theta_j`` is constant,
which motivates us to restrict the spectral transform to grids with iso-latitude rings, see [Grids](@ref).
Furthermore, this term can be written as a fast Fourier transform, if the ``\phi_i`` are
equally spaced on the latitude ring ``j``. Note that the in-ring index ``i`` can depend on
the ring index ``j``, so that one can have reduced grids, which have fewer grid points towards
the poles, for example. Also the Legendre polynomials only have to be computed for the colatitudes
``\theta_j`` (and in fact only one hemisphere, due to the north-south symmetry discussed in the
[Synthesis](@ref synthesis)). It is therefore practical and efficient to design a spectral transform
implementation for ring grids, but there is no need to hardcode a specific grid.

## Spectral packing

Spectral packing is the way how the coefficients ``a_{lm}`` of the spherical harmonics of a given spectral
field are stored in an array. SpeedyWeather.jl uses the conventional spectral packing of degree ``l``
and order ``m`` as illustrated in the following image
[(Cyp, CC BY-SA 3.0, via Wikimedia Commons)](https://commons.wikimedia.org/wiki/File:Rotating_spherical_harmonics.gif)
```@raw html
<img src="https://upload.wikimedia.org/wikipedia/commons/1/12/Rotating_spherical_harmonics.gif">
```
Every row represents an order ``l \geq 0``, starting from ``l=0`` at the top.
Every column represents an order ``m \geq 0``, starting from ``m=0`` on the left. The coefficients of these
spherical harmonics are directly mapped into a matrix ``a_{lm}`` as 

|     |``m``     |          |          |          |
| :-: | :------: | :------: | :------: | :------: | 
|``l``|``a_{00}``|          |          |          |
|     |``a_{10}``|``a_{11}``|          |          |
|     |``a_{20}``|``a_{12}``|``a_{22}``|          |
|     |``a_{30}``|``a_{13}``|``a_{23}``|``a_{33}``|

which is consistently extended for higher degrees and orders. Consequently, all spectral fields are lower-triangular matrices
with complex entries. The upper triangle excluding the diagonal are zero. Note that internally vector fields
include an additional degree, such that ``l_{max} = m_{max} + 1`` (see [Derivatives in spherical coordinates](@ref) for more information).
The harmonics with ``a_{l0}`` (the first column) are also called _zonal_ harmonics as they are constant with longitude ``\phi``.
The harmonics with ``a_{ll}`` (the main diagonal) are also called _sectoral_ harmonics as they essentially split the sphere
into ``2l`` sectors in longitude ``\phi`` without a zero-crossing in latitude.

For correctness it is mentioned here that SpeedyWeather.jl uses a `LowerTriangularMatrix` type to store
the spherical harmonic coefficients. By doing so, the upper triangle is actually *not* explicitly stored
and the data technically unravelled into a vector, but this is hidden as much as possible from the user.
For more details see [`LowerTriangularMatrices`](@ref lowertriangularmatrices).

!!! info "Array indices"
    For a spectral field `a` note that due to Julia's 1-based indexing the coefficient ``a_{lm}`` is obtained via
    `a[l+1,m+1]`. Alternatively, we may index over 1-based `l`,`m` but a comment is usually added for clarification.

[Fortran SPEEDY](https://users.ictp.it/~kucharsk/speedy-net.html) does not use the same spectral packing as
SpeedyWeather.jl. The alternative packing ``l',m'`` therein uses ``l'=m`` and ``m'=l-m`` as summarized in the
following table.

| degree ``l`` | order ``m`` |  ``l'=m`` |  ``m'=l-m`` |
| :----------: | :---------: | :-------: | :---------: |
|0             |0            |0          |0            |
|1             |0            |0          |1            |
|1             |1            |1          |0            |
|2             |0            |0          |2            |
|2             |1            |1          |1            |
|2             |2            |2          |0            |
|3             |0            |0          |3            |
|...           |...          |...        |...          |

This alternative packing uses the top-left triangle of a coefficient matrix, and the degrees and orders from above are
stored at the following indices

|      |``m'``    |          |          |          |
| :--: | :-:      | :-------:| :-------:| :-------:| 
|``l'``|``a_{00}``|``a_{10}``|``a_{20}``|``a_{30}``|
|      |``a_{11}``|``a_{21}``|``a_{31}``|          |
|      |``a_{22}``|``a_{32}``|          |          |
|      |``a_{33}``|          |          |          |

This spectral packing is not used in SpeedyWeather.jl but illustrated here for completeness and comparison with
Fortran SPEEDY.

SpeedyWeather.jl uses triangular truncation such that only spherical harmonics with ``l \leq l_{max}`` and ``|m| \leq m_{max}``
are explicitly represented. This is usually described as ``Tm_{max}``, with ``l_{max} = m_{max}`` (although in vector quantities
require one more degree ``l`` in the recursion relation of meridional gradients). For example, T31 is the spectral resolution
with ``l_{max} = m_{max} = 31``. Note that the degree ``l`` and order ``m`` are mathematically 0-based, such that the
corresponding coefficient matrix is of size 32x32.

## Available horizontal resolutions

Technically, SpeedyWeather.jl supports arbitrarily chosen resolution parameter `trunc` when
creating the [SpectralGrid](@ref) that refers to the highest non-zero degree ``l_{max}``
that is resolved in spectral space. SpeedyWeather.jl will always try to choose an
easily-Fourier transformable[^FFT] size of the grid, but as we use
[FFTW.jl](https://github.com/JuliaMath/FFTW.jl) there is quite some flexibility without
performance sacrifice. However, this has traditionally lead to typical resolutions that
we also use for testing we therefore recommend to use.
They are as follows with more details below

| `trunc`       | nlon | nlat | ``\Delta x`` |
| ------------- | ---- | ---- | ------------ |
| 31 (default)  | 96   | 48   | 400 km       |
| 42            | 128  | 64   | 312 km       |
| 63            | 192  | 96   | 216 km       |
| 85            | 256  | 128  | 165 km       |
| 127           | 384  | 192  | 112 km       |
| 170           | 512  | 256  | 85 km        |
| 255           | 768  | 384  | 58 km        |
| 341           | 1024 | 512  | 43 km        |
| 511           | 1536 | 768  | 29 km        |
| 682           | 2048 | 1024 | 22 km        |
| 1024          | 3072 | 1536 | 14 km        |
| 1365          | 4092 | 2048 | 11 km        |

Some remarks on this table
- This assumes the default quadratic truncation, you can always adapt the grid resolution via the `dealiasing` option, see [Matching spectral and grid resolution](@ref)
- `nlat` refers to the total number of latitude rings, see [Grids](@ref). With non-Gaussian grids, `nlat` will be one one less, e.g. 47 instead of 48 rings.
- `nlon` is the number of longitude points on the [Full Gaussian Grid](@ref FullGaussianGrid), for other grids there will be at most these number of points around the Equator.
- ``\Delta x`` is the horizontal resolution. For a spectral model there are many ways of estimating this[^Randall2021]. We use here the square root of the average area a grid cell covers, see [Effective grid resolution](@ref)

## Effective grid resolution

There are many ways to estimate the effective grid resolution of spectral models[^Randall2021].
Some of them are based on the wavelength a given spectral resolution allows to
represent, others on the total number of real variables per area.
However, as many atmospheric models do represent a considerable amount of physics
on the grid (see [Parameterizations](@ref parameterizations)) there is also a good argument to
include the actual grid resolution into this estimate and not just the spectral
resolution. We therefore use the average grid cell area to estimate the
resolution
```math
\Delta x = \sqrt{\frac{4\pi R^2}{N}}
```
with ``N`` number of grid points over a sphere with radius ``R``. However, we have
to acknowledge that this usually gives higher resolution compared to other methods
of estimating the effective resolution, see [^Randall2021] for a discussion. You may therefore
need to be careful to make claims that, e.g. `trunc=85` can resolve the
atmospheric dynamics at a scale of 165km.

## Derivatives in spherical coordinates

Horizontal gradients in spherical coordinates are defined for a scalar field ``A`` and the latitudes ``\theta``
and longitudes ``\lambda`` as

```math
\nabla A = \left(\frac{1}{R\cos\theta}\frac{\partial A}{\partial \lambda}, \frac{1}{R}\frac{\partial A}{\partial \theta} \right).
```

However, the divergence of a vector field ``\mathbf{u} = (u,v)`` includes additional ``\cos(\theta)`` scalings

```math
\nabla \cdot \mathbf{u} = \frac{1}{R\cos\theta}\frac{\partial u}{\partial \lambda} +
\frac{1}{R\cos\theta}\frac{\partial (v \cos\theta)}{\partial \theta},
```

and similar for the curl

```math
\nabla \times \mathbf{u} = \frac{1}{R\cos\theta}\frac{\partial v}{\partial \lambda} -
\frac{1}{R\cos\theta}\frac{\partial (u \cos\theta)}{\partial \theta}.
```

The radius of the sphere (i.e. Earth) is ``R``. The zonal gradient scales with ``1/\cos(\theta)`` as the 
longitudes converge towards the poles (note that ``\theta`` describes latitudes here, definitions using colatitudes
replace the ``\cos`` with a ``\sin``.)

Starting with a spectral field of vorticity ``\zeta`` and divergence ``\mathcal{D}`` one can obtain stream function ``\Psi``
and velocity potential ``\Phi`` by inverting the Laplace operator ``\nabla^2``:

```math
\Psi = \nabla^{-2}\zeta, \quad \Phi = \nabla^{-2}\mathcal{D}.
```

The velocities ``u,v`` are then obtained from ``(u,v) = \nabla^\bot\Psi + \nabla\Phi`` following the definition from above
and ``\nabla^\bot = (-R^{-1}\partial_\theta, (R\cos\theta)^{-1}\partial_\lambda)``

```math
\begin{aligned}
u &= -\frac{1}{R}\partial_\theta\Psi + \frac{1}{R\cos\theta}\partial_\lambda\Phi \\
v &= +\frac{1}{R}\partial_\theta\Phi + \frac{1}{R\cos\theta}\partial_\lambda\Psi.
\end{aligned}
```

Alternatively, we can use the velocities ``U = u\cos\theta, V = v\cos\theta``, which we do as the meridional gradient
for spherical harmonics is easier implemented with a ``\cos\theta``-scaling included, and because the divergence and 
curl in spherical coordinates evaluates the meridional gradient with ``U,V`` and not ``u,v``. From ``u,v`` we can
return to ``\zeta, \mathcal{D}`` via

```math
\begin{aligned}
\zeta &= \frac{1}{R\cos\theta}\partial_\lambda v - \frac{1}{R\cos\theta}\partial_\theta (u \cos\theta) \\
\mathcal{D} &= \frac{1}{R\cos\theta}\partial_\lambda u + \frac{1}{R\cos\theta}\partial_\theta (v \cos\theta).
\end{aligned}
```

Equivalently, we have

```math
\begin{aligned}
U &= -\frac{\cos\theta}{R}\partial_\theta\Psi + \frac{1}{R}\partial_\lambda\Phi \\
V &= +\frac{\cos\theta}{R}\partial_\theta\Phi + \frac{1}{R}\partial_\lambda\Psi \\
\zeta &= \frac{1}{R}\partial_\lambda \left( \frac{V}{\cos^2\theta} \right) -
\frac{\cos\theta}{R}\partial_\theta \left( \frac{U}{\cos^2\theta} \right) \\
\mathcal{D} &= \frac{1}{R}\partial_\lambda \left( \frac{U}{\cos^2\theta} \right) +
\frac{\cos\theta}{R}\partial_\theta \left( \frac{V}{\cos^2\theta} \right).
\end{aligned}
```

which is a more convenient formulation because of the way how the [Meridional derivative](@ref)
is implemented with a recursion relation, actually computing ``\cos\theta \partial_\theta``
rather than ``\partial_\theta`` directly. The remaining cosine scalings in
``(U,V)*\cos^{-2}\theta`` are done in grid-point space.
If one wanted to get back to ``\zeta, \mathcal{D}`` this is how it would be done, but
it is often more convenient to unscale ``U,V`` on the fly in the spectral transform
to obtain ``u,v`` and then divide again by ``\cos\theta`` when any gradient (or divergence or
curl) is taken. This is because other terms would need that single ``\cos\theta`` unscaling
too before a gradient is taken. How the operators ``\nabla, \nabla \times, \nabla \cdot`` can
be implemented with spherical harmonics is presented in the following sections.

Also note that SpeedyWeather.jl scales the equations with the radius `R` (see [Radius scaling](@ref scaling))
such that the divisions by `R` drop out in this last formulation too.

### Zonal derivative

The zonal derivative of a scalar field ``\Psi`` in spectral space is the zonal derivative of all its respective
spherical harmonics ``\Psi_{lm}(\phi,\theta)`` (now we use ``\phi`` for longitudes to avoid confusion with the
Legendre polynomials ``\lambda_{lm}``)

```math
v_{lm} = \frac{1}{R \cos(\theta)} \frac{\partial}{\partial \phi} \left( \lambda_l^m(\cos\theta) e^{im\phi} \right) =
\frac{im}{R \cos(\theta)} \lambda_l^m(\cos\theta) e^{im\phi} = \frac{im}{R \cos(\theta)} \Psi_{lm}
```

So for every spectral harmonic, ``\cos(\theta)v_{lm}`` is obtained from ``\Psi_{lm}`` via a multiplication
with ``im/R``. Unscaling the ``\cos(\theta)``-factor is done after transforming
the spectral coefficients ``v_{lm}`` into grid-point space. As discussed in [Radius scaling](@ref scaling),
SpeedyWeather.jl scales the stream function as ``\tilde{\Psi} = R^{-1}\Psi`` such that the division by
radius ``R`` in the gradients can be omitted. The zonal derivative becomes therefore effectively for
each spherical harmonic a scaling with its (imaginary) order ``im``. The spherical harmonics are essentially
just a Fourier transform in zonal direction and the derivative a multiplication with the respective wave
number ``m`` times imaginary ``i``.

### Meridional derivative

The meridional derivative of the spherical harmonics is a derivative of the Legendre polynomials for which the following
recursion relation applies[^Randall2021],[^Durran2010],[^GFDL]

```math
\cos\theta \frac{dP_{l,m}}{d\theta} = -l\epsilon_{l+1,m}P_{l+1,m} + (l+1)\epsilon_{l,m}P_{l-1,m}.
```
with recursion factors
```math
\epsilon_{l,m} = \sqrt{\frac{l^2-m^2}{4l^2-1}}
```

In the following we use the example of obtaining the zonal velocity ``u`` from the stream function ``\Psi``,
which is through the negative meridional gradient. For the meridional derivative itself the leading minus sign
has to be omitted. Starting with the spectral expansion

```math
\Psi(\lambda,\theta) = \sum_{l,m}\Psi_{l,m}P_{l,m}(\sin\theta)e^{im\lambda}
```
we multiply with ``-R^{-1}\cos\theta\partial_\theta`` to obtain

```math
\cos\theta\left(-\frac{1}{R}\partial_\theta\Psi \right) = -\frac{1}{R}\sum_{l,m}\Psi_{l,m}e^{im\lambda}\cos\theta\partial_\theta P_{l,m}
```
at which point the recursion from above can be applied. Collecting terms proportional to ``P_{l,m}`` then yields

```math
(\cos(\theta)u)_{l,m} = -\frac{1}{R}(-(l-1)\epsilon_{l,m}\Psi_{l-1,m} + (l+2)\epsilon_{l+1,m}\Psi_{l+1,m})
```

To obtain the coefficient of each spherical harmonic ``l,m`` of the meridional gradient of a spectral field, two 
coefficients at ``l-1,m`` and ``l+1,m`` have to be combined. This means that the coefficient of a gradient
``((\cos\theta) u)_{lm}`` is a linear combination of the coefficients of one higher and one lower degree
``\Psi_{l+1,m},\Psi_{l-1,m}``. As the coefficient ``\Psi_{lm}`` with ``m<l`` are zero, the sectoral harmonics
(``l=m``) of the gradients are obtained from the first off-diagonal only. However, the ``l=l_{max}`` harmonics of
the gradients require the ``l_{max}-1`` as well as the ``l_{max}+1`` harmonics. As a consequence
vector quantities like velocity components ``u,v`` require one more degree ``l`` than scalar quantities like
vorticity[^Bourke72]. However, for easier compatibility all spectral fields in SpeedyWeather.jl use one more
degree ``l``, but scalar quantities should not make use of it. Equivalently, the last degree ``l`` is 
set to zero before the time integration, which only advances scalar quantities.


In SpeedyWeather.jl vector quantities
like ``u,v`` use therefore one more meridional mode than scalar quantities such as vorticity ``\zeta`` or stream
function ``\Psi``. The meridional derivative in SpeedyWeather.jl also omits the ``1/R``-scaling as explained for
the [Zonal derivative](@ref) and in [Radius scaling](@ref scaling).

### Divergence and curl in spherical harmonics

The meridional gradient as described above can be applied to scalars, such as ``\Psi`` and ``\Phi`` in the conversion
to velocities ``(u,v) = \nabla^\bot\Psi + \nabla\Phi``, however, the operators curl ``\nabla \times`` and divergence
``\nabla \cdot`` in spherical coordinates involve a ``\cos\theta`` scaling _before_ the meridional gradient is applied.
How to translate this to spectral coefficients has to be derived separately[^Randall2021],[^Durran2010].

The spectral transform of vorticity ``\zeta`` is
```math
\zeta_{l,m} = \frac{1}{2\pi}\int_{-\tfrac{\pi}{2}}^\tfrac{\pi}{2}\int_0^{2\pi} \zeta(\lambda,\theta)
P_{l,m}(\sin\theta) e^{im\lambda} d\lambda \cos\theta d\theta
```
Given that ``R\zeta = \cos^{-1}\partial_\lambda v - \cos^{-1}\partial_\theta (u \cos\theta)``,
we therefore have to evaluate a meridional integral of the form
```math
\int P_{l,m} \frac{1}{\cos \theta} \partial_\theta(u \cos\theta) \cos \theta d\theta
```
which can be solved through integration by parts. As ``u\cos\theta = 0`` at ``\theta = \pm \tfrac{\pi}{2}`` only the integral
```math
= -\int \partial_\theta P_{l,m} (u \cos\theta) d\theta = -\int \cos\theta \partial_\theta P_{l,m}
(\frac{u}{\cos\theta}) \cos\theta d\theta
```
remains. Inserting the recurrence relation from the [Meridional derivative](@ref) turns this into
```math
= -\int \left(-l \epsilon_{l+1,m}P_{l+1,m} + (l+1)\epsilon_{l,m} P_{l-1,m} \right) (\frac{u}{\cos\theta})
\cos \theta d\theta
```
Now we expand ``(\tfrac{u}{\cos\theta})`` but only the ``l,m`` harmonic will project onto``P_{l,m}``. Let
``u^* = u\cos^{-1}\theta, v^* = v\cos^{-1}\theta`` we then have in total
```math
\begin{aligned}
R\zeta_{l,m} &= imv^*_{l,m} + (l+1)\epsilon_{l,m}u^*_{l-1,m} - l\epsilon_{l+1,m}u^*_{l+1,m} \\
RD_{l,m} &= imu^*_{l,m} - (l+1)\epsilon_{l,m}v^*_{l-1,m} + l\epsilon_{l+1,m}v^*_{l+1,m} \\
\end{aligned}
```
And the divergence ``D`` is similar, but ``(u,v) \to (-v,u)``. We have moved the scaling with the
radius ``R`` directly into ``\zeta,D`` as further described in [Radius scaling](@ref scaling).

### Laplacian

The spectral Laplacian is easily applied to the coefficients ``\Psi_{lm}`` of a spectral field
as the spherical harmonics are eigenfunctions of the Laplace operator ``\nabla^2`` in spherical
coordinates with eigenvalues ``-l(l+1)`` divided by the radius squared ``R^2``, i.e.
``\nabla^2 \Psi`` becomes ``\tfrac{-l(l+1)}{R^2}\Psi_{lm}`` in spectral space. For example,
vorticity ``\zeta`` and streamfunction ``\Psi`` are related by ``\zeta = \nabla^2\Psi``
in the barotropic vorticity model. Hence, in spectral space this is equivalent for every
spectral mode of degree ``l`` and order ``m`` to

```math
\zeta_{l,m} = \frac{-l(l+1)}{R^2}\Psi_{l,m}
```

This can be easily inverted to obtain the stream function ``\Psi`` from vorticity ``\zeta`` instead. In order to avoid
division by zero, we set ``\Psi_{0,0}`` here, given that the stream function is only defined up to a constant anyway.
```math
\begin{aligned}
\Psi_{l,m} &= \frac{R^2}{-l(l+1)}\zeta_{l,m} \quad \forall~l,m > 0,\\
\Psi_{0,0} &= 0.
\end{aligned}
```

See also [Horizontal diffusion](@ref diffusion) and [Normalization of diffusion](@ref).

### U,V from vorticity and divergence

After having discussed the zonal and meridional derivatives with spherical harmonics as well as the Laplace operator,
we can derive the conversion from vorticity ``\zeta`` and divergence ``D`` (which are prognostic variables) to
``U=u\cos\theta, V=v\cos\theta``. Both are linear operations that act either solely on a given harmonic
(the zonal gradient and the Laplace operator) or are linear combinations between one lower and one higher degree ``l``
(the meridional gradient). It is therefore computationally more efficient to compute ``U,V`` directly from ``\zeta,D``
instead of calculating stream function and velocity potential first. In total we have

```math
\begin{aligned}
U_{l,m} &= -\frac{im}{l(l+1)}(RD)_{l,m} + \frac{\epsilon_{l+1,m}}{l+1}(R\zeta)_{l+1,m} -
\frac{\epsilon_{l,m}}{l}(R\zeta)_{l-1,m} \\
V_{l,m} &= -\frac{im}{l(l+1)}(R\zeta)_{l,m} - \frac{\epsilon_{l+1,m}}{l+1}(RD)_{l+1,m} +
\frac{\epsilon_{l,m}}{l}(RD)_{l-1,m} \\
\end{aligned}
```

We have moved the scaling with the radius ``R`` directly into ``\zeta,D``
as further described in [Radius scaling](@ref scaling).

## References

[^Malardel2016]: Malardel S, Wedi N, Deconinck W, Diamantakis M, Kühnlein C, Mozdzynski G, Hamrud M, Smolarkiewicz P. A new grid for the IFS. ECMWF newsletter. 2016;146(23-28):321. doi: [10.21957/zwdu9u5i](https://doi.org/10.21957/zwdu9u5i)
[^Gorski2004]: Górski, Hivon, Banday, Wandelt, Hansen, Reinecke, Bartelmann, 2004. _HEALPix: A FRAMEWORK FOR HIGH-RESOLUTION DISCRETIZATION AND FAST ANALYSIS OF DATA DISTRIBUTED ON THE SPHERE_, The Astrophysical Journal. doi:[10.1086/427976](https://doi.org/10.1086/427976)
[^Willmert2020]: Justin Willmert, 2020. [justinwillmert.com](https://justinwillmert.com/)
    - [Introduction to Associated Legendre Polynomials (Legendre.jl Series, Part I)](https://justinwillmert.com/articles/2020/introduction-to-associated-legendre-polynomials/)
    - [Calculating Legendre Polynomials (Legendre.jl Series, Part II)](https://justinwillmert.com/articles/2020/calculating-legendre-polynomials/)
    - [Pre-normalizing Legendre Polynomials (Legendre.jl Series, Part III)](https://justinwillmert.com/articles/2020/pre-normalizing-legendre-polynomials/)
    - [Maintaining numerical accuracy in the Legendre recurrences (Legendre.jl Series, Part IV)](https://justinwillmert.com/articles/2020/maintaining-numerical-accuracy-in-the-legendre-recurrences/)
    - [Introducing Legendre.jl (Legendre.jl Series, Part V)](https://justinwillmert.com/articles/2020/introducing-legendre.jl/)
    - [Numerical Accuracy of the Spherical Harmonic Recurrence Coefficient (Legendre.jl Series Addendum)](https://justinwillmert.com/posts/2020/pre-normalizing-legendre-polynomials-addendum/)
    - [Notes on Calculating the Spherical Harmonics](https://justinwillmert.com/articles/2020/notes-on-calculating-the-spherical-harmonics)
    - [More Notes on Calculating the Spherical Harmonics: Analysis of maps to harmonic coefficients](https://justinwillmert.com/articles/2022/more-notes-on-calculating-the-spherical-harmonics/)
[^Daley78]:  Roger Daley & Yvon Bourassa (1978) Rhomboidal versus triangular spherical harmonic truncation: Some verification statistics, Atmosphere-Ocean, 16:2, 187-196, DOI: [10.1080/07055900.1978.9649026](https://doi.org/10.1080/07055900.1978.9649026)
[^Randall2021]: David Randall, 2021. [An Introduction to Numerical Modeling of the Atmosphere](http://hogback.atmos.colostate.edu/group/dave/at604pdf/An_Introduction_to_Numerical_Modeling_of_the_Atmosphere.pdf), Chapter 22.
[^Durran2010]: Dale Durran, 2010. [Numerical Methods for Fluid Dynamics](https://link.springer.com/book/10.1007/978-1-4419-6412-0), Springer. In particular section 6.2, 6.4.
[^GFDL]: Geophysical Fluid Dynamics Laboratory, [The barotropic vorticity equation](https://www.gfdl.noaa.gov/wp-content/uploads/files/user_files/pjp/barotropic.pdf).
[^FFT]: Depending on the implementation of the Fast Fourier Transform ([Cooley-Tukey algorithm](https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm), or or the [Bluestein algorithm](https://en.wikipedia.org/wiki/Chirp_Z-transform#Bluestein.27s_algorithm)) *easily Fourier-transformable* can mean different things: Vectors of the length ``n`` that is a power of two, i.e. ``n = 2^i`` is certainly easily Fourier-transformable, but for most FFT implementations so are ``n = 2^i3^j5^k`` with ``i,j,k`` some positive integers. In fact, [FFTW](http://fftw.org/) uses ``O(n \log n)`` algorithms even for prime sizes.
[^Bourke72]: Bourke, W. An Efficient, One-Level, Primitive-Equation Spectral Model. Mon. Wea. Rev. 100, 683–689 (1972). doi:[10.1175/1520-0493(1972)100<0683:AEOPSM>2.3.CO;2](https://doi.org/10.1175/1520-0493(1972)100<0683:AEOPSM>2.3.CO;2)
