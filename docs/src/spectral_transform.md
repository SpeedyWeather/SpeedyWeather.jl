# Spherical Harmonic Transform

The following sections outline the implementation of the spherical harmonic transform (in short _spectral_ transform)
between the coefficients of the spherical harmonics (the _spectral_ space) and the grid space on a longitude-latitude
[regular Gaussian grid](https://confluence.ecmwf.int/display/FCST/Gaussian+grids).

## Inspiration

The spectral transform implemented by SpeedyWeather.jl follows largely Justin Willmert's
[CMB.jl](https://github.com/jmert/CMB.jl) package and makes use of
[AssociatedLegendrePolynomials.jl](https://github.com/jmert/AssociatedLegendrePolynomials.jl) and
[FFTW.jl](https://github.com/JuliaMath/FFTW.jl) (for `Float32/64`) or [GenericFFT.jl](https://github.com/JuliaApproximation/GenericFFT.jl) (for generic) for the Fourier transform. Justin described his work in a Blog series [^1][^2][^3][^4][^5][^6][^7][^8].

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
    The implementations of the spherical transforms in SpeedyWeather.jl use colatitudes ``\theta = (0,\pi)``
    (0 at the north pole) but the dynamical core uses latitudes ``\theta = (-\pi/2,\pi/2)`` (``\pi/2`` at the north pole).
    However, all arrays are always sorted north to south such that `[i,1]` will access the northern-most grid points.
    Note: We may also use latitudes in the spherical harmonic transfom in the future for consistency. 

## Synthesis (spectral to grid)

The synthesis (or inverse transform) takes the spectral coefficients ``a_{lm}`` and transforms them to grid-point values
``f(\phi,\theta)`` (for the sake of simplicity first regarded as continuous). The synthesis is a linear combination of
the spherical harmonics ``Y_{lm}`` with non-zero coefficients.

```math
f(\phi,\theta) = \sum_{l=0}^{\infty} \sum_{m=-l}^l a_{lm} Y_{lm}(\phi,\theta)
```

We obtain an approximation with a finite set of ``a_{l,m}`` by truncating the series after ``l = l_{max}``.

## Analysis (grid to spectral)

Starting in grid-point space we can transform a field ``f(\lambda,\theta)`` into the spectral space of the spherical harmonics by

```math
a_{l,m} = \int_0^{2\pi} \int_{-\tfrac{\pi}{2}}^\tfrac{\pi}{2} f(\lambda,\theta) Y_{l,m}(\lambda,\theta) \cos \theta d\theta d\lambda
```

This integral has to be discretized to when grid-point values ``f(\lambda_i,\theta_i)`` are used. For more details, see [^7],[^8].


## Spectral packing

Spectral packing is the way how the coefficients ``a_{lm}`` of the spherical harmonics of a given spectral field are
stored in an array. SpeedyWeather.jl uses the conventional spectral packing of degree ``l`` and order ``m`` as
illustrated in the following image
[(Cyp, CC BY-SA 3.0, via Wikimedia Commons)](https://commons.wikimedia.org/wiki/File:Rotating_spherical_harmonics.gif)
```@raw html
<img src="https://upload.wikimedia.org/wikipedia/commons/1/12/Rotating_spherical_harmonics.gif">
```
Every row represents an order ``l \geq 0``, starting from ``l=0`` at the top. Every column represents an order
``m \geq 0``, starting from ``m=0`` on the left. The coefficients of these spherical harmonics are directly
mapped into a matrix ``a_{lm}`` as 

|     |``m``     |          |          |          |
| :-: | :------: | :------: | :------: | :------: | 
|``l``|``a_{00}``|          |          |          |
|     |``a_{10}``|``a_{11}``|          |          |
|     |``a_{20}``|``a_{12}``|``a_{22}``|          |
|     |``a_{30}``|``a_{13}``|``a_{23}``|``a_{33}``|

which is consistently extended for higher degrees and orders. Consequently, all spectral fields are lower-triangular matrices
with complex entries. The upper triangle excluding the diagonal explicitly stores zeros. Note that internally vector fields
include an additional degree, such that ``l_{max} = m_{max} + 1`` (see [Gradients in spherical coordinates](@ref) for more information).
The harmonics with ``a_{l0}`` (the first column) are also called _zonal_ harmonics as they are constant with longitude ``\phi``.
The harmonics with ``a_{ll}`` (the main diagonal) are also called _sectoral_ harmonics as they essentially split the sphere
into ``2l`` sectors in longitude ``\phi`` without a zero-crossing in latitude.

!!! info "Array indices"
    For a spectral field `alms` note that due to Julia's 1-based indexing the coefficient ``a_{lm}`` is obtained via
    `alms[l+1,m+1]`.

Fortran speedy does not use the same spectral packing as SpeedyWeather.jl. The alternative packing ``l',m'`` therein
uses ``l'=m`` and ``m'=l-m`` as summarized in the following table.

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
Fortran-speedy.

### Example transforms

```julia
julia> using SpeedyWeather
julia> alms = zeros(ComplexF64,3,3)    # spectral coefficients
julia> alms[2,2] = 1                   # only l=1,m=1 harmonic
julia> map = gridded(alms)             # convert to grid space
8×4 Matrix{Float64}:
 -0.324541  -0.600363  -0.600363  -0.324541
 -0.134429  -0.248678  -0.248678  -0.134429
  0.134429   0.248678   0.248678   0.134429
  0.324541   0.600363   0.600363   0.324541
  0.324541   0.600363   0.600363   0.324541
  0.134429   0.248678   0.248678   0.134429
 -0.134429  -0.248678  -0.248678  -0.134429
 -0.324541  -0.600363  -0.600363  -0.324541
 
julia> spectral(map)                   # back to spectral space
3×3 Matrix{ComplexF64}:
 0.0+0.0im  0.0+0.0im          0.0+0.0im
 0.0+0.0im  1.0+3.60727e-17im  0.0+0.0im
 0.0+0.0im  0.0+0.0im          0.0+0.0im
```

and we have successfully reobtained the ``l=m=1`` spherical harmonic.

## Available horizontal resolutions

SpeedyWeather.jl uses triangular truncation such that only spherical harmonics with ``l \leq l_{max}`` and ``|m| \leq m_{max}``
are explicitly represented. This is usually described as ``Tm_{max}``, with ``l_{max} = m_{max}`` (although in vector quantities
require one more degree ``l`` in the recursion relation of meridional gradients). For example, T31 is the spectral resolution
with ``l_{max} = m_{max} = 31``. Note that the degree ``l`` and order ``m`` are mathematically 0-based, such that the
corresponding coefficient matrix is of size 32x32.

Using triangular truncation[^9], there are constraints on the corresponding grid resolution. Let `nlon`, `nlat` be the number of
longitudes, latitudes on a regular Gaussian grid. Then spectral and grid resolution have to be chosen such that

- ``nlon \geq 3l_{max}+1``
- ``nlat \geq (3l_{max}+1)/2``

In general, we choose ``nlon = 2nlat``, and ideally ``nlon`` is easily Fourier-transformable, e.g. ``nlon = 2^i3^j5^k`` with some
integers ``i,j,k \geq 0``. SpeedyWeather.jl is tested at the following horizontal resolutions, with
``\Delta x = \tfrac{2\pi R}{nlon}`` as the approximate grid spacing at the Equator

| ``l_{max}``   | nlon | nlat | ``\Delta x`` |
| ------------- | ---- | ---- | ------------ |
| 31 (default)  | 96   | 48   | 400 km       |
| 42            | 128  | 64   | 300 km       |
| 85            | 256  | 128  | 160 km       |
| 170           | 512  | 256  | 80 km        |
| 341           | 1024 | 512  | 40 km        |
| 682           | 2048 | 1024 | 20 km        |

Choosing `trunc` as argument in `run_speedy` will automatically choose `nlon`,`nlat` as presented in the table.
Other common choices are T63 (192x96), T127 (384x192), T255 (768x384), T511 (1536x768), among others.

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
longitudes converge towards the poles (note that ``\theta`` describes latitudes here, defintions using colatitudes
replace the ``\cos`` with a ``\sin``.)

Starting with a spectral field of vorticity ``\zeta`` and divergence ``\mathcal{D}`` one can obtain stream function ``\Psi``
and velocity potential ``\Phi`` by inverting the Laplace operator ``\nabla^2``:

```math
\Psi = \nabla^{-2}\zeta, \quad \Phi = \nabla^{-2}\mathcal{D}.
```

The velocities ``u,v`` are then obtained from ``(u,v) = \nabla^\bot\Psi + \nabla\Phi`` following the defintion from above
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

which is a more convenient formulation as required ``\cos\theta`` scalings are reduced to a minimum.
The remaining ``(U,V)*\cos^{-2}\theta`` are done in grid-point space and usually in combination with
other operations like the computation of the vorticity flux. But also note that SpeedyWeather.jl scales
the equations with the radius `R` (see [Radius scaling](@ref)) such that the divisions by `R` drop out too.
As described in [Meridional derivative](@ref), it is more convenient to implement ``\cos\theta \partial_\theta``
via a recursion relation for the Legendre polynomials than ``\partial_\theta`` directly.
How the operators ``\nabla, \nabla \times, \nabla \cdot`` can be implemented with spherical harmonics is
presented in the following sections.

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
the spectral coefficients ``v_{lm}`` into grid-point space. As discussed in [Radius scaling](@ref), SpeedyWeather.jl
scales the stream function as ``\tilde{\Psi} = R^{-1}\Psi`` such that the division by radius ``R`` in the gradients
can be omitted. The zonal derivative becomes therefore effectively for each spherical harmonic a scaling with its
(imaginary) order ``im``. The spherical harmonics are essentially just a Fourier transform
in zonal direction and the derivative a multiplication with the respective wave number ``m`` times imaginary ``i``.

### Meridional derivative

The meridioinal derivative of the spherical harmonics is a derivative of the Legendre polynomials for which the following
recursion relation applies[^10],[^11]

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
the gradients require the ``l_{max}-1`` as well as the ``l_{max}+1`` harmonics. In SpeedyWeather.jl vector quantitie
like ``u,v`` use therefore one more meridional mode than scalar quantities such as vorticity ``\zeta`` or stream
function ``\Psi``. The meridional derivative in SpeedyWeather.jl also omits the ``1/R``-scaling as explained for
the [Zonal derivative](@ref) and in [Radius scaling](@ref).

### Divergence and curl in spherical harmonics

The meridional gradient as described above can be applied to scalars, such as ``\Psi`` and ``\Phi`` in the conversion
to velocities ``(u,v) = \nabla^\bot\Psi + \nabla\Phi``, however, the operators curl ``\nabla \times`` and divergence
``\nabla \cdot`` in spherical coordinates involve a ``\cos\theta`` scaling _before_ the meridional gradient is applied.
How to translate this to spectral coefficients has to be derived separately[^10],[^11].

The spectral transform of vorticity ``\zeta`` is
```math
\zeta_{l,m} = \frac{1}{2\pi}\int_{-\tfrac{\pi}{2}}^\tfrac{\pi}{2}\int_0^{2\pi} \zeta(\lambda,\theta) P_{l,m}(\sin\theta) e^{im\lambda} d\lambda \cos\theta d\theta
```
Given that ``R\zeta = \cos^{-1}\partial_\lambda v - \cos^{-1}\partial_\theta (u \cos\theta)``, we therefore have to evaluate a meridional integral of the form
```math
\int P_{l,m} \frac{1}{\cos \theta} \partial_\theta(u \cos\theta)) \cos \theta d\theta
```
which can be solved through integration by parts. As ``u\cos\theta = 0`` at ``\theta = \pm \tfrac{\pi}{2}`` only the integral
```math
= -\int \partial_\theta P_{l,m} (u \cos\theta) d\theta = -\int \cos\theta \partial_\theta P_{l,m} (\frac{u}{\cos\theta}) \cos\theta d\theta
```
remains. Inserting the recurrence relation from the [Meridional derivative](@ref) turns this into
```math
= -\int \left(-l \epsilon_{l+1,m}P_{l+1,m} + (l+1)\epsilon_{l,m} P_{l-1,m} \right) (\frac{u}{\cos\theta}) \cos \theta d\theta
```
Now we expand ``(\tfrac{u}{\cos\theta})`` but only the ``l,m`` harmonic will project onto``P_{l,m}``. Let
``u^* = u\cos^{-1}\theta, v^* = v\cos^{-1}\theta`` we then have in total
```math
\begin{aligned}
R\zeta_{l,m} &= imv^*_{l,m} + (l+1)\epsilon_{l,m}u^*_{l-1,m} - l\epsilon_{l+1,m}u^*_{l+1,m} \\
RD_{l,m} &= imu^*_{l,m} - (l+1)\epsilon_{l,m}v^*_{l-1,m} + l\epsilon_{l+1,m}v^*_{l+1,m} \\
\end{aligned}
```
And the divergence ``D`` is similar, but ``(u,v) \to (-v,u)``. We have moved the scaling with the radius ``R`` directly into ``\zeta,D`` as further described in [Radius scaling](@ref).

### Laplacian

The spectral Laplacian is easily applied to the coefficients ``\Psi_{lm}`` of a spectral field
as the spherical harmonics are eigenfunctions of the Laplace operator ``\nabla^2`` in spherical coordinates with
eigenvalues ``-l(l+1)`` divided by the radius squared ``R^2``, i.e. ``\nabla^2 \Psi`` becomes ``\tfrac{-l(l+1)}{R^2}\Psi_{lm}``
in spectral space. For example, vorticity ``\zeta`` and streamfunction ``\Psi`` are related by ``\zeta = \nabla^2\Psi``
in the barotropic vorticity model. Hence, in spectral space this is equivalent for every spectral mode of
degree ``l`` and order ``m`` to

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

See also [Horizontal diffusion](@ref) and [Normalization of diffusion](@ref).

### U,V from vorticity and divergence

After having discussed the zonal and meridional derivatives with spherical harmonics as well as the Laplace operator,
we can derive the conversion from vorticity ``\zeta`` and divergence ``D`` (which are prognostic variables) to
``U=u\cos\theta, V=v\cos\theta``. Both are linear operations that act either solely on a given harmonic
(the zonal gradient and the Laplace operator) or are linear combinations between one lower and one higher degree ``l``
(the meridional gradient). It is therefore computationally more efficient to compute ``U,V`` directly from ``\zeta,D``
instead of calculating stream function and velocity potential first. In total we have

```math
\begin{aligned}
U_{l,m} &= -\frac{im}{l(l+1)}(RD)_{l,m} + \frac{\epsilon_{l+1,m}}{l+1}(R\zeta)_{l+1,m} - \frac{\epsilon_{l,m}}{l}(R\zeta)_{l-1,m} \\
V_{l,m} &= -\frac{im}{l(l+1)}(R\zeta)_{l,m} - \frac{\epsilon_{l+1,m}}{l+1}(RD)_{l+1,m} + \frac{\epsilon_{l,m}}{l}(RD)_{l-1,m} \\
\end{aligned}
```

We have moved the scaling with the radius ``R`` directly into ``\zeta,D`` as further described in [Radius scaling](@ref).

## References
[^1]: Justin Willmert, 2020. [Introduction to Associated Legendre Polynomials (Legendre.jl Series, Part I)](https://justinwillmert.com/articles/2020/introduction-to-associated-legendre-polynomials/)
[^2]: Justin Willmert, 2020. [Calculating Legendre Polynomials (Legendre.jl Series, Part II)](https://justinwillmert.com/articles/2020/calculating-legendre-polynomials/)
[^3]: Justin Willmert, 2020. [Pre-normalizing Legendre Polynomials (Legendre.jl Series, Part III)](https://justinwillmert.com/articles/2020/pre-normalizing-legendre-polynomials/)
[^4]: Justin Willmert, 2020. [Maintaining numerical accuracy in the Legendre recurrences (Legendre.jl Series, Part IV)](https://justinwillmert.com/articles/2020/maintaining-numerical-accuracy-in-the-legendre-recurrences/)
[^5]: Justin Willmert, 2020. [Introducing Legendre.jl (Legendre.jl Series, Part V)](https://justinwillmert.com/articles/2020/introducing-legendre.jl/)
[^6]: Justin Willmert, 2020. [Numerical Accuracy of the Spherical Harmonic Recurrence Coefficient (Legendre.jl Series Addendum)](https://justinwillmert.com/posts/2020/pre-normalizing-legendre-polynomials-addendum/)
[^7]: Justin Willmert, 2020. [Notes on Calculating the Spherical Harmonics](https://justinwillmert.com/articles/2020/notes-on-calculating-the-spherical-harmonics)
[^8]: Justin Willmert, 2022. [More Notes on Calculating the Spherical Harmonics: Analysis of maps to harmonic coefficients](https://justinwillmert.com/articles/2022/more-notes-on-calculating-the-spherical-harmonics/)
[^9]: David Randall, 2021. [An Introduction to Numerical Modeling of the Atmosphere](http://hogback.atmos.colostate.edu/group/dave/at604pdf/An_Introduction_to_Numerical_Modeling_of_the_Atmosphere.pdf), Chapter 22.
[^10]: Dale Durran, 2010. [Numerical Methods for Fluid Dynamics](https://link.springer.com/book/10.1007/978-1-4419-6412-0), Springer. In particular section 6.2, 6.4.
[^11]: Geophysical Fluid Dynamics Laboratory, [The barotropic vorticity equation](https://www.gfdl.noaa.gov/wp-content/uploads/files/user_files/pjp/barotropic.pdf).