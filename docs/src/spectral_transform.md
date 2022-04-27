# Spherical Harmonic Transform

The following sections outline the implementation of the spherical harmonic transform (in short _spectral_ transform)
between the coefficients of the spherical harmonics (the _spectral_ space) and the grid space on a longitude-latitude
[regular Gaussian grid](https://confluence.ecmwf.int/display/FCST/Gaussian+grids).

## Inspiration

The spectral transform implemented by SpeedyWeather.jl follows largely Justin Willmert's
[CMB.jl](https://github.com/jmert/CMB.jl) package and makes use of
[AssociatedLegendrePolynomials.jl](https://github.com/jmert/AssociatedLegendrePolynomials.jl) and
[FFTW.jl](https://github.com/JuliaMath/FFTW.jl)/[FastTransforms.jl](https://github.com/JuliaApproximation/FastTransforms.jl) for the Fourier transform. Justin described his work in a Blog series [^1],[^2],[^3],[^4],[^5],[^6],[^7],[^8]

## Spherical harmonics

The [spherical harmonics](https://en.wikipedia.org/wiki/Spherical_harmonics) ``Y_{lm}`` of degree ``l`` and order ``m``
over the longitude ``\theta = (0,2\pi)`` and latitude ``\phi = (-\tfrac{\pi}{2},\tfrac{\pi}{2})`` (or
using colatitudes ``\phi = (0,\pi)``), are

```math
Y_{lm}(\theta,\phi) = \lambda_l^m(\cos\theta) e^{im\phi}
```

with ``\lambda_l^m`` being the pre-normalized associated Legendre polynomials, and ``e^{im\phi}`` are the
complex exponentials (the Fourier modes). Together they form a set of orthogonal basis functions on the sphere.
For an interactive visualisation of the spherical harmonics, see
[here](https://justinwillmert.com/posts/2020/plots-of-the-spherical-harmonics-eigenmodes/).

## Synthesis or inverse spectral transform (spectral to grid)

```math
f(\theta,\phi) = \sum_{l=0}^{l_{max}} \sum_{m=-l}^l a_{lm} Y_{lm}(\theta,\phi)
```

## Analysis or forward spectral transform (grid to spectral)

```math
\hat{a}_{lm} = \sum_{i=1}^N f(\theta_i,\phi_i) Y_{lm}(\theta_i,\phi_i) \sin(\theta_i) \Delta \theta_i \Delta \phi_i
```

## Spectral packing

Conventional packing ``l,m`` versus alternative packing ``l',m'`` and arbitrary numbering ``i``.

| degree ``l`` | order ``m`` |  ``l'=m`` |  ``m'=l-m`` | ``i`` |
|-|-|-|-|-|
|0|0|0|0|1|
|1|0|0|1|2|
|1|1|1|0|3|
|2|0|0|2|4|
|2|1|1|1|5|
|2|2|2|0|6|
|3|0|0|3|7|
|...|...|...|...|...|

Degree $l$, order $m$

| |``m``| | | |
|-|-|-|-|-|
|l|1| | | |
| |2|3| | |
| |4|5|6| |
| |7|8|9|10|

Alternative packing

| |``m'``| | | |
|-|-|-|-|-|
|``l'``|1|2|4|7|
| |3|5|8| |
| |6|9| | |
| |10| | | |

## Examples

```julia
julia> using SpeedyWeather
julia> alms = zeros(ComplexF64,3,3)    # spectral coefficients
julia> alms[2,2] = 1                   # only l=1,m=1 Legendre polynomial
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

and we have reobtained the ``l=m=1`` spheri

## References

[^1]: Justin Willmert, 2020. [Introduction to Associated Legendre Polynomials (Legendre.jl Series, Part I)](https://justinwillmert.com/articles/2020/introduction-to-associated-legendre-polynomials/)  
[^2]: Justin Willmert, 2020. [Calculating Legendre Polynomials (Legendre.jl Series, Part II)](https://justinwillmert.com/articles/2020/calculating-legendre-polynomials/)  
[^3]: Justin Willmert, 2020. [Pre-normalizing Legendre Polynomials (Legendre.jl Series, Part III)](https://justinwillmert.com/articles/2020/pre-normalizing-legendre-polynomials/)  
[^4]: Justin Willmert, 2020. [Maintaining numerical accuracy in the Legendre recurrences (Legendre.jl Series, Part IV)](https://justinwillmert.com/articles/2020/maintaining-numerical-accuracy-in-the-legendre-recurrences/)  
[^5]: Justin Willmert, 2020. [Introducing Legendre.jl (Legendre.jl Series, Part V)](https://justinwillmert.com/articles/2020/introducing-legendre.jl/)  
[^6]: Justin Willmert, 2020. [Numerical Accuracy of the Spherical Harmonic Recurrence Coefficient (Legendre.jl Series Addendum)](https://justinwillmert.com/posts/2020/pre-normalizing-legendre-polynomials-addendum/)  
[^7]: Justin Willmert, 2020. [Notes on Calculating the Spherical Harmonics](https://justinwillmert.com/articles/2020/notes-on-calculating-the-spherical-harmonics)  
[^8]: Justin Willmert, 2022. [More Notes on Calculating the Spherical Harmonics: Analysis of maps to harmonic coefficients](https://justinwillmert.com/articles/2022/more-notes-on-calculating-the-spherical-harmonics/)  