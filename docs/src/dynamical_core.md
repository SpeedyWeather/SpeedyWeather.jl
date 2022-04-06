# Dynamical core

A mathematical and implementation-specific description of the dynamical core used in SpeedyWeather.jl

## Mathematical background

The [primitive equations](https://en.wikipedia.org/wiki/Primitive_equations) solved by SpeedyWeather.jl are

```math
\begin{aligned}
\partial_t u = ... \\
\partial_t v = ... \\
\partial_t T = ... \\ 
\partial_t Q = ... \\
\end{aligned}
```

more to come

## Horizontal diffusion

The horizontal diffusion in SpeedyWeather.jl is implemented as an ``n``-th power Laplacian ``\nabla^{2n}`` in spectral space
that is applied implicitly.

Given that the spherical harmonics are eigenfunctions of the Laplace operator ``\nabla^2`` in spherical coordinates with
eigenvalues ``-l(l+1)``, applying hyperdiffusion (hyper when ``n>1``) can be implemented as a multiplication of the
spectral coefficients ``a_{lm}`` with ``(-l(l+1))^n``. For example, vorticity ``\zeta`` and streamfunction ``\Psi`` are related
by
```math
\zeta = \nabla^2\Psi
```
Hence, in spectral space this is equivalent for every spectral mode of degree ``l`` and order ``m`` to
```math
\zeta_{l,m} = \frac{-l(l+1)}{R^2}\Psi_{l,m}
```
with ``R`` the radius of the sphere (i.e. Earth). It is therefore computationally not more expensive to apply hyperdiffusion
as the ``(-l(l+1))^nR^{-2n}`` can be precomputed. In SpeedyWeather.jl the diffusion is applied implicitly.
For that, consider a [leapfrog scheme](https://en.wikipedia.org/wiki/Leapfrog_integration)
with time step ``\Delta t`` of variable ``\zeta`` to obtain from time steps ``i-1`` and ``i``, the next time step ``i+1``
```math
\zeta_{i+1} = \zeta_{i-1} + 2\Delta t d\zeta,
```
with ``d\zeta`` being some tendency evaluated from ``\zeta_i``. Now we want to add a diffusion term ``-\nu \nabla^{2n}\zeta``
with viscosity ``\nu``, wich however, is implicitly calculated from ``\zeta_{i+1}``, then
```math
\zeta_{i+1} = \zeta_{i-1} + 2\Delta t (d\zeta - \nu\nabla^{2n}\zeta_{i+1})
```
As the application of ``\nu\nabla^{2n}`` is, for every spectral mode, equivalent to a multiplication of
a constant, we can rewrite this to
```math
\zeta_{i+1} = \frac{\zeta_{i-1} + 2\Delta t d\zeta}{1 + 2\Delta \nu\nabla^{2n}},
```
and expand the numerator to
```math
\zeta_{i+1} = \zeta_{i-1} + 2\Delta t \frac{d\zeta - \nu\nabla^{2n}\zeta_{i-1}}{1+2\Delta t \nu \nabla^{2n}},
```
Hence the diffusion can be applied implicitly by updating the tendency ``d\zeta`` as
```math
d\zeta \to \frac{d\zeta - \nu\nabla^{2n}\zeta_{i-1}}{1+2\Delta t \nu \nabla^{2n}}
```
which only depends on ``\zeta_{i-1}``. Now let ``D_\text{explicit} = \nu\nabla^{2n}`` be the explicit part and
``D_\text{implicit} = 1 + 2\Delta t \nu\nabla^{2n}`` the implicit part. Both parts can be precomputed and are
only an element-wise multiplication in spectral space. For every spectral harmonic ``l,m`` we do
```math
d\zeta \to D_\text{implicit}^{-1}(d\zeta - D_\text{explicit}\zeta_{i-1}).
```
Hence 2 multiplications and 1 subtraction with precomputed constants.
However, we will normalize the (hyper-)Laplacians as described in the following.

### Normalization of diffusion

In physics, the Laplace operator ``\nabla^2`` is often used to represent diffusion due to viscosity in a fluid. In that case,
the viscosity coefficient is ``\nu`` of units ``\text{m}^2\text{s}^{-1}`` and the full operator reads as ``\nu \nabla^2`` with units
``(\text{m}^2\text{s}^{-1})(\text{m}^{-2}) = \text{s}^{-1}``. This motivates us to normalize the Laplace operator by a constant
of units ``\text{m}^{-2}`` and the viscosity coefficient by its inverse such that the viscosity coefficient becomes a
damping timescale of unit ``\text{s}^{-1}``. Given the application in spectral space we decide to normalize by the
largest eigenvalue ``l_\text{max}(l_\text{max}+1)`` such that all entries in the discrete spectral Laplace operator are
in ``[0,1]``. The normalized viscosity coefficient ``\nu^* = l_\text{max}(l_\text{max}+1)\nu`` is therefore reinterpreted
as the time scale at which the highest wavenumber is dampened to zero due to diffusion. Together we have 
```math
D^\text{explicit}_{l,m} = \nu^* \frac{l(l+1)}{l_\text{max}(l_\text{max}+1)}
```
and the hyper-Laplacian of power ``n`` follows as
```math
D^\text{explicit,n}_{l,m} = \nu^* \left(\frac{l(l+1)}{l_\text{max}(l_\text{max}+1)}\right)^n
```
and the implicit part is accordingly ``D^\text{implicit,n}_{l,m} = 1 + 2\Delta t D^\text{explicit,n}_{l,m}``.

## Time integration

SpeedyWeather.jl uses a leapfrog time scheme with a Robert's and William's filter
to dampen the computational mode and achieve 3rd order accuracy.

### Oscillation equation

```math
\frac{dF}{dt} = i\omega F
```