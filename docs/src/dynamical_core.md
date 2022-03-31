# Dynamical core

A mathematical and implementation-specific description of the dynamical core used in SpeedyWeather.jl

## Mathematical background

The [primitive equations](https://en.wikipedia.org/wiki/Primitive_equations) solved by SpeedyWeather.jl are

$
\begin{aligned}
\partial_t u = ... \\
\partial_t v = ... \\
\partial_t T = ... \\ 
\partial_t Q = ... \\
\end{aligned}
$

more to come

## Horizontal diffusion

The horizontal diffusion in SpeedyWeather.jl is implemented as an $n$-th power Laplacian $\nabla^{2n}$ in spectral space
that is applied implicitly.

Given that the spherical harmonics are eigenfunctions of the Laplace operator $\nabla^2$ in spherical coordinates with
eigenvalues $-l(l+1)$, applying hyperdiffusion (hyper when $n>1$) can be implemented as a multiplication of the
spectral coefficients $a_{lm}$ with $(-l(l+1))^n$. For example, vorticity $\zeta$ and streamfunction $\Psi$ are related
by
$
\zeta = \nabla^2\Psi
$
Hence, in spectral space this is equivalent for every spectral mode of degree $l$ and order $m$ to
$
\zeta_{l,m} = \frac{-l(l+1)}{R^2}\Psi_{l,m}
$
with $R$ the radius of the sphere (i.e. Earth). It is therefore computationally not more expensive to apply hyperdiffusion
as the $(-l(l+1))^nR^{-2n}$ can be precomputed. In SpeedyWeather.jl the diffusion is applied implicitly.
For that, consider a [leapfrog scheme](https://en.wikipedia.org/wiki/Leapfrog_integration)
with time step $\Delta t$ of variable $\zeta$ to obtain from time steps $i-1$ and $i$, the next time step $i+1$
$
\zeta_{i+1} = \zeta_{i-1} + 2\Delta t d\zeta,
$ 
with $d\zeta$ being some tendency evaluated from $\zeta_i$. Now we want to add a diffusion term $-\nu \nabla^{2n}\zeta$
with viscosity $\nu$, wich however, is implicitly calculated from $\zeta_{i+1}$, then
$
\zeta_{i+1} = \zeta_{i-1} + 2\Delta t (d\zeta - \nu\nabla^{2n}\zeta_{i+1})
$
As the application of $\nu\nabla^{2n}$ is, for every spectral mode, equivalent to a multiplication of a constant, we can
rewrite this to
$
\zeta_{i+1} = \frac{\zeta_{i-1} + 2\Delta t d\zeta}{1 + 2\Delta \nu\nabla^{2n}},
$
and expand the numerator to
$
\zeta_{i+1} = \zeta_{i-1} + 2\Delta t \frac{d\zeta - \nu\nabla^{2n}\zeta_{i-1}}{1+2\Delta t \nu \nabla^{2n}},
$
Hence the diffusion can be applied implicitly by updating the tendency $d\zeta$ as
$
d\zeta \to \frac{d\zeta - \nu\nabla^{2n}\zeta_{i-1}}{1+2\Delta t \nu \nabla^{2n}}
$
which only depends on $\zeta_{i-1}$. 

### Implementation details

```julia
using SpeedyWeather

P = Params(T=Float64)
G = GeoSpectral{P.T}(P)
B = Boundaries{P.T}(P,G)

fourier(B.ϕ0trunc,G),G)
```

## Time integration

SpeedyWeather.jl uses a leapfrog time scheme with a Robert's and William's filter
to dampen the computational mode and achieve 3rd order accuracy.

### Oscillation equation

$
\frac{dF}{dt} = i\omega F
$