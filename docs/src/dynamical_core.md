# Dynamical core

A mathematical and implementation-specific description of the dynamical core used in SpeedyWeather.jl.
We start by describing the barotropic vorticity equations which is one set of equations that
SpeedyWeather.jl can solve (see [How to run SpeedyWeather.jl](@ref)) as many details therein also
apply to the [Shallow water equations](@ref) and [Primitive equations](@ref) explained thereafter.

## Barotropic vorticity equation

The barotropic vorticity equation is the prognostic equation that describes the time evolution of
relative vorticity ``\zeta`` with advection, Coriolis force and diffusion in a single global layer.

```math
\frac{\partial \zeta}{\partial t} + \nabla \cdot (\mathbf{u}(\zeta + f)) = \nu\nabla^{2n}\zeta
```
with time ``t``, velocity vector ``\mathbf{u} = (u,v)``, Coriolis parameter ``f`` and hyperdiffusion
``\nu\nabla^{2n}\zeta`` (see [Horizontal diffusion](@ref)). Starting with some relative vorticity
``\zeta``, the [Laplacian](@ref) is inverted to obtain the stream function ``\Psi``

```math
\Psi = \nabla^{-2}\zeta
```

The zonal velocity ``u`` and meridional velocity ``v`` are then the (negative) meridional gradient
and zonal gradient of ``\Psi``

```math
\begin{aligned}
u &= -\frac{1}{R} \frac{\partial \Psi}{\partial \theta} \\
v &= \frac{1}{R\cos(\theta)} \frac{\partial \Psi}{\partial \phi} \\
\end{aligned}
```

which is described in [Gradients in spectral space](@ref).

### Algorithm

We briefly outline the algorithm that SpeedyWeather.jl uses in order to integrate the barotropic vorticity equation

1. Start with initial conditions of ``\zeta_{lm}`` in spectral space
2. Use ``\zeta_{lm}`` to
    - Invert the [Laplacian](@ref) to obtain the stream function ``\Psi_{lm}`` in spectral space
    - Transform ``\zeta_{lm}`` to ``\zeta`` in grid-point space
3. Use ``\Psi_lm`` to
    - obtain zonal velocity ``(\cos(\theta)u)_{lm}`` through a [Meridional derivative](@ref)
    - obtain meridional velocity ``(\cos(\theta)v)_{lm}`` through a [Zonal derivative](@ref)
4. Transform zonal and meridional velocity ``(\cos(\theta)u)_{lm}``, ``(\cos(\theta)v)_{lm}`` to grid-point space and unscale the ``\cos(\theta)`` factor to obtain ``u,v``.
5. Multiply ``u,v`` with ``\zeta+f`` in grid-point space
6. Transform ``u(\zeta + f)`` and ``v(\zeta+f)`` to spectral space
7. Compute the divergence of ``(\mathbf{u}(\zeta + f))_{lm}`` in spectral space through a [Meridional derivative](@ref) and [Zonal derivative](@ref) which will be the tendency of ``\zeta_{lm}``
8. Compute the [Horizontal diffusion](@ref) based on that tendency
9. Compute a leapfrog time step as described in [Time integration](@ref)
10. Repeat from 1.

## Shallow water equations

more to come

## Primitive equations

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

In SpeedyWeather.jl we use hyerdiffusion through an ``n``-th power Laplacian ``\nabla^{2n}`` (hyper when ``n>1``) which
can be implemented as a multiplication of the spectral coefficients ``\Psi_{lm}`` with ``(-l(l+1))^nR^{-2n}`` (see
spectral [Laplacian](@ref)) It is therefore computationally not more expensive to apply hyperdiffusion over diffusion
as the ``(-l(l+1))^nR^{-2n}`` can be precomputed.

In SpeedyWeather.jl the diffusion is applied _implicitly_. For that, consider a
[leapfrog scheme](https://en.wikipedia.org/wiki/Leapfrog_integration) with time step ``\Delta t`` of variable ``\zeta`` to
obtain from time steps ``i-1`` and ``i``, the next time step ``i+1``
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

## Radius scaling

SpeedyWeather.jl uses a scaling for vorticity ``\zeta`` and stream function ``\Psi`` that is
```math
\tilde{\zeta} = \zeta R, \tilde{\Psi} = \Psi R^{-1}.
```
In the barotropic voriticity equation model the inversion of the Laplcians in order to obtain
``\Psi`` from ``\zeta`` therefore becomes
```math
\tilde{\zeta} = \tilde{\nabla}^2 \tilde{\Psi}
```
where the dimensionless gradients simply omit the scaling with ``1/R``, ``\tilde{\nabla} = R\nabla``.
The scaled [Barotropic vorticity equation](@ref) is

```math
\partial_{\tilde{t}}\tilde{\zeta} + \tilde{\nabla} \cdot (\mathbf{u}(\tilde{\zeta} + \tilde{f})) = \tilde{\nu}\tilde{\nabla}^{2n}\tilde{\zeta}
```
with

- ``\tilde{t} = tR``, the scaled time ``t``
- ``\mathbf{u} = (u,v)``, the velocity vector (no scaling applied)
- ``\tilde{f} = fR``, the scaled Coriolis parameter ``f``
- ``\tilde{\nu} = \nu^* R``, the scaled viscosity ``\nu^*``, which itself is normalized to a damping time scale, see [Normalization of diffusion](@ref).

## Time integration

SpeedyWeather.jl uses a leapfrog time scheme with a Robert's and William's filter
to dampen the computational mode and achieve 3rd order accuracy.

### Oscillation equation

```math
\frac{dF}{dt} = i\omega F
```