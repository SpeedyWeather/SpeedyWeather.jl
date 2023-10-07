# Barotropic vorticity model

The barotropic vorticity model describes the evolution of a 2D non-divergent flow with
velocity components ``\mathbf{u} = (u,v)`` through self-advection, forces and dissipation.
Due to the non-divergent nature of the flow, it can be described by (the vertical component)
of the relative vorticity ``\zeta = \nabla \times \mathbf{u}``.

The dynamical core presented here to solve the barotropic vorticity equations largely follows
the idealized models with spectral dynamics developed at the
Geophysical Fluid Dynamics Laboratory[^1]: A barotropic vorticity model[^2].

Many concepts of the [Shallow water model](@ref) and the [Primitive equation model](@ref) are similar,
such that for example [horizontal diffusion](@ref diffusion) and the [Time integration](@ref leapfrog)
are only explained here.

## Barotropic vorticity equation

The barotropic vorticity equation is the prognostic equation that describes the time evolution of
relative vorticity ``\zeta`` with advection, Coriolis force, forcing and diffusion in a
single global layer on the sphere.

```math
\frac{\partial \zeta}{\partial t} + \nabla \cdot (\mathbf{u}(\zeta + f)) =
F_\zeta + \nabla \times \mathbf{F}_\mathbf{u} + (-1)^{n+1}\nu\nabla^{2n}\zeta
```
We denote time``t``, velocity vector ``\mathbf{u} = (u, v)``, Coriolis parameter ``f``,
and hyperdiffusion ``(-1)^{n+1} \nu \nabla^{2n} \zeta``
(``n`` is the hyperdiffusion order,  see [Horizontal diffusion](@ref diffusion)).
We also include possible forcing terms
``F_\zeta, \mathbf{F}_\mathbf{u} = (F_u,F_v)`` which act on the vorticity and/or on the
zonal velocity ``u`` and the meridional velocity ``v`` and hence the curl
``\nabla \times \mathbf{F}_\mathbf{u}`` is a tendency for relative vorticity ``\zeta``.
See [Extending SpeedyWeather](@ref) how to define these.

Starting with some relative vorticity ``\zeta``, the [Laplacian](@ref) is
inverted to obtain the stream function ``\Psi``

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

which is described in [Derivatives in spherical coordinates](@ref). Using ``u`` and ``v`` we can then
advect the absolute vorticity ``\zeta + f``. In order to avoid to calculate both the curl and the
divergence of a flux we rewrite the barotropic vorticity equation as
```math
\frac{\partial \zeta}{\partial t} = F_\zeta +
\nabla \times (\mathbf{F} + \mathbf{u}_\perp(\zeta + f)) + (-1)^{n+1}\nu\nabla^{2n}\zeta
```
with ``\mathbf{u}_\perp = (v,-u)`` the rotated velocity vector, because
``-\nabla\cdot\mathbf{u} = \nabla \times \mathbf{u}_\perp``. This is the form that is solved
in the `BarotropicModel`, as outlined in the following section.

## Algorithm

We briefly outline the algorithm that SpeedyWeather.jl uses in order to integrate the barotropic vorticity equation.
As an initial step

0\. Start with initial conditions of ``\zeta_{lm}`` in spectral space and
transform this model state to grid-point space:
- Invert the [Laplacian](@ref) of vorticity ``\zeta_{lm}`` to obtain the stream function ``\Psi_{lm}`` in spectral space
- obtain zonal velocity ``(\cos(\theta)u)_{lm}`` through a [Meridional derivative](@ref)
- obtain meridional velocity ``(\cos(\theta)v)_{lm}`` through a [Zonal derivative](@ref)
- Transform zonal and meridional velocity ``(\cos(\theta)u)_{lm}``, ``(\cos(\theta)v)_{lm}`` to grid-point space
- Unscale the ``\cos(\theta)`` factor to obtain ``u,v``
- Transform ``\zeta_{lm}`` to ``\zeta`` in grid-point space

Now loop over

1. Compute the forcing (or drag) terms ``F_\zeta, \mathbf{F}_\mathbf{u}``
2. Multiply ``u,v`` with ``\zeta+f`` in grid-point space
3. Add ``A = F_u + v(\zeta + f)`` and ``B = F_v - u(\zeta + f)``
4. Transform these vector components to spectral space ``A_{lm}``, ``B_{lm}``
5. Compute the curl of ``(A,B)_{lm}`` in spectral space, add to ``F_\zeta`` to accumulate the tendency of ``\zeta_{lm}``
6. Compute the [horizontal diffusion](@ref diffusion) based on that tendency
7. Compute a leapfrog time step as described in [Time integration](@ref leapfrog) with a [Robert-Asselin and Williams filter](@ref)
8. Transform the new spectral state of ``\zeta_{lm}`` to grid-point ``u,v,\zeta`` as described in 0.
9. Possibly do some output
10. Repeat from 1.



## [Horizontal diffusion](@id diffusion)

In SpeedyWeather.jl we use hyperdiffusion through an ``n``-th power Laplacian ``(-1)^{n+1}\nabla^{2n}`` (hyper when ``n>1``) which
can be implemented as a multiplication of the spectral coefficients ``\Psi_{lm}`` with ``(-l(l+1))^nR^{-2n}`` (see
spectral [Laplacian](@ref)) It is therefore computationally not more expensive to apply hyperdiffusion over diffusion
as the ``(-l(l+1))^nR^{-2n}`` can be precomputed. Note the sign change ``(-1)^{n+1}`` here is such that the dissipative nature of the diffusion operator is retained for ``n`` odd and even.

In SpeedyWeather.jl the diffusion is applied _implicitly_. For that, consider a
[leapfrog scheme](https://en.wikipedia.org/wiki/Leapfrog_integration) with time step ``\Delta t`` of variable ``\zeta`` to
obtain from time steps ``i-1`` and ``i``, the next time step ``i+1``
```math
\zeta_{i+1} = \zeta_{i-1} + 2\Delta t d\zeta,
```
with ``d\zeta`` being some tendency evaluated from ``\zeta_i``. Now we want to add a diffusion term ``(-1)^{n+1}\nu \nabla^{2n}\zeta``
with coefficient ``\nu``, which however, is implicitly calculated from ``\zeta_{i+1}``, then

```math
\zeta_{i+1} = \zeta_{i-1} + 2\Delta t (d\zeta + (-1)^{n+1} \nu\nabla^{2n}\zeta_{i+1})
```
As the application of ``(-1)^{n+1}\nu\nabla^{2n}`` is, for every spectral mode, equivalent to a multiplication of
a constant, we can rewrite this to
```math
\zeta_{i+1} = \frac{\zeta_{i-1} + 2\Delta t d\zeta}{1 - 2\Delta (-1)^{n+1}\nu\nabla^{2n}},
```
and expand the numerator to
```math
\zeta_{i+1} = \zeta_{i-1} + 2\Delta t \frac{d\zeta + (-1)^{n+1} \nu\nabla^{2n}\zeta_{i-1}}{1 - 2\Delta t (-1)^{n+1}\nu \nabla^{2n}},
```
Hence the diffusion can be applied implicitly by updating the tendency ``d\zeta`` as
```math
d\zeta \to \frac{d\zeta + (-1)^{n+1}\nu\nabla^{2n}\zeta_{i-1}}{1 - 2\Delta t \nu \nabla^{2n}}
```
which only depends on ``\zeta_{i-1}``. Now let ``D_\text{explicit} = (-1)^{n+1}\nu\nabla^{2n}`` be the explicit part and
``D_\text{implicit} = 1 - (-1)^{n+1} 2\Delta t \nu\nabla^{2n}`` the implicit part. Both parts can be precomputed and are
``D_\text{implicit} = 1 - 2\Delta t \nu\nabla^{2n}`` the implicit part. Both parts can be precomputed and are
only an element-wise multiplication in spectral space. For every spectral harmonic ``l,m`` we do
```math
d\zeta \to D_\text{implicit}^{-1}(d\zeta + D_\text{explicit}\zeta_{i-1}).
```
Hence 2 multiplications and 1 subtraction with precomputed constants.
However, we will normalize the (hyper-)Laplacians as described in the following. This also will take care of the alternating sign such that the diffusion operation is dissipative regardless the power ``n``.

### Normalization of diffusion

In physics, the Laplace operator ``\nabla^2`` is often used to represent diffusion due to viscosity in a fluid or
diffusion that needs to be added to retain numerical stability. In both cases,
the coefficient is ``\nu`` of units ``\text{m}^2\text{s}^{-1}`` and the full operator reads as ``\nu \nabla^2``
with units ``(\text{m}^2\text{s}^{-1})(\text{m}^{-2}) = \text{s}^{-1}``.
This motivates us to normalize the Laplace operator by a constant of units ``\text{m}^{-2}`` and the coefficient
by its inverse such that it becomes a damping timescale of unit ``\text{s}^{-1}``. Given the application in
spectral space we decide to normalize by the largest eigenvalue ``-l_\text{max}(l_\text{max}+1)`` such that
all entries in the discrete spectral Laplace operator are in ``[0,1]``. This also has the effect that the
alternating sign drops out, such that higher wavenumbers are always dampened and not amplified.
The normalized coefficient ``\nu^* = l_\text{max}(l_\text{max}+1)\nu`` (always positive) is
therefore reinterpreted as the (inverse) time scale at which the highest wavenumber is dampened
to zero due to diffusion. Together we have 
```math
D^\text{explicit}_{l,m} = -\nu^* \frac{l(l+1)}{l_\text{max}(l_\text{max}+1)}
```
and the hyper-Laplacian of power ``n`` follows as
```math
D^\text{explicit,n}_{l,m} = -\nu^* \left(\frac{l(l+1)}{l_\text{max}(l_\text{max}+1)}\right)^n
```
and the implicit part is accordingly ``D^\text{implicit,n}_{l,m} = 1 - 2\Delta t D^\text{explicit,n}_{l,m}``.
Note that the diffusion time scale ``\nu^*`` is then also scaled by the radius, see next section.

## [Radius scaling](@id scaling)

Similar to a non-dimensionalization of the equations, SpeedyWeather.jl scales the barotropic vorticity
equation with ``R^2`` to obtain normalized gradient operators as follows.
A scaling for vorticity ``\zeta`` and stream function ``\Psi`` is used that is
```math
\tilde{\zeta} = \zeta R, \tilde{\Psi} = \Psi R^{-1}.
```
This is also convenient as vorticity is often ``10^{-5}\text{ s}^{-1}`` in the atmosphere,
but the stream function more like ``10^5\text{ m}^2\text{ s}^{-1}`` and so this scaling
brings both closer to 1 with a typical radius of the Earth of 6371km.
The inversion of the Laplacians in order to obtain ``\Psi`` from ``\zeta`` therefore becomes
```math
\tilde{\zeta} = \tilde{\nabla}^2 \tilde{\Psi}
```
where the dimensionless gradients simply omit the scaling with ``1/R``, ``\tilde{\nabla} = R\nabla``.
The [Barotropic vorticity equation](@ref) scaled with ``R^2`` is
```math
\partial_{\tilde{t}}\tilde{\zeta} + \tilde{\nabla} \cdot (\mathbf{u}(\tilde{\zeta} + \tilde{f})) =
\nabla \times \tilde{\mathbf{F}} + (-1)^{n+1}\tilde{\nu}\tilde{\nabla}^{2n}\tilde{\zeta}
```
with
- ``\tilde{t} = tR^{-1}``, the scaled time ``t``
- ``\mathbf{u} = (u,v)``, the velocity vector (no scaling applied)
- ``\tilde{f} = fR``, the scaled Coriolis parameter ``f``
- ``\tilde{\mathbf{F}} = R\mathbf{F}``, the scaled forcing vector ``\mathbf{F}``
- ``\tilde{\nu} = \nu^* R``, the scaled diffusion coefficient ``\nu^*``, which itself is normalized to a damping time scale, see [Normalization of diffusion](@ref).

So scaling with the radius squared means we can use dimensionless operators, however, this comes at the
cost of needing to deal with both a time step in seconds as well as a scaled time step in seconds per
meter, which can be confusing. Furthermore, some constants like Coriolis or the diffusion coefficient
need to be scaled too during initialization, which may be confusing too because values are not
what users expect them to be. SpeedyWeather.jl follows the logic that the scaling to the prognostic
variables is only applied just before the time integration and variables are unscaled for output
and after the time integration finished. That way, the scaling is hidden as much as possible from
the user. In hopefully many other cases it is clearly denoted that a variable or constant is
*scaled*.

## [Time integration](@id leapfrog)

SpeedyWeather.jl is based on the [Leapfrog time integration](https://en.wikipedia.org/wiki/Leapfrog_integration]),
which, for relative vorticity ``\zeta``, is
in its simplest form
```math
\frac{\zeta_{i+1} - \zeta_{i-1}}{2\Delta t} = RHS(\zeta_i),
```
meaning we step from the previous time step ``i-1``, leapfrogging over the current time step``i``
to the next time step ``i+1`` by evaluating the tendencies on the right-hand side ``RHS``
at the current time step ``i``. The time stepping is done in spectral space.
Once the right-hand side ``RHS`` is evaluated, leapfrogging is a linear operation, meaning
that its simply applied to every spectral coefficient ``\zeta_{lm}`` as one would evaluate
it on every grid point in grid-point models.

For the Leapfrog time integration two time steps of the prognostic variables have to be stored,
``i-1`` and ``i``. Time step ``i`` is used to evaluate the tendencies which are then added
to ``i-1`` in a step that also swaps the indices for the next time step ``i \to i-1`` and ``i+1 \to i``,
so that no additional memory than two time steps have to be stored at the same time.

The Leapfrog time integration has to be initialized with an Euler forward step in order
to have a second time step ``i+1`` available when starting from ``i`` to actually leapfrog over.
SpeedyWeather.jl therefore does two initial time steps that are different from
the leapfrog time steps that follow and that have been described above.

1) an Euler forward step with ``\Delta t/2``, then
2) one leapfrog time step with ``\Delta t``, then
3) leapfrog with ``2 \Delta t`` till the end

This is particularly done in a way that after 2. we have ``t=0`` at ``i-1`` and ``t=\Delta t`` at ``i``
available so that 3. can start the leapfrogging without any offset from the intuitive spacing
``0,\Delta  t, 2\Delta t, 3\Delta t,...``. The following schematic can be useful

|                    | time at step ``i-1`` | time at step ``i`` | time step at ``i+1`` |
| ------------------ | -------------------- | ------------------ | -------------------- |
| Initial conditions | ``t = 0``            |                    |                      |
| 1: Euler           | (T) ``\quad t = 0``  |  ``t=\Delta t/2``  |                      |
| 2: Leapfrog with ``\Delta t``|``t = 0``|(T) ``\quad t = \Delta t/2``| ``t = \Delta t``|
| 3 to ``n``: Leapfrog with ``2\Delta t``|``t-\Delta t``|(T) ``\qquad \quad \quad t``| ``t+\Delta t`` |

The time step that is used to evaluate the tendencies is denoted with (T).
It is always the time step furthest in time that is available.


## Robert-Asselin and Williams filter

The standard leapfrog time integration is often combined with a Robert-Asselin filter[^Robert66][^Asselin72]
to dampen a computational mode. The idea is to start with a standard leapfrog step to obtain
the next time step ``i+1`` but then to correct the current time step ``i`` by applying a filter
which dampens the computational mode. The filter looks like a discrete Laplacian in time
with a ``(1, -2, 1)`` stencil, and so, maybe unsurprisingly, is efficient to filter out
a "grid-scale oscillation" in time, aka the computational mode. Let ``v`` be the unfiltered
variable and ``u`` be the filtered variable, ``F`` the right-hand side tendency,
then the standard leapfrog step is
```math
v_{i+1} = u_{i-1} + 2\Delta tF(v_i)
```
Meaning we start with a filtered variable ``u`` at the previous time step ``i-1``, evaluate
the tendency ``F(v_i)`` based on the current time step ``i`` to obtain an
unfiltered next time step ``v_{i+1}``. We then filter the current time step ``i``
(which will become ``i-1`` on the next iteration)
```math
u_i = v_i + \frac{\nu}{2}(v_{i+1} - 2v_i + u_{i-1})
```
by adding a discrete Laplacian with coefficient ``\tfrac{\nu}{2}`` to it, evaluated
from the available filtered and unfiltered time steps centred around ``i``:
``v_{i-1}`` is not available anymore because it was overwritten by the filtering
at the previous iteration, ``u_i, u_{i+1}`` are not filtered yet when applying
the Laplacian. The filter parameter ``\nu`` is typically chosen between 0.01-0.2,
with stronger filtering for higher values.

Williams[^Williams2009] then proposed an additional filter step to regain accuracy
that is otherwise lost with a strong Robert-Asselin filter[^Amezcua2011][^Williams2011].
Now let ``w`` be unfiltered, ``v`` be once filtered, and ``u`` twice filtered, then
```math
\begin{aligned}
w_{i+1} &= u_{i-1} + 2\Delta tF(v_i) \\
u_i &= v_i + \frac{\nu\alpha}{2}(w_{i+1} - 2v_i + u_{i-1}) \\
v_{i+1} &= w_{i+1} - \frac{\nu(1-\alpha)}{2}(w_{i+1} - 2v_i + u_{i-1})
\end{aligned}
```
with the Williams filter parameter ``\alpha \in [0.5,1]``. For ``\alpha=1``
we're back with the Robert-Asselin filter (the first two lines).

The Laplacian in the parentheses is often called a *displacement*,
meaning that the filtered value is displaced (or corrected) in the direction
of the two surrounding time steps. The Williams filter now also applies
the same displacement, but in the opposite direction to the next time
step ``i+1`` as a correction step (line 3 above) for a once-filtered
value ``v_{i+1}`` which will then be twice-filtered by the Robert-Asselin
filter on the next iteration. For more details see the referenced publications.

The initial Euler step (see [Time integration](@ref leapfrog), Table) is not filtered.
Both the the Robert-Asselin and Williams filter are then switched on for all
following leapfrog time steps.

## References

[^1]: Geophysical Fluid Dynamics Laboratory, [Idealized models with spectral dynamics](https://www.gfdl.noaa.gov/idealized-models-with-spectral-dynamics/)
[^2]: Geophysical Fluid Dynamics Laboratory, [The barotropic vorticity equation](https://www.gfdl.noaa.gov/wp-content/uploads/files/user_files/pjp/barotropic.pdf).
[^Robert66]: Robert, André. “The Integration of a Low Order Spectral Form of the Primitive Meteorological Equations.” Journal of the Meteorological Society of Japan 44 (1966): 237-245.
[^Asselin72]: ASSELIN, R., 1972: Frequency Filter for Time Integrations. Mon. Wea. Rev., 100, 487–490, doi:[10.1175/1520-0493(1972)100<0487:FFFTI>2.3.CO;2](https://doi.org/10.1175/1520-0493(1972)100<0487:FFFTI>2.3.CO;2.)
[^Williams2009]: Williams, P. D., 2009: A Proposed Modification to the Robert–Asselin Time Filter. Mon. Wea. Rev., 137, 2538–2546, [10.1175/2009MWR2724.1](https://doi.org/10.1175/2009MWR2724.1).
[^Amezcua2011]: Amezcua, J., E. Kalnay, and P. D. Williams, 2011: The Effects of the RAW Filter on the Climatology and Forecast Skill of the SPEEDY Model. Mon. Wea. Rev., 139, 608–619, doi:[10.1175/2010MWR3530.1](https://doi.org/10.1175/2010MWR3530.1). 
[^Williams2011]: Williams, P. D., 2011: The RAW Filter: An Improvement to the Robert–Asselin Filter in Semi-Implicit Integrations. Mon. Wea. Rev., 139, 1996–2007, doi:[10.1175/2010MWR3601.1](https://doi.org/10.1175/2010MWR3601.1). 
