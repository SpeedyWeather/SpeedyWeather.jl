# Shallow water model

The shallow water model describes the evolution of a 2D flow described by its velocity and
an interface height that conceptually represents pressure.
A divergent flow affects the interface height which in turn can impose a pressure gradient
force onto the flow. The dynamics include advection, forces, dissipation, and continuity.

The following description of the shallow water model largely follows the idealized models with
spectral dynamics developed at the Geophysical Fluid Dynamics Laboratory[^1]:
The Shallow Water Equations[^2].

## Shallow water equations

The shallow water equations of velocity ``\mathbf{u} = (u,v)`` and interface height ``\eta``
(i.e. the deviation from the fluid's rest height ``H``) are,
formulated in terms of relative vorticity ``\zeta = \nabla \times \mathbf{u}``,
divergence ``\mathcal{D} = \nabla \cdot \mathbf{u}``

```math
\begin{aligned}
\frac{\partial \zeta}{\partial t} + \nabla \cdot (\mathbf{u}(\zeta + f)) &=
\nabla \times \mathbf{F} + (-1)^{n+1}\nu\nabla^{2n}\zeta, \\
\frac{\partial \mathcal{D}}{\partial t} - \nabla \times (\mathbf{u}(\zeta + f)) &=
\nabla \cdot \mathbf{F} -\nabla^2(\tfrac{1}{2}(u^2 + v^2) + g\eta) + (-1)^{n+1}\nu\nabla^{2n}\mathcal{D}, \\
\frac{\partial \eta}{\partial t} + \nabla \cdot (\mathbf{u}h) &= F_\eta,
\end{aligned}
```

We denote time``t``, Coriolis parameter ``f``, a forcing vector ``\mathbf{F} = (F_u,F_v)``,
hyperdiffusion ``(-1)^{n+1} \nu \nabla^{2n}`` (``n`` is the hyperdiffusion order,  see [Horizontal diffusion](@ref diffusion)),
gravitational acceleration ``g``, dynamic layer thickness ``h``, and a forcing for
the interface height ``F_\eta``. In the shallow water model the dynamics layer thickness ``h`` is
```math
h = \eta + H - H_b
```
that is, the layer thickness at rest ``H`` plus the interface height ``\eta`` minus orography ``H_b``.

In the shallow water system the flow can be described through ``u,v`` or ``\zeta,\mathcal{D}`` which
are related through the stream function ``\Psi`` and the velocity potential ``\Phi`` (which is zero
in the [Barotropic vorticity equation](@ref)).
```math
\begin{aligned}
\zeta &= \nabla^2 \Psi \\
\mathcal{D} &= \nabla^2 \Phi \\
\mathbf{u} &= \nabla^\perp \Psi + \nabla \Phi
\end{aligned}
```
With ``\nabla^\perp`` being the rotated gradient operator, in cartesian coordinates ``x,y``:
``\nabla^\perp = (-\partial_y, \partial_x)``. See [Derivatives in spherical coordinates](@ref)
for further details. Especially because the inversion of the [Laplacian](@ref) and the
gradients of ``\Psi, \Phi`` can be computed in a single pass, see
[U,V from vorticity and divergence](@ref).

The divergence/curl of the vorticity flux ``\mathbf{u}(\zeta + f)`` are combined with the
divergence/curl of the forcing vector ``\mathbf{F}``, as
```math
\begin{aligned}
- \nabla \cdot (\mathbf{u}(\zeta + f)) + \nabla \times \mathbf{F} &=
\nabla \times (\mathbf{F} + \mathbf{u}_\perp(\zeta + f)) \\
\nabla \times (\mathbf{u}(\zeta + f)) + \nabla \cdot \mathbf{F} &=
\nabla \cdot (\mathbf{F} + \mathbf{u}_\perp(\zeta + f))
\end{aligned}
```
equivalently to how this is done in the [Barotropic vorticity equation](@ref)
with ``\mathbf{u}_\perp = (v,-u)``.

## Algorithm

0\. Start with initial conditions of relative vorticity ``\zeta_{lm}``, divergence ``D_{lm}``,
and interface height ``\eta_{lm}`` in spectral space and transform this model state to grid-point space:
- Invert the [Laplacian](@ref) of ``\zeta_{lm}`` to obtain the stream function ``\Psi_{lm}`` in spectral space
- Invert the [Laplacian](@ref) of ``D_{lm}`` to obtain the velocity potential ``\Phi_{lm}`` in spectral space
- obtain velocities ``U_{lm} = (\cos(\theta)u)_{lm}, V_{lm} = (\cos(\theta)v)_{lm}`` from ``\nabla^\perp\Psi_{lm} + \nabla\Phi_{lm}``
- Transform velocities ``U_{lm}``, ``V_{lm}`` to grid-point space ``U,V``
- Unscale the ``\cos(\theta)`` factor to obtain ``u,v``
- Transform ``\zeta_{lm}``, ``D_{lm}``, ``\eta_{lm}`` to ``\zeta, D, \eta`` in grid-point space

Now loop over

1. Compute the forcing vector ``\mathbf{F} = (F_u,F_v)`` for ``u`` and ``v``
2. Multiply ``u,v`` with ``\zeta+f`` in grid-point space
3. Add ``A = F_u + v(\zeta + f)`` and ``B = F_v - u(\zeta + f)``
4. Transform these vector components to spectral space ``A_{lm}``, ``B_{lm}``
5. Compute the curl of ``(A,B)_{lm}`` in spectral space which is the tendency of ``\zeta_{lm}``
5. Compute the divergence of ``(A,B)_{lm}`` in spectral space which is the tendency of ``\mathcal{D}_{lm}``
6. Compute the kinetic energy ``\frac{1}{2}(u^2 + v^2)`` and transform to spectral space
6. Add to the kinetic energy the "geopotential" ``g\eta_{lm}`` in spectral space to obtain the Bernoulli potential
7. Take the Laplacian of the Bernoulli potential and subtract from the divergence tendency
7. Compute the volume fluxes ``uh,vh`` in grid-point space via ``h = \eta + H - H_b``
7. Transform to spectral space and take the divergence for ``-\nabla \cdot (\mathbf{u}h)`` which is the tendency for ``\eta``
7. Add possibly forcing ``F_\eta`` for ``\eta`` in spectral space
7. Correct the tendencies following the [semi-implicit time integration](@ref implicit_swm) to prevent fast gravity waves from causing numerical instabilities
6. Compute the [horizontal diffusion](@ref diffusion) based on the ``\zeta,\mathcal{D}`` tendencies
7. Compute a leapfrog time step as described in [Time integration](@ref leapfrog) with a [Robert-Asselin and Williams filter](@ref)
8. Transform the new spectral state of ``\zeta_{lm}``, ``\mathcal{D}_{lm}``, ``\eta_{lm}`` to grid-point ``u,v,\zeta,\mathcal{D},\eta`` as described in 0.
9. Possibly do some output
10. Repeat from 1.

## [Semi-implicit time integration](@id implicit_swm)

Probably the biggest advantage of a spectral model is its ability to solve (parts of) the equations implicitly
a low computational cost. The reason is that a linear operator can be easily inverted in spectral space,
removing the necessity to solve large equation systems. An operation like ``\Psi = \nabla^{-2}\zeta``
in grid-point space is costly because it requires a global communication, coupling all grid points.
In spectral space ``\nabla^2`` is a diagonal operator, meaning that there is no communication between
harmonics and its inversion is therefore easily done on a mode-by-mode basis of the harmonics.

This can be made use of when facing time stepping constraints with explicit schemes, where
ridiculously small time steps to resolve fast waves would otherwise result in a horribly slow simulation. 
In the shallow water system there are gravity waves that propagate at a wave speed of ``\sqrt{gH}``
(typically 300m/s), which, in order to not violate the
[CFL criterion](https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition)
for explicit time stepping, would need to be resolved.
Therefore, treating the terms that are responsible for gravity waves implicitly would remove
that time stepping constraint and allows us to run the simulation at the time step needed to
resolve the advective motion of the atmosphere, which is usually one or two orders of magnitude
longer than gravity waves.

In the following we will describe how the semi implicit time integration can be combined
with the [Leapfrog time stepping](@ref leapfrog) and the [Robert-Asselin and Williams filter](@ref)
for a large increase in numerical stability with gravity waves. Let ``V_i`` be the model
state of all prognostic variables at time step ``i``, the [leapfrog time stepping](@ref leapfrog)
is then

```math
\frac{V_{i+1} - V_{i-1}}{2\Delta t} = N(V_{i})
```
with the right-hand side operator ``N`` evaluated at the current time step ``i``. Now the
idea is to split the terms in ``N`` into non-linear terms that are evaluated explicitly
in ``N_E`` and into the linear terms ``N_I``, solved implicitly, that are responsible for
the gravity waves. Linearization happens around a state of rest without orography.

We could already assume to evaluate ``N_I`` at ``i+1``, but in fact,
we can introduce ``\alpha \in [0,1]`` so that for ``\alpha=0`` we use ``i-1`` (i.e. explicit),
for ``\alpha=1/2`` it is centred implicit ``\tfrac{1}{2}N_I(V_{i-1}) + \tfrac{1}{2}N_I(V_{i+1})``,
and for ``\alpha=1`` a fully backwards scheme ``N_I(V_{i+1})`` evaluated at ``i+1``.

```math
\frac{V_{i+1} - V_{i-1}}{2\Delta t} = N_E(V_{i}) + \alpha N_I(V_{i+1}) + (1-\alpha)N_I(V_{i-1})
```
Let ``\delta V = \tfrac{V_{i+1} - V_{i-1}}{2\Delta t}`` be the tendency we need for the Leapfrog
time stepping. Introducing ``\xi = 2\alpha\Delta t`` we have

```math
\delta V = N_E(V_i) + N_I(V_{i-1}) + \xi N_I(\delta V)
```
because ``N_I`` is a linear operator. This is done so that we can solve for ``\delta V``
by inverting ``N_I``, but let us gather the other terms as ``G`` first.
```math
G = N_E(V_i) + N_I(V_{i-1}) = N(V_i) + N_I(V_{i-1} - V_i)
```
For the shallow water equations we will only make use of the last formulation, meaning 
we first evaluate the whole right-hand side ``N(V_i)`` at the current time step
as we would do with fully explicit time stepping but then add the implicit terms
``N_I(V_{i-1} - V_i)`` afterwards to move those terms from ``i`` to ``i-1``.
Note that we could also directly evaluate the implicit terms at ``i-1``
as it is suggested in the previous formulation ``N_E(V_i) + N_I(V_{i-1})``, the result would be the same.
But in general it can be more efficient to do it one or the other way, and in fact
it is also possible to combine both ways. This will be discussed in the
[semi-implicit time stepping for the primitive equations](@ref implicit_primitive).

We can now implicitly solve for ``\delta V`` by
```math
\delta V = (1-\xi N_I)^{-1}G
```
So what is ``N_I``? In the shallow water system the gravity waves are caused by
```math
\begin{aligned}
\frac{\partial \mathcal{D}}{\partial t} &= -g\nabla^2\eta \\
\frac{\partial \eta}{\partial t} &= -H\mathcal{D},
\end{aligned}
```
which is a linearization of the equations around a state of rest with uniform constant layer
thickness ``h = H``. The continuity equation with the ``-\nabla(\mathbf{u}h)``
term, for example, is linearized to ``-\nabla(\mathbf{u}H) = -H\mathcal{D}``.
The divergence and continuity equations can now be written following the
``\delta V = G + \xi N_I(\delta V)`` formulation from above as
a coupled system (The vorticity equation is zero for the linear gravity wave equation in the shallow water equations,
hence no semi-implicit correction has to be made to the vorticity tendency).

```math
\begin{aligned}
\delta \mathcal{D} &= G_\mathcal{D} - \xi g \nabla^2 \delta \eta \\
\delta \eta &= G_\mathcal{\eta} - \xi H \delta\mathcal{D}
\end{aligned}
```
with

```math
\begin{aligned}
G_\mathcal{D} &= N_\mathcal{D} - \xi g \nabla^2 (\eta_{i-1} - \eta_i) \\
G_\mathcal{\eta} &= N_\eta - \xi H (\mathcal{D}_{i-1} - \mathcal{D}_i)
\end{aligned}
```
Inserting the second equation into the first, we can first solve for ``\delta \mathcal{D}``,
and then for ``\delta \eta``. Reminder that we do this in spectral space to every harmonic
independently, so the Laplace operator ``\nabla^2 = -l(l+1)`` takes the form of its eigenvalue
``-l(l+1)`` (normalized to unit sphere, as are the [scaled shallow water equations](@ref scaled_swm))
and its inversion is therefore just the inversion of this scalar.
```math
\delta D = \frac{G_\mathcal{D} - \xi g\nabla^2 G_\eta}{1 - \xi^2 H \nabla^2} =: S^{-1}(G_\mathcal{D} - \xi g\nabla^2 G_\eta) 
```
Where the last formulation just makes it clear that ``S = 1 - \xi^2 H \nabla^2`` is the
operator to be inverted. ``\delta \eta`` is then obtained via insertion as written above.
Equivalently, by adding a superscript ``l`` for every degree of the spherical harmonics,
we have
```math
\delta \mathcal{D}^l = \frac{G_\mathcal{D}^l + \xi g l(l+1) G_\eta^l}{1 + \xi^2 H l(l+1)}
```

The idea of the semi-implicit time stepping is now as follows:
1) Evaluate the right-hand side explicitly at time step ``i`` to obtain the explicit, preliminary tendencies ``N_\mathcal{D},N_\eta`` (and ``N_\zeta`` without a need for semi-implicit correction)
2) Move the implicit terms from ``i`` to ``i-1`` when calculating ``G_\mathcal{D}, G_\eta``
3) Solve for ``\delta \mathcal{D}``, the new, corrected tendency for divergence.
4) With ``\delta \mathcal{D}`` obtain ``\delta \eta``, the new, corrected tendency for ``\eta``.
5) Apply horizontal diffusion as a correction to ``N_\zeta, \delta \mathcal{D}`` as outlined in [Horizontal diffusion](@ref diffusion).
6) Leapfrog with tendencies that have been corrected for both semi-implicit and diffusion.

Some notes on the semi-implicit time stepping

- The inversion of the semi-implicit time stepping depends on ``\delta t``, that means every time the time step changes, the inversion has to be recalculated.
- You may choose ``\alpha = 1/2`` to dampen gravity waves but initialization shocks still usually kick off many gravity waves that propagate around the sphere for many days.
- With increasing ``\alpha > 1/2`` these waves are also slowed down, such that for ``\alpha = 1`` they quickly disappear in several hours.
- Using the [scaled shallow water equations](@ref scaled_swm) the time step ``\delta t`` has to be the scaled time step ``\tilde{\Delta t} = \delta t/R`` which is divided by the radius ``R``. Then we use the normalized eigenvalues ``-l(l+1)`` which also omit the ``1/R^2`` scaling, see [scaled shallow water equations](@ref scaled_swm) for more details.

## [Scaled shallow water equations](@id scaled_swm)

Similar to the [scaled barotropic vorticity equations](@ref scaling),
SpeedyWeather.jl scales in the shallow water equations.
The vorticity and the divergence equation are scaled with ``R^2``, the radius of the sphere squared,
but the continuity equation is scaled with ``R``. We also combine the vorticity flux and forcing into
a single divergence/curl operation as mentioned in [Shallow water equations](@ref) above

```math
\begin{aligned}
\frac{\partial \tilde{\zeta}}{\partial \tilde{t}} &=
\tilde{\nabla} \times (\tilde{\mathbf{F}} + \mathbf{u}_\perp(\tilde{\zeta} + \tilde{f})) +
(-1)^{n+1}\tilde{\nu}\tilde{\nabla}^{2n}\tilde{\zeta} \\
\frac{\partial \tilde{\mathcal{D}}}{\partial \tilde{t}} &=
\tilde{\nabla} \cdot (\tilde{\mathbf{F}} + \mathbf{u}_\perp(\tilde{\zeta} + \tilde{f})) -
\tilde{\nabla}^2\left(\tfrac{1}{2}(u^2 + v^2) + g\eta \right) +
(-1)^{n+1}\tilde{\nu}\tilde{\nabla}^{2n}\tilde{\mathcal{D}} \\
\frac{\partial \eta}{\partial \tilde{t}} &=
- \tilde{\nabla} \cdot (\mathbf{u}h) + \tilde{F}_\eta.
\end{aligned}
```

As in the [scaled barotropic vorticity equations](@ref scaling), one needs to scale
the time step, the Coriolis force, the forcing and the diffusion coefficient, but then
enjoys the luxury of working with dimensionless gradient operators. As before,
SpeedyWeather.jl will scale vorticity and divergence just before the model integration
starts and unscale them upon completion and for output. In the 
[semi-implicit time integration](@ref implicit_swm) we solve an equation that also
has to be scaled. It is with radius squared scaling (because it is the tendency for
the divergence equation which is also scaled with ``R^2``)

```math
R^2 \delta D = R^2\frac{G_\mathcal{D} - \xi g\nabla^2 G_\eta}{1 - \xi^2 H \nabla^2}
```
As ``G_\eta`` is only scaled with ``R`` we have
```math
\tilde{\delta D} = \frac{\tilde{G_\mathcal{D}} - \tilde{\xi} g\tilde{\nabla}^2 \tilde{G_\eta}}{1 - \tilde{\xi}^2 H \tilde{\nabla}^2}
```
The ``R^2`` normalizes the Laplace operator in the numerator, but using the scaled ``G_\eta`` we also scale ``\xi``
(which is convenient, because the time step within is the one we use anyway). The denominator ``S``
does not actually change because ``\xi^2\nabla^2 = \tilde{\xi}^2\tilde{\nabla}^2`` as ``\xi^2`` is scaled with
``1/R^2``, but the Laplace operator with ``R^2``. So overall we just have to use the scaled time step ``\tilde{\Delta t}``
and normalized eigenvalues for ``\tilde{\nabla}^2``.


## References

[^1]: Geophysical Fluid Dynamics Laboratory, [Idealized models with spectral dynamics](https://www.gfdl.noaa.gov/idealized-models-with-spectral-dynamics/)
[^2]: Geophysical Fluid Dynamics Laboratory, [The Shallow Water Equations](https://www.gfdl.noaa.gov/wp-content/uploads/files/user_files/pjp/shallow.pdf).