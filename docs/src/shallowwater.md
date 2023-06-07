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
6. Compute the Bernoulli potential ``P`` by adding the kinetic energy ``\frac{1}{2}(u^2 + v^2)`` and the "geopotential" ``g\eta`` in grid-point space
7. Transform the Bernoulli potential to spectral space ``P_{lm}``, take the Laplacian and subtract from the divergence tendency
7. Compute the volume fluxes ``uh,vh`` in grid-point space via ``h = \eta + H - H_b``
7. Transform to spectral space and take the divergence for ``-\nabla \cdot (\mathbf{u}h)`` which is the tendency for ``\eta``
7. Add possibly forcing ``F_\eta`` for ``\eta`` in spectral space
7. Correct the tendencies following the [Semi implicit time integration](@ref) to prevent fast gravity waves from causing numerical instabilities
6. Compute the [horizontal diffusion](@ref diffusion) based on the ``\zeta,\mathcal{D}`` tendencies
7. Compute a leapfrog time step as described in [Time integration](@ref leapfrog)
8. Transform the new spectral state of ``\zeta_{lm}``, ``\mathcal{D}_{lm}``, ``\eta_{lm}`` to grid-point ``u,v,\zeta,\mathcal{D},\eta`` as described in 0.
9. Possibly do some output
10. Repeat from 1.

## Semi implicit time integration 

Probably the biggest advantage of a spectral model is its ability to solve (parts of) the equations implicitly
a low computational cost. The reason is that a linear operator can be easily inverted in spectral space,
removing the necessity to solve large equation systems. An operation like ``\Psi = \nabla^{-2}\zeta``
in grid-point space is costly because it requires a global communication, coupling all grid points.
In spectral space ``\nabla^2`` is a diagonal operator, meaning that there is no communication between
harmonics and its inversion is therefore easily done on a mode-by-mode basis of the harmonics.

This can be made use of when facing time stepping constraints with explicit schemes, such that
ridiculuously small time steps to resolve fast waves would otherwise result in a horribly slow simulation. 
In the shallow water system there are gravity waves that propagate at a wave speed of ``\sqrt{gH}``
(typically 300m/s), which, in order to not violate the
[CFL criterion](https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition)
for explicit time stepping, would need to be resolved.
Therefore, treating the terms that are responsible for gravity waves implicitly would remove
that time stepping constraint and allows us to run the simulation at the time step needed to
resolve the advective motion of the atmosphere, which is usually one or two orders of magnitude
longer than gravity waves.

In the following we will describe how the semi implicit time integration can be combined
with the [Leapfrog time stepping](@ref leapfrog) for a large increase in numerical
stability with gravity waves.



## Scaled shallow water equations

Similar to the scaled barotropic vorticity equations, the scaled shallow water equations scale the vorticity and the divergence equation with ``R^2``, but the continuity equation with ``R``

```math
\begin{aligned}
\frac{\partial \tilde{\zeta}}{\partial \tilde{t}} + \tilde{\nabla} \cdot (\mathbf{u}(\tilde{\zeta} + \tilde{f})) &=
\tilde{\nu}\tilde{\nabla}^{2n}\tilde{\zeta} \\
\frac{\partial \tilde{\mathcal{D}}}{\partial \tilde{t}} - \tilde{\nabla} \times (\mathbf{u}(\tilde{\zeta} + \tilde{f})) &=
-\tilde{\nabla}^2\left(\tfrac{1}{2}(u^2 + v^2) + g\eta \right) + \tilde{\nu}\tilde{\nabla}^{2n}\tilde{\mathcal{D}} \\
\frac{\partial \eta}{\partial \tilde{t}} + \tilde{\nabla} \cdot (\mathbf{u}h) &= 0.
\end{aligned}
```

## References

[^1]: Geophysical Fluid Dynamics Laboratory, [Idealized models with spectral dynamics](https://www.gfdl.noaa.gov/idealized-models-with-spectral-dynamics/)
[^2]: Geophysical Fluid Dynamics Laboratory, [The Shallow Water Equations](https://www.gfdl.noaa.gov/wp-content/uploads/files/user_files/pjp/shallow.pdf).