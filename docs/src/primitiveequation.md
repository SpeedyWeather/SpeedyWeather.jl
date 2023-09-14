# Primitive equation model

The [primitive equations](https://en.wikipedia.org/wiki/Primitive_equations) are a hydrostatic approximation
of the compressible Navier-Stokes equations for an ideal gas on a rotating sphere. We largely follow
the idealised spectral dynamical core developed by GFDL[^GFDL1] and documented therein[^GFDL2].

The primitive equations solved by SpeedyWeather.jl for relative vorticity ``\zeta``, divergence ``\mathcal{D}``,
logarithm of surface pressure ``\ln p_s``, temperature ``T`` and specific humidity ``q`` are

```math
\begin{aligned}
\frac{\partial \zeta}{\partial t} &= \nabla \times (\mathbf{\mathcal{P}}_\mathbf{u}
+ (f+\zeta)\mathbf{u}_\perp - W(\mathbf{u}) - R_dT_v\nabla \ln p_s) \\
\frac{\partial \mathcal{D}}{\partial t} &= \nabla \cdot (\mathcal{P}_\mathbf{u}
+ (f+\zeta)\mathbf{u}_\perp - W(\mathbf{u}) - R_dT_v\nabla \ln p_s) - \nabla^2(\frac{1}{2}(u^2 + v^2) + \Phi) \\
\frac{\partial \ln p_s}{\partial t} &= -\frac{1}{p_s} \nabla \cdot \int_0^{p_s} \mathbf{u}~dp \\
\frac{\partial T}{\partial t} &= \mathcal{P}_T -\nabla\cdot(\mathbf{u}T) + T\mathcal{D} - W(T) + \kappa T_v \frac{D \ln p}{Dt} \\
\frac{\partial q}{\partial t} &= \mathcal{P}_q -\nabla\cdot(\mathbf{u}q) + q\mathcal{D} - W(q)\\
\end{aligned}
```

with velocity ``\mathbf{u} = (u,v)``, rotated velocity ``\mathbf{u}_\perp = (v,-u)``,
Coriolis parameter ``f``, ``W`` the [Vertical advection](@ref) operator, dry air gas constant ``R_d``,
[Virtual temperature](@ref) ``T_v``, [Geopotential](@ref) ``\Phi``, pressure ``p``
and surface pressure ``p_s``, thermodynamic ``\kappa = R\_d/c_p``
with ``c_p`` the heat capacity at constant pressure. Horizontal hyper diffusion of the
form ``(-1)^{n+1}\nu\nabla^{2n}`` with coefficient ``\nu`` and power ``n``  is added for
every variable that is advected, meaning ``\zeta, \mathcal{D}, T, q``, but left out
here for clarity, see [Horizontal diffusion](@ref diffusion).

The parameterizations for the tendencies of ``u,v,T,q`` from physical processes are denoted as
``\mathcal{P}_\mathbf{u} = (\mathcal{P}_u, \mathcal{P}_v), \mathcal{P}_T, \mathcal{P}_q``
and are further described in the corresponding sections, see [Parameterizations](@ref parameterizations).

SpeedyWeather.jl implements a `PrimitiveWet` and a `PrimitiveDry` dynamical core.
For a dry atmosphere, we have ``q = 0`` and the virtual temperature ``T_v = T``
equals the temperature (often called _absolute_ to distinguish from the virtual temperature).
The terms in the primitive equations and their discretizations are discussed in
the following sections. 

## Virtual temperature

!!! info "In short: Virtual temperature"
    Virtual temperature is the temperature dry air would need to have to be
    as light as moist air. It is used in the dynamical core to include the
    effect of humidity on the density while replacing density through the
    ideal gas law with temperature.

We assume the atmosphere to be composed of two ideal gases: Dry air and water vapour.
Given a specific humidity ``q`` both gases mix, their pressures ``p_d``, ``p_w``
(``d`` for dry, ``w`` for water vapour), and densities ``\rho_d, \rho_w`` add in a given
air parcel that has temperature ``T``. The ideal gas law then holds for both gases
```math
\begin{aligned}
p_d &= \rho_d R_d T \\
p_w &= \rho_w R_w T \\
\end{aligned}
```
with the respective specific gas constants ``R_d = R/m_d`` and ``R_w = R/m_w``
obtained from the universal gas constant ``R`` divided by the molecular masses
of the gas. The total pressure ``p`` in the air parcel is
```math
p = p_d + p_w = (\rho_d R_d + \rho_w R_w)T
```
We ultimately want to replace the density ``\rho = \rho_w + \rho_d`` in the dynamical core,
using the ideal gas law, with the temperature ``T``, so that we never have
to calculate the density explicitly. However, in order to not deal with
two densities (dry air and water vapour) we would like to replace
temperature with a virtual temperature that includes the effect of
humidity on the density. So, wherever we use the ideal gas law
to replace density with temperature, we would use the virtual temperature,
which is a function of the absolute temperature and specific humidity,
instead. A higher specific humidity in an air parcel lowers 
the density as water vapour is lighter than dry air. Consequently,
the virtual temperature of moist air is higher than its absolute temperature
because warmer air is lighter too at constant pressure. We therefore
think of the virtual temperature as the temperature dry air would need to have
to be as light as moist air.

Starting with the last equation, with some manipulation we can write
the ideal gas law as total density ``\rho`` times a gas constant
times the virtual temperature that is supposed to be a function
of absolute temperature, humidity and some constants
```math
p  = (\rho R_d + \rho_w (R_w - R_d)) T = \rho R_d (1 +
\frac{1 - \tfrac{R_d}{R_w}}{\tfrac{R_d}{R_w}} \frac{\rho_w}{\rho_w + \rho_d})T
```
Now we identify
```math
\mu = \frac{1 - \tfrac{R_d}{R_w}}{\tfrac{R_d}{R_w}}
```
as some constant that is positive for water vapour being lighter than dry air
(``\tfrac{R_d}{R_w} = \tfrac{m_w}{m_d} < 1``) and
```math
q = \frac{\rho_w}{\rho_w + \rho_d}
```
as the specific humidity. Given temperature ``T`` and specific humidity ``q``, we
can therefore calculate the virtual temperature ``T_v`` as
```math
T_v = (1 + \mu q)T
```

For completeness we want to mention here that the above product, because it is a 
product of two variables ``q,T`` has to be computed in grid-point space,
see [Spherical Harmonic Transform](@ref).
To obtain an approximation to the virtual temperature in spectral space without
expensive transforms one can linearize
```math
T_v \approx T + \mu q\bar{T}
```
with a global constant temperature ``\bar{T}``, for example obtained from the
``l=m=0`` mode, ``\bar{T} = T_{0,0}\frac{1}{\sqrt{4\pi}}`` but depending on the
normalization of the spherical harmonics that factor needs adjustment.
We call this the _linear virtual temperature_ which is used for the geopotential
calculation, see [#254](https://github.com/SpeedyWeather/SpeedyWeather.jl/issues/254).

## Vertical coordinates

We start with some general considerations that apply when changing the vertical
coordinate from height ``z`` to something else. Let ``\Psi(x,y,z,t)``
be some variable that depends on space and time. Now we want to express
``\Psi`` using some other coordinate ``\eta`` in the vertical. Regardless of
the coordinate system the value of ``\Psi`` at the to ``z`` corresponding ``\eta``
(and vice versa) has to be the same as we only want to change the coordinate,
not ``\Psi`` itself.

```math
\Psi(x,y,\eta,t) = \Psi(x,y,z(x,y,\eta,t),t)
```
So you can think of ``z`` as a function of ``\eta`` and ``\eta`` as a function of ``z``.
The chain rule lets us differentiate ``\Psi`` with respect to ``z`` or ``\eta``
```math
\frac{\partial \Psi}{\partial z} = \frac{\partial \Psi}{\partial \eta}\frac{\partial \eta}{\partial z},
\qquad \frac{\partial \Psi}{\partial \eta} = \frac{\partial \Psi}{\partial z}\frac{\partial z}{\partial \eta}
```
But for derivatives with respect to ``x,y,t`` we have to apply the multi-variable
chain-rule as both ``\Psi`` and ``\eta`` depend on it. So a derivative with respect to
``x`` on ``\eta`` levels (where ``\eta`` constant) becomes
```math
\left. \frac{\partial \Psi}{\partial x}\right\vert_\eta = 
\left. \frac{\partial \Psi}{\partial x}\right\vert_z +
\frac{\partial \Psi}{\partial z}
\left. \frac{\partial z}{\partial x}\right\vert_\eta
```
So we first take the derivative of ``\Psi`` with respect to ``x``, but then also have to
account for the fact that, at a given ``\eta``, ``z`` depends on ``x`` which is again
dealt with using the univariate chain rule from above. We will make use of that
for the [Pressure gradient](@ref).

### Sigma coordinates

The problem with pure pressure coordinates is that they are not terrain-following.
For example, the 1000 hPa level in the Earth's atmosphere cuts through mountains.
A flow field on such a level is therefore not continuous and one would need to deal with
boundaries. Especially with spherical harmonics we need a terrain-following vertical
coordinate to transform between continuous fields in grid-point space and spectral space.

SpeedyWeather.jl currently uses so-called sigma coordinates for the vertical. 
This coordinate system uses fraction of surface pressure in the vertical, i.e.
```math
\sigma = \frac{p}{p_s}
```
with ``\sigma = [0,1]`` and ``\sigma = 0`` being the top (zero pressure) and ``\sigma = 1``
the surface (at surface pressure). As a consequence the vertical dimension is also
indexed from top to surface.

!!! info "Vertical indexing"
    Pressure, sigma, or hybrid coordinates in the vertical range from lowest values at the top
    to highest values at the surface. Consistently, we also index the vertical dimension top to
    surface. This means that ``k=1`` is the top-most layer, and ``k=N_{lev}`` (or similar)
    is the layer that sits directly above the surface.

Sigma coordinates are therefore terrain-following, as ``\sigma = 1`` is always at surface pressure
and so this level bends itself around every mountain, although the actual pressure on this
level can vary. For a visualisation see [#329](https://github.com/SpeedyWeather/SpeedyWeather.jl/issues/329).

One chooses ``\sigma`` levels associated with the ``k``-th layer and the pressure
can be reobtained from the surface pressure ``p_s``
```math
p_k = \sigma_kp_s
```
The layer thickness in terms of pressure is
```math
\Delta p_k = p_{k+\tfrac{1}{2}} - p_{k-\tfrac{1}{2}} =
(\sigma_{k+\tfrac{1}{2}} - \sigma_{k-\tfrac{1}{2}}) p_s = \Delta \sigma_k p_s
``` 
which can also be expressed with the layer thickness in sigma coordinates ``\Delta \sigma_k``
times the surface pressure. In SpeedyWeather.jl one chooses the half levels
``\sigma_{k+\tfrac{1}{2}}`` first and then obtains the full levels through averaging
```math
\sigma_k = \frac{\sigma_{k+\tfrac{1}{2}} + \sigma_{k-\tfrac{1}{2}}}{2}
```

## Geopotential

In the hydrostatic approximation the vertical momentum equation becomes

```math
\frac{\partial p}{\partial z} = -\rho g,
```
meaning that the (negative) vertical pressure gradient is given by the density in 
that layer times the gravitational acceleration. The heavier the fluid the more
the pressure will increase below. Inserting the ideal gas
law
```math
\frac{\partial gz}{\partial p} = -\frac{R_dT_v}{p},
```
with the geopotential ``\Phi = gz`` we can write this in terms of the logarithm
of pressure
```math
\frac{\partial \Phi}{\partial \ln p} = -R_dT_v.
```
Note that we use the [Virtual temperature](@ref) here as we replaced the density
through the ideal gas law with temperature. Given a vertical temperature profile
``T_v`` and the (constant) surface geopotential ``\Phi_s = gz_s`` where ``z_s``
is the orography, we can integrate this equation from the surface to the top
to obtain ``\Phi_k`` on every layer ``k``.
The surface is at ``k = N+\tfrac{1}{2}`` (see [Vertical coordinates](@ref)) 
with ``N`` vertical levels. We can integrate the geopotential onto half levels as
(``T_k^v`` is the virtual temperature at layer ``k``, the subscript ``v`` has been
moved to be a superscript)
```math
\Phi_{k-\tfrac{1}{2}} = \Phi_{k+\tfrac{1}{2}} + R_dT^v_k(\ln p_{k+1/2} - \ln p_{k-1/2})
```
or onto full levels with
```math
\Phi_{k} = \Phi_{k+\tfrac{1}{2}} + R_dT^v_k(\ln p_{k+1/2} - \ln p_k).
```
We use this last formula first to get from ``\Phi_s`` to ``\Phi_N``, and then for
every ``k`` twice to get from ``\Phi_k`` to ``\Phi_{k-1}`` via ``\Phi_{k-\tfrac{1}{2}}``.
For the first half-level integration we use ``T_k`` for the second ``T_{k-1}``.

!!! warning "Semi-implicit time integration: Geopotential"
    With the semi-implicit time integration in SpeedyWeather the
    Geopotential is not calculated from the spectral temperature 
    at the current, but at the previous time step.
    This is because this is a linear term that we solve implicitly
    to avoid instabilities from gravity waves.
    For details see section [Semi-implicit time stepping](@ref implicit_primitive).

## Surface pressure tendency

The surface pressure increases with a convergence of the flow above. Written in terms
of the surface pressure directly, and not its logarithm
```math
\frac{\partial p_s}{\partial t} = -\nabla \cdot \int_0^{p_s} \mathbf{u}~dp
```
For ``k`` discrete layers from 1
at the top to ``N`` at the surface layer this can be written as
```math
\frac{\partial p_s}{\partial t} = - \sum_{k=1}^N \nabla \cdot (\mathbf{u}_k \Delta p_k)
```
which can be thought of as a vertical integration of the pressure thickness-weighted divergence.
In ``\sigma``-coordinates with ``\Delta p_k = \Delta \sigma_k p_s`` (see [Vertical coordinates](@ref))
this becomes
```math
\frac{\partial p_s}{\partial t} = - \sum_{k=1}^N \sigma_k \nabla \cdot (\mathbf{u}_k p_s)
= -\sum_{k=1}^N \sigma_k (\mathbf{u}_k \cdot \nabla p_s + p_s \nabla \cdot \mathbf{u}_k)
```
Using the logarithm of pressure ``\ln p`` as the vertical coordinate this becomes
```math
\frac{\partial \ln p_s}{\partial t} = 
-\sum_{k=1}^N \sigma_k (\mathbf{u}_k \cdot \nabla \ln p_s + \nabla \cdot \mathbf{u}_k)
```
The second term is the divergence ``\mathcal{D}_k`` at layer ``k``.
We introduce ``\bar{a} = \sum_k \Delta \sigma_k a_k``, the ``\sigma``-weighted vertical integration operator
applied to some variable ``a``. This is essentially an average as ``\sum_k \Delta \sigma_k = 1``.
The surface pressure tendency can then be written as
```math
\frac{\partial \ln p_s}{\partial t} = 
-\mathbf{\bar{u}} \cdot \nabla \ln p_s - \bar{\mathcal{D}}
```
which is form used by SpeedyWeather.jl to calculate the tendency of (the logarithm of) surface pressure.

As we will have ``\ln p_s`` available in spectral space at the beginning of a time step, the
gradient can be easily computed (see [Derivatives in spherical coordinates](@ref)). However,
we then need to transform both gradients to grid-point space for the scalar product with 
the (vertically ``\sigma``-averaged) velocity vector ``\mathbf{\bar{u}}`` before transforming it
back to spectral space where the tendency is needed. In general, we can do the ``\sigma``-weighted
average in spectral or in grid-point space, although it is computationally cheaper in spectral space.
We therefore compute ``- \bar{\mathcal{D}}`` entirely in spectral space. With ``()`` denoting
spectral space and ``[]`` grid-point space (hence, ``([])`` and ``[()]`` are the transforms in the
respective directions) we therefore do
```math
\left(\frac{\partial \ln p_s}{\partial t}\right) = 
\left(-\mathbf{\overline{[u]}} \cdot [\nabla (\ln p_s)]\right) - \overline{(\mathcal{D})}
```
But note that it would also be possible to do
```math
\left(\frac{\partial \ln p_s}{\partial t}\right) = 
\left(-\mathbf{\overline{[u]}} \cdot [\nabla (\ln p_s)] - \overline{[\mathcal{D}]}\right)
```
Meaning that we would compute the vertical average in grid-point space, subtract from the
pressure gradient flux before transforming to spectral space. The same amount of transforms
are performed but in the latter, the vertical averaging is done in grid-point space.

!!! warning "Semi-implicit time integration: Surface pressure tendency"
    With the semi-implicit time integration in SpeedyWeather the
    ``- \overline{(\mathcal{D})}`` term is not evaluated from the spectral divergence
    ``\mathcal{D}`` at the current, but at the previous time step.
    This is because this is a linear term that we solve implicitly
    to avoid instabilities from gravity waves.
    For details see section [Semi-implicit time stepping](@ref implicit_primitive).

## Vertical advection

The advection equation ``\tfrac{DT}{Dt} = 0`` for a tracer ``T`` is, in flux form,
for layer ``k``
```math
\frac{\partial (T_k \Delta p_k)}{\partial t} = - \nabla \cdot (\mathbf{u}_k T_k \Delta p_k)
- (M_{k+\tfrac{1}{2}}T_{k+\tfrac{1}{2}} - M_{k-\tfrac{1}{2}}T_{k-\tfrac{1}{2}})
```
which can be through the gradient product rule, and using the conservation of mass
(see [Vertical velocity](@ref)) transformed into an advective form. In sigma coordinates this simplifies to
```math
\frac{\partial T_k}{\partial t} = - \mathbf{u}_k \cdot \nabla T_k
- \frac{1}{\Delta \sigma_k}\left(\dot{\sigma}_{k+\tfrac{1}{2}}(T_{k+\tfrac{1}{2}} - T_k) - \dot{\sigma}_{k-\tfrac{1}{2}}(T_k - T_{k-\tfrac{1}{2}})\right)
```
With the reconstruction at the faces, ``T_{k+\tfrac{1}{2}}``, and ``T_{k-\tfrac{1}{2}}`` depending on one's choice
of the advection scheme. For a second order centered scheme, we choose ``T_{k+\tfrac{1}{2}} = \tfrac{1}{2}(T_k + T_{k+1})``
and obtain
```math
\frac{\partial T_k}{\partial t} = - \mathbf{u}_k \cdot \nabla T_k
- \frac{1}{2\Delta \sigma_k}\left(\dot{\sigma}_{k+\tfrac{1}{2}}(T_{k+1} - T_k) + \dot{\sigma}_{k-\tfrac{1}{2}}(T_k - T_{k-1})\right)
```
However, note that this scheme is dispersive and easily leads to instabilities at higher resolution, where a
more advanced vertical advection scheme becomes necessary. For convenience, we may write ``W(T)``
to denote the vertical advection term ``\dot{\sigma}\partial_\sigma T``, without specifying which schemes is used.
The vertical velocity ``\dot{\sigma}`` is calculated as described in the following.

### Vertical velocity

In the section [Surface pressure tendency](@ref) we used that the surface pressure changes
with the convergence of the flow above, which derives from the conservation of mass.
Similarly, the conservation of mass for layer ``k`` can be expressed as
(setting ``T=1`` in the advection equation in section [Vertical advection](@ref))
```math
\frac{\partial \Delta p_k}{\partial t} = -\nabla \cdot (\mathbf{u}_k \Delta p_k)
- (M_{k+\tfrac{1}{2}} - M_{k-\tfrac{1}{2}})
```
Meaning that the pressure thickness ``\Delta p_k`` of layer ``k`` changes with
a horizontal divergence ``-\nabla \cdot (\mathbf{u}_k \Delta p_k)`` if not
balanced by a net vertical mass flux ``M`` into of the layer through the
bottom and top boundaries of ``k`` at ``k\pm\tfrac{1}{2}``. ``M`` is defined positive
downward as this is the direction in which both pressure and sigma coordinates increase.
The boundary conditions are ``M_\tfrac{1}{2} = M_{N+\tfrac{1}{2}} = 0``, such that there
is no mass flux into the top layer from above or out of the surface layer ``N`` and into the ground
or ocean.

When integrating from the top down to layer ``k`` we obtain the mass flux downwards out of layer ``k``
```math
M_{k+\tfrac{1}{2}} = - \sum_{r=1}^k \nabla \cdot (\mathbf{u}_k \Delta p_k) - \frac{\partial p_{k+\tfrac{1}{2}}}{\partial t}
```
In sigma coordinates we have ``M_{k+\tfrac{1}{2}} = p_s \dot{\sigma}_{k+\tfrac{1}{2}}`` with
``\dot{\sigma}`` being the vertical velocity in sigma coordinates, also defined at interfaces
between layers. To calculate ``\dot{\sigma}`` we therefore compute
```math
\dot{\sigma}_{k+\tfrac{1}{2}} = \frac{M_{k+\tfrac{1}{2}}}{p_s} = 
- \sum_{r=1}^k \Delta \sigma_r (\mathbf{u}_k \cdot \nabla \ln p_s + \mathcal{D}_r) 
+ \sigma_{k+\tfrac{1}{2}}(-\mathbf{\bar{u}} \cdot \nabla \ln p_s - \bar{\mathcal{D}})
```
With ``\bar{A}`` denoting a sigma thickness-weighted vertical average as in section [Surface pressure tendency](@ref).
Now let ``\bar{A_k}`` be that average from ``r=1`` to ``r=k`` only and not necessarily down to the surface, as required in the
equation above, then we can also write
```math
\dot{\sigma}_{k+\tfrac{1}{2}} = 
- \overline{\mathbf{u}_k \cdot \nabla \ln p_s} - \bar{\mathcal{D}}_k
+ \sigma_{k+\tfrac{1}{2}}(-\mathbf{\bar{u}} \cdot \nabla \ln p_s - \bar{\mathcal{D}})
```
See also Hoskins and Simmons, 1975[^HS75]. These vertical averages are the same as required by the 
[Surface pressure tendency](@ref) and in the [Temperature equation](@ref), they are therefore all calculated
at once, storing the partial averages ``\overline{\mathbf{u}_k \cdot \nabla \ln p_s}`` and ``\bar{\mathcal{D}}_k`` on the fly.

## Pressure gradient

The pressure gradient term in the primitive equations is
```math
-\frac{1}{\rho}\nabla_z p
```
with density ``\rho`` and pressure ``p``. The gradient here is taken at constant ``z`` hence the
subscript. If we move to a pressure-based vertical coordinate system we will need to evaluate
gradients on constant levels of pressure though, i.e. ``\nabla_p``. There is, by definition,
no gradient of pressure on constant levels of pressure, but we can use the chain rule (see 
[Vertical coordinates](@ref)) to rewrite this as (use only ``x`` but ``y`` is equivalent)
```math
0 = \left. \frac{\partial p}{\partial x} \right\vert_p =
\left. \frac{\partial p}{\partial x} \right\vert_z +
\frac{\partial p}{\partial z}\left. \frac{\partial z}{\partial x} \right\vert_p
```
Using the hydrostatic equation ``\partial_z p = -\rho g`` this becomes
```math
\left. \frac{\partial p}{\partial x} \right\vert_z = \rho g \left. \frac{\partial z}{\partial x} \right\vert_p
```
Or, in terms of the geopotential ``\Phi = gz``
```math
\frac{1}{\rho}\nabla_z p = \nabla_p \Phi
```
which is the actual reason why we use pressure coordinates: As density ``\rho`` also depends on
the pressure ``p`` the left-hand side means an implicit system when solving for pressure ``p``.
To go from pressure to sigma coordinates we apply the chain rule from section
[Vertical coordinates](@ref) again and obtain
```math
\nabla_p \Phi = \nabla_\sigma \Phi - \frac{\partial \Phi}{\partial p}\nabla_\sigma p
= \nabla_\sigma \Phi + \frac{1}{\rho}\nabla_\sigma p
```
where the last step inserts the hydrostatic equation again. With the ideal gas law, and note
that we use [Virtual temperature](@ref) ``T_v`` everywhere where the ideal gas law is used,
but in combination with the dry gas constant ``R_d``
```math
\nabla_p \Phi = \nabla_\sigma \Phi + \frac{R_dT_v}{p} \nabla_\sigma p
```
Combining the pressure in denominator and gradient to the logarithm and with
``\nabla \ln p = \nabla \ln p_s`` in [Sigma coordinates](@ref) (the logarithm of
``\sigma_k`` adds a constant that drops out in the gradient) we therefore
have
```math
- \frac{1}{\rho}\nabla_z p = -\nabla_p \Phi = -\nabla_\sigma \Phi - R_dT_v \nabla_\sigma \ln p_s
```
From left to right: The pressure gradient force in ``z``-coordinates; in pressure coordinates;
and in sigma coordinates. Each denoted with the respective subscript on gradients. 
SpeedyWeather.jl uses the latter.
In sigma coordinates we may drop the ``\sigma`` subscript on gradients, but still meaning
that the gradient is evaluated on a surface of our vertical coordinate.
In vorticity-divergence formulation of the momentum equations the ``\nabla_\sigma \Phi``
drops out in the vorticity equation (``\nabla \times \nabla \Phi = 0``),
but becomes a ``-\nabla^2 \Phi`` in the divergence equation,
which is therefore combined with the kinetic energy term
``-\nabla^2(\tfrac{1}{2}(u^2 + v^2))`` similar as it is done in the [Shallow water equations](@ref).
You can think of ``\tfrac{1}{2}(u^2 + v^2) + \Phi`` as the Bernoulli potential in
the primitive equations. However, due to the change into sigma coordinates the surface pressure
gradient also has to be accounted for. Now highlighting only the pressure gradient force, we
have in total
```math
\begin{aligned}
\frac{\partial \zeta}{\partial t} &= \nabla \times (... - R_dT_v\nabla \ln p_s) + ... \\
\frac{\partial \mathcal{D}}{\partial t} &= \nabla \cdot (... - R_dT_v\nabla \ln p_s) - \nabla^2\Phi + ...
\end{aligned}
```
In our vorticity-divergence formulation and with sigma coordinates.

### Semi-implicit pressure gradient

With the [semi-implicit time integration](@ref implicit_primitive) in SpeedyWeather.jl
the pressure gradient terms are further modified as follows. See that section for details
why, but here is just to mention that we need to split the terms into linear and non-linear
terms. The linear terms are then evaluated at the previous time step for the implicit
scheme such that we can avoid instabilities from gravity waves.

We split the (virtual) temperature into a reference vertical profile ``T_k`` and its anomaly,
``T_v = T_k + T_v'``. The reference profile ``T_k`` has to be a global
constant for the spectral transform but can depend on the vertical.
With this, the previous equation becomes
```math
\begin{aligned}
\frac{\partial \zeta}{\partial t} &= \nabla \times (... - R_dT_v'\nabla \ln p_s) + ... \\
\frac{\partial \mathcal{D}}{\partial t} &= \nabla \cdot (... - R_dT_v'\nabla \ln p_s) - \nabla^2(\Phi + R_d T_k \ln p_s) + ...
\end{aligned}
```
In the vorticity equation the term with the reference profile drops out as ``\nabla \times \nabla = 0``,
and in the divergence equation we move it into the Laplace operator. Now the linear terms
are gathered with the Laplace operator and for the semi-implicit scheme we calculate both the
[Geopotential](@ref) ``\Phi`` and the contribution to the "linear pressure gradient"
``R_dT_k \ln p_s`` at the previous time step for the
[semi-implicit time integration](@ref implicit_primitive) for details see therein.

## Vorticity advection

Vorticity advection in the primitive equation takes the form
```math
\begin{aligned}
\frac{\partial u}{\partial t} &= (f+\zeta)v \\
\frac{\partial v}{\partial t} &= -(f+\zeta)u \\
\end{aligned}
```
Meaning that we add the Coriolis parameter ``f`` and the relative vorticity ``\zeta``
and multiply by the respective velocity component. While the primitive equations here
are written with vorticity and divergence, we use ``u,v`` here as other tendencies
will be added and the curl and divergence are only taken once after transform into
spectral space. To obtain a tendency for vorticity and divergence, we rewrite this as
```math
\begin{aligned}
\frac{\partial \zeta}{\partial t} &= \nabla \times (f+\zeta)\mathbf{u}_\perp \\
\frac{\partial \mathcal{D}}{\partial t} &= \nabla \cdot (f+\zeta)\mathbf{u}_\perp \\
\end{aligned}
```
with ``\mathbf{u}_\perp = (v,-u)`` the rotated velocity vector, see [Barotropic vorticity equation](@ref).

## Humidity equation

The dynamical core treats humidity as an (active) tracer, meaning that after the physical
parameterizations for humidity ``\mathcal{P}`` are calculated in grid-point space,
humidity is only advected with the flow. The only exception is the [Virtual temperature](@ref)
as high levels of humidity will lower the effective density, which is why we use the
virtual instead of the absolute temperature. The equation to be solved for humidity
is therefore,
```math
\left( \frac{\partial q}{\partial t} \right) = \left(\left[\mathcal{P}_q - W_q +
q\mathcal{D} \right]\right) -\nabla\cdot([\mathbf{u}q])
```
With ``()`` denoting spectral space and ``[]`` grid-point space, so that
``([])`` and ``[()]`` are the transforms in the respective directions.
To avoid confusion with that notation, we write the tendency of humidity due
to [Vertical advection](@ref) as ``W_q``. This equation is identical to a tracer equation,
with ``\mathcal{P}_q`` denoting sources and sinks. Note that [Horizontal diffusion](@ref)
should be applied to every advected variable.

A very similar equation is solved for (absolute) temperature
as described in the following.

## Temperature equation

The first law of thermodynamic states that the internal energy ``I`` is increased by
the heat ``Q`` applied minus the work ``W`` done by the system. We neglect changes
in chemical composition ([^Vallis], chapter 1.5). For an ideal gas, the internal 
energy is ``c_vT`` with ``c_v`` the heat capacity at constant volume and temperature
``T``. The work done is ``pV``, with pressure ``p`` and the specific volume ``V``
```math
dI = Q - p dV.
```
For fluids we replace the differential ``d`` here with the material derivative ``\tfrac{D}{Dt}``.
With ``V = \tfrac{1}{\rho}`` and density ``\rho`` we then have
```math
c_v \frac{DT}{Dt} = -p \frac{D (1/\rho)}{Dt} + Q
```
Using the ideal gas law to replace ``\tfrac{1}{\rho}`` with ``\tfrac{RT_v}{p}`` (we are
using the [Virtual temperature](@ref) again), and using
```math
p\frac{D (1/p)}{Dt} = -\frac{1}{p} \frac{Dp}{Dt}
```
we have
```math
(c_v + R)\frac{DT}{Dt} = \frac{RT_v}{p}\frac{Dp}{Dt} + Q
```
And further, with ``c_p = c_v + R`` the heat capacity at constant pressure,
``\kappa = \tfrac{R}{c_p}``, and using the logarithm of pressure
```math
\frac{DT}{Dt} = \kappa T_v\frac{D \ln p}{Dt} + \frac{Q}{c_p}
```
This is the form of the temperature equation that SpeedyWeather.jl uses. Temperature
is advected through the material derivative and first term on the right-hand side
represents an adiabatic conversion term describing how the temperature changes with
changes in pressure. Recall that this term originated from the work term in
the first law of thermodynamics. The forcing term ``\tfrac{Q}{c_p}`` is here
identified as the physical parameterizations changing the temperature, for example
radiation, and hence we will call it ``P_T``.

Similar to the [Humidity equation](@ref) we write the equation for (absolute) temperature ``T`` as

```math
\left( \frac{\partial T}{\partial t} \right) = \left(\left[\mathcal{P}_T - W_T +
T\mathcal{D} + \kappa T_v \frac{D \ln p}{Dt} \right]\right) -\nabla\cdot([\mathbf{u}T])
```
``W_T`` is the [Vertical advection](@ref) of temperature. We evaluate the adiabatic conversion
term completely in grid-point space following Simmons and Burridge, 1981[^SB81] Equation 3.12 and 3.13.
Leaving out the ``\kappa T_v`` for clarity, the term at level ``k`` is
```math
\left(\frac{D \ln p}{D t}\right)_k = \mathbf{u}_k \cdot \nabla \ln p_k
- \frac{1}{\Delta p_k} \left[\left( \ln \frac{p_{k+\tfrac{1}{2}}}{p_{k-\tfrac{1}{2}}}\right)
\sum_{r=1}^{k-1}\nabla \cdot (\mathbf{u}_k \Delta p_k) + \alpha_k \nabla \cdot (\mathbf{u}_k \Delta p_k) \right]
```
with
```math
\alpha_k = 1 - \frac{p_{k-\tfrac{1}{2}}}{\Delta p_k} \ln \frac{p_{k+\tfrac{1}{2}}}{p_{k-\tfrac{1}{2}}}
```
In sigma coordinates this simplifies to, following similar steps as in [Surface pressure tendency](@ref)
```math
\begin{aligned}
\left(\frac{D \ln p}{D t}\right)_k &= \mathbf{u}_k \cdot \nabla \ln p_s \\
&- \frac{1}{\Delta \sigma_k} \left( \ln \frac{\sigma_{k+\tfrac{1}{2}}}{\sigma_{k-\tfrac{1}{2}}}\right)
\sum_{r=1}^{k-1}\Delta \sigma_r (\mathcal{D}_r + \mathbf{u}_r \cdot \nabla \ln p_s) -
\alpha_k (\mathcal{D}_k + \mathbf{u}_k \cdot \nabla \ln p_s)
\end{aligned}
```
Let ``A_k = \mathcal{D}_k + \mathbf{u}_k \cdot \nabla \ln p_s`` and
``\beta_k = \tfrac{1}{\Delta \sigma_k} \left( \ln \tfrac{\sigma_{k+\tfrac{1}{2}}}{\sigma_{k-\tfrac{1}{2}}}\right)``,
then this can also be summarised as
```math
\left(\frac{D \ln p}{D t}\right)_k = \mathbf{u}_k \cdot \nabla \ln p_s
- \beta_k \sum_{r=1}^{k-1}\Delta \sigma_r A_r - \alpha_k A_k
```
The ``\alpha_k, \beta_k`` are constants and can be precomputed. The surface pressure flux
``\mathbf{u}_k \cdot \nabla \ln p_s`` has to be computed, so does the vertical sigma-weighted
average from top to ``k-1``, which is done when computing other vertical averages for the
[Surface pressure tendency](@ref).

### Semi-implicit temperature equation

For the [semi-implicit scheme](@ref implicit_primitive) we need to split the temperature
equation into linear and non-linear terms, as the linear terms need to be evaluated
at the previous time step. Decomposing temperature ``T`` into ``T = T_k + T'``
with the reference profile ``T_k`` and its anomaly ``T'``, the temperature equation
becomes
```math
\left( \frac{\partial T}{\partial t} \right) = \mathcal{P}_T - W_T +
T'\mathcal{D} + \kappa T_v \frac{D \ln p}{Dt} -\nabla\cdot(\mathbf{u}T')
```
Note that we do not change the adiabatic conversion term. While its linear
component ``\kappa T_k^v \tfrac{D \ln p_s}{D t}`` (the subscript ``v`` for
[Virtual temperature](@ref) as been raised)
would need to be evaluated at the previous time step, we still evaluate this
term at the current time step and move it within the semi-implicit corrections
to the previous time step afterwards.

## [Semi-implicit time stepping](@id implicit_primitive)

Conceptually, the semi-implicit time stepping in the [Primitive equation model](@ref) is
the same as in the [Shallow water model](@ref implicit_swm),
but

- tendencies for divergence ``\mathcal{D}``, logarithm of surface pressure ``\ln p_s`` but also temperature ``T`` are computed semi-implicitly,
- the vertical layers are coupled, creating a linear equation system that is solved via matrix inversion.

The linear terms of the primitive equations follow a linearization around a state of rest without
orography and a reference vertical temperature profile. The scheme described here largely follows
Hoskins and Simmons [^HS75], which has also been used in Simmons and Burridge [^SB81].

As before, let ``\delta V = \tfrac{V_{i+1} - V_{i-1}}{2\Delta t}`` be the tendency we need for the Leapfrog
time stepping. With the implicit time step ``\xi = 2\alpha\Delta t``, ``\alpha \in [\tfrac{1}{2},1]`` we have

```math
\delta V = N_E(V_i) + N_I(V_{i-1}) + \xi N_I(\delta V)
```
with ``N_E`` being the explicitly-treated non-linear terms and ``N_I`` the implicitly-treated linear terms, such that
``N_I`` is a linear operator. We can therefore solve for ``\delta V``
by inverting ``N_I``, 
```math
\delta V = (1-\xi N_I)^{-1}G
```
where we gathered the uncorrected right-hand side as ``G``
```math
G = N_E(V_i) + N_I(V_{i-1}) = N(V_i) + N_I(V_{i-1} - V_i).
```
So for every linear term in ``N_I`` we have two options corresponding to two sides of this equation

1. Evaluate it at the previous time step ``i-1``
2. Or, evaluate it at the current time step ``i`` as ``N(V_i)``, but then move it back to the previous time step ``i-1`` by adding (in spectral space) the linear operator ``N_I`` evaluated with the difference between the two time steps.

If there is a tendency that is easily evaluated in spectral space it is easier to follow 1. However,
a term that is costly to evaluate in grid-point space should usually follow the latter. The reason is that
the previous time step is generally not available in grid-point space (unless recalculated through a costly transform
or stored with additional memory requirements) so it is easier to follow 2 where the ``N_I`` is
available in spectral space. For the adiabatic conversion term in the [Temperature equation](@ref)
we follow 2 as one would otherwise need to split this term into a non-linear and linear term,
evaluating it essentially twice in grid-point space. 

So what is ``G`` in the [Primitive equation model](@ref)? 

```math
\begin{aligned}
G_\mathcal{D} &= N^E_\mathcal{D} - \nabla^2(\Phi^{i-1} + R_dT_k^v (\ln p_s)^{i-1})
= N^E_\mathcal{D} - \nabla^2( \mathbf{R}T^{i-1} + \mathbf{U}\ln p_s^{i-1}) \\
G_{\ln p_s} &= N_{\ln p_s}^E - \overline{\mathcal{D}^{i-1}}
= N_{\ln p_s}^E + \mathbf{W}\mathcal{D}^{i-1} \\
G_T &= N_T + \mathbf{L}(\mathcal{D}^{i-1} - \mathcal{D}^i) \\
\end{aligned}
```
``G`` is for the divergence, pressure and temperature equation the "uncorrected" tendency.
Moving time step ``i - 1 \to i`` we would be back with a fully explicit scheme. In the divergence equation
the [Geopotential](@ref) ``\Phi`` is calculated from temperature ``T`` at the previous time step
``i-1`` (denoted as superscript) and the "linear" [Pressure gradient](@ref) from the logarithm of
surface pressure at the previous time step.
One can think of these two calculations as linear operators, ``\mathbf{R}`` and ``\mathbf{U}``.
We will shortly discuss their properties. While we could combine them with the Laplace operator
``\nabla^2`` (which is also linear) we do not do this as ``\mathbf{R, U}`` do not depend on
the degree and order of the spherical harmonics (their *wavenumber*) but on the vertical,
but ``\nabla^2`` does not depend on the vertical, only on the wavenumber.
All other terms are gathered in ``N_\mathcal{D}^E`` (subscript ``E`` has been raised)
and calculated as described in the respective section at the current time step ``i``.

For the pressure tendency, the subtraction with the thickness-weighted vertical average ``\bar{\mathcal{D}}``
is the linear term that is treated implicitly. We call this operator ``\mathbf{W}``.
For the temperature tendency, we evaluate all terms explicitly at the current time step in ``N_T``
but then move the linear term in the adiabatic conversion term with the operator ``\mathbf{L}``
back to the previous time step. For details see [Semi-implicit temperature equation](@ref).

The operators ``\mathbf{R, U, L, W}`` are all linear, meaning that we can apply them 
in spectral space to each spherical harmonic independently -- the vertical is coupled however.
With ``N`` being the number of vertical levels and the prognostic variables like
temperature for a given degree ``l`` and order ``m`` being a column vector in the vertical,
``T_{l,m} \in \mathbb{R}^N``, these operators have the following shapes

```math
\begin{aligned}
\mathbf{R} &\in \mathbb{R}^{N\times N} \\
\mathbf{U} &\in \mathbb{R}^{N\times 1} \\
\mathbf{L} &\in \mathbb{R}^{N\times N} \\
\mathbf{W} &\in \mathbb{R}^{1\times N} \\
\end{aligned}
```

``\mathbf{R}`` is an integration in the vertical hence it is an upper triangular matrix such that
the first (an top-most) ``k=1`` element of the resulting vector depends on all vertical levels
of the temperature mode ``T_{l,m}``, but the surface ``k=N`` only on the temperature mode at the surface.
``\mathbf{U}`` takes the surface value of the ``l,m`` mode of the logarithm of surface pressure
``(\ln p_s)_{l,m}`` and multiplies it element-wise with the reference temperature profile and
the dry gas constant. So the result is a column vector.
``\mathbf{L}`` is an ``N \times N`` matrix as the adiabatic conversion term couples all layers.
``\mathbf{W}`` is a row vector as it represents the vertical averaging of the spherical harmonics
of a divergence profile. So, ``\mathbf{W}\mathcal{D}`` is a scalar product for every ``l,m``
giving a contribution of all vertical layers in divergence to the (single-layer!) logarithm of
surface pressure tendency.

With the ``G``s defined we can now write the semi-implicit tendencies ``\delta \mathcal{D}``,
``\delta T``, ``\delta \ln p_s`` as (first equation in this section)
```math
\begin{aligned}
\delta \mathcal{D} &= G_D - \xi \nabla^2(\mathbf{R}\delta T + \mathbf{U} \delta \ln p_s)\\
\delta T &= G_T + \xi \mathbf{L}\delta \mathcal{D} \\
\delta \ln p_s &= G_{\ln p_s} + \xi \mathbf{W}\delta \mathcal{D}
\end{aligned}
```
Solving for ``\delta \mathcal{D}`` with the "combined" tendency
```math
G = G_D - \xi \nabla^2(\mathbf{R}G_T + \mathbf{U}G_{\ln p_s})
```
via
```math
\delta \mathcal{D} = G - \xi^2\nabla^2(\mathbf{RL + UW})\delta \mathcal{D}
```
(``\mathbf{UW}`` is a matrix of size ``N \times N``) yields
```math
\delta D = \left( 1 + \xi^2\nabla^2(\mathbf{RL + UW})  \right)^{-1}G = \mathbf{S}^{-1}G
```
The other tendencies ``\delta T`` and ``\delta \ln p_s`` are then obtained
through insertion above. We may call the operator to be inverted ``\mathbf{S}``
which is of size ``l_{max} \times N \times N``, hence for every degree ``l`` of
the spherical harmonics (which the Laplace operator depends on) a
``N \times N`` matrix coupling the ``N`` vertical levels. Furthermore, ``S`` depends
- through ``\xi`` on the time step ``\Delta t``,
- through ``\mathbf{R,W,L}`` on the vertical level spacing ``\Delta \sigma_k``
- through ``\mathbf{U}`` on the reference temperature profile ``T_k``

so for any changes of these the matrix inversion of ``\mathbf{S}`` has to
be recomputed. Otherwise the algorithm for the semi-implicit scheme is as follows

0\. Precompute the linear operators ``\mathbf{R,U,L,W}`` and with them the matrix inversion ``\mathbf{S}^{-1}``.

Then for every time step

1. Compute the uncorrected tendencies evaluated at the current time step for the explicit terms and the previous time step for the implicit terms.
2. Exception in SpeedyWeather.jl is the adiabatic conversion term, which is, using ``\mathbf{L}`` moved afterwards from the current ``i`` to the previous time step ``i-1``.
3. Compute the combined tendency ``G`` from the uncorrected tendencies ``G_\mathcal{D}``, ``G_T``, ``G_{\ln p_s}``.
4. With the inverted operator get the corrected tendency for divergence, ``\delta \mathcal{D} = \mathbf{S}^{-1}G``.
5. Obtain the corrected tendencies for temperature ``\delta T`` and surface pressure ``\delta \ln p_s`` from ``\delta \mathcal{D}``.
6. Apply [Horizontal diffusion](@ref) (which is only mentioned here as it further updates the tendencies).
7. Use ``\delta \mathcal{D}``, ``\delta T`` and ``\delta \ln p_s`` in the [Leapfrog time integration](@ref leapfrog).

## Horizontal diffusion

Horizontal diffusion in the primitive equations is applied to vorticity ``\zeta``, divergence ``\mathcal{D}``,
temperature ``T`` and humidity ``q``. In short, all variables that are advected.
For the dry equations, ``q=0`` and no diffusion has to be applied.

The horizontal diffusion is applied implicitly in spectral space, as already described in
[Horizontal diffusion](@ref) for the barotropic vorticity equation.

## Algorithm

The following algorithm describes a time step of the `PrimitiveWetModel`, for
the `PrimitiveDryModel` humidity can be set to zero and respective steps skipped.

0\. Start with initial conditions of relative vorticity ``\zeta_{lm}``, divergence ``D_{lm}``,
temperature ``T_{lm}``, humidity ``q_{lm}`` and the logarithm of surface pressure ``(\ln p_s)_{lm}``
in spectral space. Variables ``\zeta, D, T, q`` are defined on all vertical levels, the logarithm
of surface pressure only at the surface. Transform this model state to grid-point space,
obtaining velocities is done as in the shallow water model

- Invert the [Laplacian](@ref) of ``\zeta_{lm}`` to obtain the stream function ``\Psi_{lm}`` in spectral space
- Invert the [Laplacian](@ref) of ``D_{lm}`` to obtain the velocity potential ``\Phi_{lm}`` in spectral space
- obtain velocities ``U_{lm} = (\cos(\theta)u)_{lm}, V_{lm} = (\cos(\theta)v)_{lm}`` from ``\nabla^\perp\Psi_{lm} + \nabla\Phi_{lm}``
- Transform velocities ``U_{lm}``, ``V_{lm}`` to grid-point space ``U,V``
- Unscale the ``\cos(\theta)`` factor to obtain ``u,v``

Additionally we

- Transform ``\zeta_{lm}``, ``D_{lm}``, ``T_{lm}, (\ln p_s)_{lm}`` to ``\zeta, D, \eta, T, \ln p_s`` in grid-point space
- Compute the (non-linearized) [Virtual temperature](@ref) in grid-point space.

Now loop over

1. Compute all tendencies of ``u, v, T, q`` due to physical parameterizations in grid-point space.
2. Compute the gradient of the logarithm of surface pressure ``\nabla (\ln p_s)_{lm}`` in spectral space and convert the two fields to grid-point space. Unscale the ``\cos(\theta)`` on the fly.
3. For every layer ``k`` compute the pressure flux ``\mathbf{u}_k \cdot \nabla \ln p_s`` in grid-point space. 
4. For every layer ``k`` compute a linearized [Virtual temperature](@ref) in spectral space.
5. For every layer ``k`` compute a temperature anomaly (virtual and absolute) relative to a vertical reference profile ``T_k`` in grid-point space.
6. Compute the [Geopotential](@ref) ``\Phi`` by integrating the virtual temperature vertically in spectral space from surface to top.
7. Integrate ``u,v,D`` vertically to obtain ``\bar{u},\bar{v},\bar{D}`` in grid-point space and also ``\bar{D}_{lm}`` in spectral space. Store on the fly also for every layer ``k`` the partial integration from 1 to ``k-1`` (top to layer above). These will be used in the adiabatic term of the [Temperature equation](@ref).
8. Compute the [Surface pressure tendency](@ref) with the vertical averages from the previous step. For the [semi-implicit time stepping](@ref implicit_primitive)
9. For every layer ``k`` compute the [Vertical velocity](@ref).
10. For every layer ``k`` add the linear contribution of the [Pressure gradient](@ref) ``RT_k (\ln p_s)_{lm}`` to the geopotential ``\Phi`` in spectral space.
11. For every layer ``k`` compute the [Vertical advection](@ref) for ``u,v,T,q`` and add it to the respective tendency.
12. For every layer ``k`` compute the tendency of ``u,v`` due to [Vorticity advection](@ref) and the [Pressure gradient](@ref) ``RT_v \nabla \ln p_s`` and add to the respective existing tendency. Unscale ``\cos(\theta)``, transform to spectral space, take curl and divergence to obtain tendencies for ``\zeta_{lm},\mathcal{D}_{lm}``.
13. For every layer ``k`` compute the adiabatic term and the horizontal advection in the [Temperature equation](@ref) in grid-point space, add to existing tendency and transform to spectral.
14. For every layer ``k`` compute the horizontal advection of humidity ``q`` in the [Humidity equation](@ref) in grid-point space, add to existing tendency and transform to spectral.
15. For every layer ``k`` compute the kinetic energy ``\tfrac{1}{2}(u^2 + v^2)``, transform to spectral and add to the [Geopotential](@ref). For the [semi-implicit time stepping](@ref implicit_primitive) also add the linear pressure gradient calculated from the previous time step. Now apply the Laplace operator and subtract from the divergence tendency.
16. Correct the tendencies following the [semi-implicit time integration](@ref implicit_swm) to prevent fast gravity waves from causing numerical instabilities.
17. Compute the [horizontal diffusion](@ref diffusion) for the advected variables ``\zeta,\mathcal{D},T,q``
18. Compute a leapfrog time step as described in [Time integration](@ref leapfrog) with a [Robert-Asselin and Williams filter](@ref)
19. Transform the new spectral state of ``\zeta_{lm}``, ``\mathcal{D}_{lm}``, ``T_{lm}``, ``q_{lm}`` and ``(\ln p_s)_{lm}`` to grid-point ``u,v,\zeta,\mathcal{D},T,q,\ln p_s`` as described in 0.
20. Possibly do some output
21. Repeat from 1.

## Scaled primitive equations

## References

[^GFDL1]: Geophysical Fluid Dynamics Laboratory, [Idealized models with spectral dynamics](https://www.gfdl.noaa.gov/idealized-models-with-spectral-dynamics/)
[^GFDL2]: Geophysical Fluid Dynamics Laboratory, [The Spectral Dynamical Core](https://www.gfdl.noaa.gov/wp-content/uploads/files/user_files/pjp/spectral_core.pdf)
[^Vallis]: GK Vallis, 2006. [Atmopsheric and Ocean Fluid Dynamics](http://vallisbook.org/), Cambridge University Press.
[^SB81]: Simmons and Burridge, 1981. *An Energy and Angular-Momentum Conserving Vertical Finite-Difference Scheme and Hybrid Vertical Coordinates*, Monthly Weather Review. DOI: [10.1175/1520-0493(1981)109<0758:AEAAMC>2.0.CO;2](https://doi.org/10.1175/1520-0493(1981)109<0758:AEAAMC>2.0.CO;2).
[^HS75]: Hoskins and Simmons, 1975. *A multi-layer spectral model and the semi-implicit method*, Quart. J. R. Met. Soc. DOI: [10.1002/qj.49710142918](https://doi.org/10.1002/qj.49710142918)