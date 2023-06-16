# Primitive equation model

The [primitive equations](https://en.wikipedia.org/wiki/Primitive_equations) are a hydrostatic approximation
of the compressible Navier-Stokes equations for an ideal gas on a rotating sphere. We largely follow
the idealised spectral dynamical core developed by GFDL[^1] and documented therein[^2].

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
Coriolis parameter ``f``, ``W`` the vertical advection operator, dry air gas constant ``R_d``,
virtual temperature ``T_v``, geopotential ``\Phi``, pressure ``p``, thermodynamic ``\kappa = R\_d/c_p``
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
obtained from the univeral gas constant ``R`` divided by the molecular masses
of the gas. The total pressure ``p`` in the air parcel is
```math
p = p_d + p_w = (\rho_d R_d + \rho_w R_w)T
```
We ultimately want to replace the density ``\rho = \rho_w + \rho_d`` in the dynamical core,
using the ideal gas law, with the temperature ``T``, so that we never have
to calculate the density explicitly. However, in order to not deal with
two densities (dry air and water vapour) we would like to replace
temperature with a virtual temperature that includes the effect of
humidity on the density. So, whereever we use the ideal gas law
to replace density with temperature, we would use the virtual temperature,
which is a function of the absolute temperature and specific humidity,
instead. A higher specific humidity in an air parcel lowers 
the density as water vapour is lighter than dry air. Consequently,
the virtual temperature of moist air is higher than its absolute temperature
because warmer air is lighter too at constant pressure. We therefore
think of the virtual temperature as the temperature dry air would need to have
to be as light as moist air.

Starting with the last equation, with some manipulation we can write
the ideal gas law as total density ``rho`` times a gas constant
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
product of two variables ``q,T`` has to be computed in grid-point space, see [Spectral Transform].
To obtain an approximation to the virtual temperature in spectral space without
expensive transforms one can linearize
```math
T_v = T + \mu q\bar{T}
```
With a global constant temperature ``\bar{T}``, for example obtained from the
``l=m=0`` mode, ``\bar{T} = T_{0,0}\frac{1}{\sqrt{4\pi}}`` but depending on the
normalization of the spherical harmonics that factor needs adjustment.

## Vertical coordinates

### General

Let ``\Psi(x,y,z,t)`` 

SpeedyWeather.jl currently uses sigma coordinates for the vertical. 


```math
\sigma = \frac{p}{p_s}
```

```math
p_k = \sigma_kp_s
```
```math
\Delta p_k = p_{k+1} - p_k = \Delta \sigma_k p_s
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
```math
\Phi_{k-\tfrac{1}{2}} = \Phi_{k+\tfrac{1}{2}} + R_dT^v_k(\ln p_{k+1/2} - \ln p_{k-1/2})
```

## Surface pressure tendency

## Vertical advection

## Pressure gradient force

## Temperature equation

## [Semi-implicit time stepping](@id implicit_primitive)

## Horizontal diffusion

## Algorithm

## Scaled primitive equations

## References

[^1]: Geophysical Fluid Dynamics Laboratory, [Idealized models with spectral dynamics](https://www.gfdl.noaa.gov/idealized-models-with-spectral-dynamics/)
[^2]: Geophysical Fluid Dynamics Laboratory, [The Spectral Dynamical Core](https://www.gfdl.noaa.gov/wp-content/uploads/files/user_files/pjp/spectral_core.pdf)