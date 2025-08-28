# Large-scale condensation

Large-scale condensation in an atmospheric general circulation represents the
micro-physics that kick in when an air parcel reaches saturation.
Subsequently, the water vapour inside it condenses, forms droplets around
condensation nuclei, which grow, become heavy and eventually fall out
as precipitation. This process is never actually representable at the resolution
of global (or even regional) atmospheric models as typical cloud droplets
have a size of micrometers. Atmospheric models therefore rely on large-scale
quantities such as specific humidity, pressure and temperature within a
given grid cell, even though there might be considerable variability of
these quantities within the area covered by that grid cell.
Higher resolution would use more grid cells within the same area,
resolving more variability that averaged out at lower resolution.

## Condensation implementations

Currently implemented are

```@example condensation
using InteractiveUtils # hide
using SpeedyWeather
subtypes(SpeedyWeather.AbstractCondensation)
```

which are described in the following.

## Explicit large-scale condensation

We parameterize this process of _large-scale condensation_ when relative humidity
in a grid cell reaches saturation and remove the excess humidity quickly
(given time integration constraints, see below) and with an implicit (in the time
integration sense) latent heat release. Vertically integrating the
tendency of specific humidity due to this process is then the _large-scale precipitation_.

Immediate condensation of humidity ``q_i > q^\star`` at time step ``i`` given its saturation ``q^\star``
humidity calculated from temperature ``T_i`` is

```math
\begin{aligned}
q_{i+1} - q_i &= q^\star(T_i) - q_i \\
T_{i+1} - T_i &= -\frac{L_v}{c_p}( q^\star(T_i) - q_i  )
\end{aligned}
```

This condensation is _explicit_ in the time integration sense, meaning that we only use
quantities at time step ``i`` to calculate the tendency. The latent heat release of that
condensation is in the second equation. However, treating this explicitly poses the problem
that because the saturation humidity is calculated from the current temperature ``T_i``,
which is increased due to the latent heat release, the humidity after this time step will be
undersaturated.

## Implicit large-scale condensation

Ideally, one would want to condense towards the _new_ saturation humidity
``q^\star(T_{i+1})`` at ``i+1`` so that condensation draws the relative humidity back down
to 100% not below it.  Taylor expansion at ``i`` of the equation above with ``q^\star(T_{i+1})``
and ``\Delta T = T_{i+1} - T_i`` (and ``\Delta q`` similarly) to first order yields

```math
q_{i+1} - q_i = q^\star(T_{i+1}) - q_i = q^\star(T_i) + (T_{i+1} - T_i)
\frac{\partial q^\star}{\partial T} (T_i) + O(\Delta T^2) - q_i
```

Now we make a linear approximation to the derivative and drop the ``O(\Delta T^2)`` term.
Inserting the (explicit) latent heat release yields

```math
\Delta q = q^\star(T_i) + -\frac{L_v}{c_p} \Delta q \frac{\partial q^\star}{\partial T} (T_i) - q_i
```

And solving for ``\Delta q`` yields

```math
\left[ 1 + \frac{L_v}{c_p} \frac{\partial q^\star}{\partial T} (T^i) \right] \Delta q = q^\star(T_i) - q_i
```

meaning that the implicit immediate condensation can be formulated as (see also [^Frierson2006])

```math
\begin{aligned}
q_{i+1} - q_i &= \frac{q^\star(T_i) - q_i}{1 + \frac{L_v}{c_p} \frac{\partial q^\star}{\partial T}(T_i)} \\
T_{i+1} - T_i &= -\frac{L_v}{c_p}( q_{i+1} - q_i )
\end{aligned}
```

With Euler forward time stepping this is great, but with our [leapfrog timestepping + RAW filter](@ref leapfrog)
this is very dispersive (see [#445](https://github.com/SpeedyWeather/SpeedyWeather.jl/issues/445))
although the implicit formulation is already much better.
We therefore introduce a time step ``\Delta t_c`` which makes the implicit condensation not immediate
anymore but over several time steps ``\Delta t`` of the leapfrogging.

```math
\begin{aligned}
\delta q = \frac{q_{i+1} - q_i}{\Delta t} &= \frac{q^\star(T_i) - q_i}{ \Delta t_c
\left( 1 + \frac{L_v}{c_p} \frac{\partial q^\star}{\partial T}(T_i) \right)} \\
\delta T = \frac{T_{i+1} - T_i}{\Delta t} &= -\frac{L_v}{c_p}( \frac{q_{i+1} - q_i}{\Delta t} )
\end{aligned}
```

For ``\Delta t = \Delta t_c`` we have an immediate condensation,
for ``n = \frac{\Delta t_c}{\Delta t}`` condensation takes place over ``n`` time steps.
One could tie this time scale for condensation to a physical unit, like 6 hours,
but because the time step here is ideally short, but cannot be too short for
numerical stability, we tie it here to the time step of the numerical integration.
This also means that at higher resolution condensation is more immediate than
at low resolution, but the dispersive time integration of this term is in all
cases similar (and not much higher at lower resolution).

The above implied that condensation takes place at 100% relative humidity as
we are directly comparing the specific humidity to the saturation specific humidity
in ``q^\star - q_i``. However, in coarse-resolution models there is a good argument
that condensation may already take place below a grid-cell average of 100% relative humidity.
Given some subgrid-scale variability parts of the grid cells may condense even though
the average humidity is below 100%. To change this threshold one can introduce
a parameter ``r``, e.g. ``r=0.95`` for condensation to take place at 95% relative humidity,
as follows

```math
\delta q = \frac{q_{i+1} - q_i}{\Delta t} = \frac{rq^\star(T_i) - q_i}{ \Delta t_c
\left( 1 + \frac{L_vr}{c_p} \frac{\partial q^\star}{\partial T}(T_i) \right)}
```
``r`` is a linear scale and therefore can be taken out of the gradient
``\frac{\partial q^\star}{\partial T}`` in the denominator.

## Re-evaporation

Reevaporation is a process that requires to compute the downward rain water flux ``F_r``
iteratively from the top layer to the bottom layer, we therefore start with ``F_r = 0``
at the top of the atmosphere. We will use it in units of ``m/s``, that is the rainfall rate.
But you also may use ``kg/m^2/s`` the mass flux of rain water per second.
Which would need to be divided by the water density ``\rho`` to convert into a rain
rate in ``m/s``, or then multiply with the time step ``\Delta t`` and to get
a rainfall amount in meters.

When rainfall is created in one layer but falls through a drier (and often warmer) layer below
then the rain water can re-evaporate, effectively causing a humidity flux into lower layers, and reducing the amount of rain that reaches the ground. We parameterize this effect
proportional to the difference of humidity ``q`` to saturation ``q^*``.

```math
\begin{aligned}
F_{e, k}   &= \min(c\max(q^* - q, 0), 1) F_{r, k} \\
F_{r, k+1} &= F_{r, k} - F_{e, k}
\end{aligned}
```

So ``F_e`` is the evaporated rain in layer ``k`` as a fraction of ``F_r`` the rain water
flux into layer ``k`` from above (``k`` increases top to bottom). ``c`` is a proportionality
constant that effectively parameterizes the effectiveness of reevaporation.
For ``c=0`` reevaporation is disabled, for ``c = 1/q^*`` (though ``c`` is a global constant)
and ``q=0`` evaporation would be immediate, i.e. when the humidity in the layer is zero
but ``c`` is chosen to be on the order of  one over the saturation humidity then all
rain water would evaporate. We add a ``\min(..., 1)`` to avoid evaporating more rain water
than is available. We then convert ``F_e`` into a humidity tendency by division with
``\Delta p/(\Delta t g \rho)``, layer pressure thickness ``\Delta p``, gravity ``g``,
but this depends on the units you use for the rain water flux.
Since reevaporation takes the same latent heat (though opposite sign) as condensation,
this tendency has to be substracted from the large-scale condensation tendency but the
implicit time stepping can be used as before if reevaporation is calculated before the
latent heat release.

The reevaportation in `ImplicitCondensation` is controlled by `reevaporation` (dimensionless),
the proportionality constant ``c \leq 0`` here. 

```@example condensation
spectral_grid = SpectralGrid()
large_scale_condensation = ImplicitCondensation(spectral_grid, reevaporation=0)
```

would disable reevaporation `reevaporation = 30` instead would be roughly equivalent
to immediate evaporation at 30˚C and zero ambient humidity near surface.

## Snow fall

We parameterize snow fall if large-scale condensation occurs below a freezing temperature ``T_f``,
we currently use ``T_f = -10˚C`` but this may change in future versions, check `model.large_scale_condensation`
for that. Freezing of rain water to snow occurs immediately, so we move the rain water flux
to the snow flux ``F_s = F_r`` and set the ``F_r = 0`` afterwards, if temperature in that layer
is below freezing ``T < T_f``. The latent heat release from freezing is then

```math
\delta T_f = -\frac{L_i}{c_p} \delta q_f
```

with latent heat of fusion ``L_i``. As ``L_i`` is an order of magnitude smaller than the latent heat
of vaporization we add this temperature tendency explicitly in time. The minus sign is to denote
that a negative humidity tendency due to freezing condensation (=snow) should release latent heat
and therefore warm the air. ``\delta q_f`` is ``\min(0, \delta q)`` if the freezing condition is
met (negative because condensation is a negative humidity tendency and to exclude a tendency
that is dominated by reevaporation), ``\delta q_f = 0`` otherwise.

Like rain water can reevaporate in the layer so can snow melt in the layer below and turn into rain again.
Precipitation often occurs in layers high up in the atmosphere where water may be subject to freezing
but if the atmopshere below is warm then there is enough energy available to melt the snow before it
reaches the ground. We want to parameterize this effect (otherwise it can easily snow in the tropics).

The available energy ``E_m`` in ``J/kg`` (of air) for melting snow is proportional to the temperature
of the air above a melt threshold ``T_m``

```math
E_m = c_p \max(T - T_m, 0)
```

We can convert this to a maximum melt rate, i.e. the amount of snow that this energy is able to melt
in one time step. Note that for rate we use units of rain water height (or depth)
in meters per second what we also use for the snow and rain water fluxes ``F_r, F_s``. So this is not
the actual snow height as it would have on the ground but if melted to water.

```math
F_m = \min(F_s, \frac{E_m}{L_i}\frac{\Delta p}{\Delta t g\rho})
```

We cap this to ``F_s`` so that one cannot melt more snow than there is. This snow melt flux
is then subtracted from ``F_s`` but added to ``F_r`` to move snow water to rain water.
The according latent heat required for melting is 

```math
\begin{aligned}
\delta q_m &= F_m \frac{\Delta t g \rho}{\Delta p} \\
\delta T_m &= - \frac{L_i}{c_p}\delta q_m
\end{aligned}
```

But note that we do not add ``\delta q_m`` to the humidity tendency as this is a
phase transition from snow to rain water and so does not increase water vapour ``q``.
We solely use this to calculate the rain water concentration in ``[kg/kg]`` from melting,
and translate it to latent heat.

We calculate the melting of a downward snow flux before [Re-evaporation](@ref).
This is such that melting snow becomes rain water and is subject to reevaporation
within one layer. This is effectively equivalent to allowing sublimation of snow
to water vapour.

Snow can be enabled/disabled with the `snow` keyword argument, and
``T_f`` and ``T_m`` (both in Kelvin) can be passed on too, e.g.

```@example condensation
spectral_grid = SpectralGrid()
large_scale_condensation = ImplicitCondensation(spectral_grid, snow=true, freezing_threshold=263)
```

## Large-scale precipitation

The tendencies ``\delta q`` in units of kg/kg/s are vertically
integrated (top to bottom, the direction of pressure ``p``) to diagnose the
large-scale precipitation ``P`` in units of meters

```math
P = - \int_{top}^{bottom} \frac{\Delta t}{g \rho} \delta q dp 
```

with gravity ``g``, water density ``\rho`` and time step ``\Delta t``.
``P`` is therefore interpreted as the amount of precipitation that falls down
during the time step ``\Delta t`` of the time integration. Note that ``\delta q``
is always negative due to the ``q > q^\star`` condition for saturation,
hence ``P`` is positive only.
It is then accumulated over several time steps, e.g. over the course of an
hour to yield a typical rain rate of mm/h.
The water density is taken as reference density of ``1000~kg/m^3``.

While the above equation to diagnose large-scale precipitation holds, this
would include the total precipitation from both rain water and snow water
(make sure _not_ to add ``\delta q_m`` the snow to rain melt rate).
Given that we already compute the rain water and snow water fluxes ``F_r, F_s``
throughout the column we can also use ``F_{r, N+1}, F_{s, N+1}`` the fluxes
out of the last layer ``N`` and have ``P = \Delta t (F_{r, N+1} + F_{s, N+1})``.

A schematic of the large-scale precipitation parameterization is illustrated below:

![image](https://github.com/user-attachments/assets/4422c0df-7a5d-40ba-9434-da6292bce489)





## References

[^Frierson2006]: Frierson, D. M. W., I. M. Held, and P. Zurita-Gotor, 2006: A Gray-Radiation Aquaplanet Moist GCM. Part I: Static Stability and Eddy Scale. J. Atmos. Sci., 63, 2548-2566, DOI:[10.1175/JAS3753.1](https://doi.org/10.1175/JAS3753.1).
