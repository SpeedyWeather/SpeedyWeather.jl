# Large-scale condensation

Large-scale condensation in an atmospheric general circulation represents the
micro-physical that kicks in when an air parcel reaches saturation.
Subsequently, the water vapour inside it condenses, forms droplets around
condensation nuclei, which grow, become heavy and eventually fall out
as precipitation. This process is never actually representable at the resolution
of global (or even regional) atmospheric models as typical cloud droplets
have a size of micrometers. Atmospheric models therefore rely on large-scale
quantities such as specific humidity, pressure and temperature within a
given grid cell, even though there might be considerably variability of
these quantities within a grid cell if the resolution was higher.

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
undersaturated. Ideally, one would want to condense towards the new saturation humidity
``q^\star(T_{i+1})`` so that condensation draws the relative humidity back down to 100% not below it. 
Taylor expansion at ``i`` with ``\delta T = T_{i+1} - T_i`` (and ``\delta q`` similarly)
to first order yields

```math
q_{i+1} - q_i = q^\star(T_{i+1}) - q_i = q^\star(T_i) + (T_{i+1} - T_i)
\frac{\partial q^\star}{\partial T} (T_i) + O(\delta T^2) - q_i
```

Now we make a linear approximation to the derivative and drop the ``O(\delta T^2)`` term.
Inserting the (explicit) latent heat release yields

```math
\delta q = q^\star(T_i) + -\frac{L_v}{c_p} \delta q \frac{\partial q^\star}{\partial T} (T_i) - q_i
```

And solving for ``\delta q`` yields

```math
\left[ 1 + \frac{L_v}{c_p} \frac{\partial q^\star}{\partial T} (T^i) \right] \delta q = q^\star(T_i) - q_i
```

meaning that the implicit immediate condensation can be formulated as

```math
\begin{aligned}
q_{i+1} - q_i &= \frac{q^\star(T_i) - q_i}{1 + \frac{L_v}{c_p} \frac{\partial q^\star}{\partial T}(T_i)} \\
T_{i+1} - T_i &= -\frac{L_v}{c_p}( q_{i+1} - q_i )
\end{aligned}
```

With Euler forward time stepping this is great, but with our [leapfrog timestepping + RAW filter](@id leapfrog)
this is very dispersive (see #445) although the implicit formulation is already much better.
We therefore introduce a time step ``\Delta t_c`` which makes the implicit condensation not immediate
anymore but over several time steps ``\Delta t`` of the leapfrogging.

```math
\begin{aligned}
\frac{q_{i+1} - q_i}{\Delta t} &= \frac{q^\star(T_i) - q_i}{ \Delta t_c
\left( 1 + \frac{L_v}{c_p} \frac{\partial q^\star}{\partial T}(T_i) \right)} \\
\frac{T_{i+1} - T_i}{\Delta t} &= -\frac{L_v}{c_p}( \frac{q_{i+1} - q_i}{\Delta t} )
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