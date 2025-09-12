# Vertical diffusion

Vertical diffusion in SpeedyWeather.jl is implemented as a Laplacian in
the vertical [Sigma coordinates](@ref) with a diffusion coefficient ``K``
that in general depends on space and time and is flow-aware, meaning
it is recalculated on every time step depending on the vertical stability
of the atmospheric column.

Vertical diffusion can be applied to velocities ``u, v``, temperature ``T``
(done via dry static energy, see below) and specific humidity ``q``.

## Implementations

The following schemes for vertical diffusion are currently implemented

```@example surface_fluxes
using InteractiveUtils # hide
using SpeedyWeather
subtypes(SpeedyWeather.AbstractVerticalDiffusion)
```

You can always set `vertical_diffusion=nothing` which will disable all vertical diffusion.
`BulkRichardsonDiffusion` is explained in the following.

## Laplacian in sigma coordinates

The vertical diffusion of a variable, say ``u``, takes the in
sigma coordinates the form

```math
\frac{\partial u}{\partial t} = \frac{\partial}{\partial \sigma} K
\frac{\partial u}{\partial \sigma}
```

as a tendency to ``u`` with no-flux boundary conditions
``\frac{\partial u}{\partial \sigma} = 0``
at ``\sigma = 0`` (top of the atmosphere) and ``\sigma = 1`` (the surface).
That way the diffusion preserves the integral of the variable ``u`` from
``\sigma = 0`` to ``\sigma = 1``.

```math
\frac{\partial}{\partial t} \int_0^1 u d\sigma  = \int_0^1 \frac{\partial}{\partial \sigma} K
\frac{\partial u}{\partial \sigma} d\sigma = 
K\frac{\partial u}{\partial \sigma} \vert_{\sigma = 1} - K\frac{\partial u}{\partial \sigma} \vert_{\sigma = 0}
= 0
```

Discretising the diffusion operator ``\partial_\sigma K \partial_\sigma`` over ``N``
vertical layers ``k = 1...N`` with ``u_k``, ``K_k`` on those layers at respective
coordinates ``\sigma_k`` that are generally not equally spaced using centred finite differences

```math
\frac{\partial}{\partial \sigma} K \frac{\partial u}{\partial \sigma} \approx
\frac{
    \frac{K_{k+1} + K_k}{2}     \frac{u_{k+1} - u_k    }{\sigma_k+1 - \sigma_k} - 
    \frac{K_{k}   + K_{k-1}}{2} \frac{u_{k}   - u_{k-1}}{\sigma_k   - \sigma_{k-1}}
}{
    \sigma_{k+1/2} - \sigma_{k-1/2}
}
```

We reconstruct ``K`` on the faces ``k+1/2`` with a simple
arithmetic average of ``K_k`` and ``K_{k+1}``. This is necessary for the
multiplication with the gradients which are only available on the faces after the
centred gradients are computed. We then take the gradient again to obtain the final
tendencies again at cell centres ``k``.

## Bulk Richardson-based diffusion coefficient

We calculate the diffusion coefficient ``K`` based on the bulk Richardson number ``Ri``
[^Frierson2006] which is computed as follows

```math
Ri = \frac{gz \left( \Theta_v(z) - \Theta_v(z_N) \right)}{|v(z)|^2 \Theta_v(z_N)}
```

(see [Bulk Richardson-based drag coefficient](@ref) in comparison).
Gravitational acceleration is ``g``, height ``z``, ``\Theta_v`` the virtual
potential temperature where we use the virtual dry static energy
``c_pT_v + gz`` with ``T_v`` the [Virtual temperature](@ref).
The boundary layer height ``h`` (vertical index ``k_h``) is defined as the height
of the lowermost layer where ``Ri_{k_h} > Ri_c`` with ``Ri_c = 1``
the critical Richardson number.

The diffusion coefficient ``K`` is for every layer ``k \geq k_h`` in the boundary
layer calculated depending on the height ``z`` of a layer and its bulk Richardson
number ``Ri``.

```math
K(z) = \begin{cases}
    K_b(z) \quad &\text{for} \quad z \leq f_b h \\
    K_b(f_b h) \frac{z}{f_b h} \left( 1 - \frac{z - f_b h}{(1 - f_b)h} \right)^2
        \quad &\text{for} \quad f_b h < z \leq h \\
\end{cases}
```

with ``f_b = 0.1`` the fraction of the boundary layer height ``h`` above which
the second case guarantees a smooth transition in ``K`` to zero at ``z = h``.
``K_b(z)`` is then defined as

```math
Kb(z) = \begin{cases}
    \kappa u_N \sqrt{C}z \quad &\text{for} \quad Ri_N \leq 0 \\
    \kappa u_N \sqrt{C}z \left( 1 + \frac{Ri}{Ri_c}\frac{\log(z/z_0)}{1 - Ri/Ri_c}\right)^{-1}
        \quad &\text{for} \quad Ri_N > 0 \\
\end{cases}
```

``C`` is the surface drag coefficient as computed in [Bulk Richardson-based drag coefficient](@ref).
The subscript ``N`` denotes the lowermost model layer ``k=N``.

As for the [Bulk Richardson-based drag coefficient](@ref) we also simplify this calculation
here by approximating ``\log(z/z_0) \approx \log(Z/z_0)`` with the height Z of the lowermost layer
given resolution and a reference surface temperature, for more details see description in that section.

## Vertical diffusion tendencies

The vertical diffusion is then computed as a tendency for ``u, v, q`` and temperature ``T`` via the
dry static energy ``SE = c_p T + gz``, i.e.

```math
\frac{\partial T}{\partial t} = \frac{1}{c_p}\frac{\partial}{\partial \sigma} K
\frac{\partial SE}{\partial \sigma} 
```

where we just fold the heat capacity ``c_p`` into the diffusion coefficient ``K \to K/c_p``.
The other variables are diffused straight-forwardly as
``\partial_t u = \partial_\sigma K \partial_\sigma u``, etc.

## References

[^Frierson2006]: Frierson, D. M. W., I. M. Held, and P. Zurita-Gotor, 2006: A Gray-Radiation Aquaplanet Moist GCM. Part I: Static Stability and Eddy Scale. J. Atmos. Sci., 63, 2548-2566. DOI: [10.1175/JAS3753.1](https://doi.org/10.1175/JAS3753.1).