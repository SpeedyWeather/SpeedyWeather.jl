# Vertical diffusion

Vertical diffusion represents the turbulent mixing in the vertical of momentum ``u, v``,
temperature ``T`` (via the dry static energy, see below) and [specific humidity](@ref "Humidity equation") ``q``.
In the boundary layer it spreads the surface fluxes of momentum, heat and moisture
over the depth of the boundary layer; above, it mixes layers with strong wind shear and
removes superadiabatic lapse rates. The diffusivity ``K`` depends on space and time and is
flow-aware: it is recalculated on every time step from the stability of the atmospheric column.

## Implementations

The following schemes for vertical diffusion are currently implemented

```@example vertical_diffusion
using InteractiveUtils # hide
using SpeedyWeather
subtypes(SpeedyWeather.AbstractVerticalDiffusion)
```

`BulkRichardsonDiffusion` is the default in both `PrimitiveDryModel` and `PrimitiveWetModel`.
Its options are passed on as keyword arguments, and `vertical_diffusion = nothing` disables it

```@example vertical_diffusion
spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
vertical_diffusion = BulkRichardsonDiffusion(spectral_grid, background_diffusivity = 0.5)
model = PrimitiveWetModel(spectral_grid; vertical_diffusion)
model_without = PrimitiveWetModel(spectral_grid; vertical_diffusion = nothing)
nothing # hide
```

`BulkRichardsonDiffusion` is explained in the following. It uses the surface drag and wind of the
boundary layer (`model.boundary_layer`, the default) and is called after the surface fluxes.
In a model without a boundary layer only its local and background diffusivity act. See [Parameterizations](@ref) for the
general parameterization interface these implement.

## Diffusion as a flux divergence

The vertical diffusion of a variable ``X`` with diffusivity ``K`` [m²/s]
takes the form

```math
\frac{\partial X}{\partial t} = \frac{1}{\rho} \frac{\partial}{\partial z} \left( \rho K \frac{\partial X}{\partial z} \right)
= -g \frac{\partial F}{\partial p}, \qquad F = \rho K \frac{\partial X}{\partial z}
```

with air density ``\rho``, height ``z``, pressure ``p`` and gravity ``g`` using the hydrostatic
balance ``\partial p = -\rho g \partial z``. ``F`` is the (negative) turbulent flux of ``X``.
We discretise this over the ``N`` layers ``k = 1 \dots N`` (top to surface, see
[Vertical coordinates](@ref vertical_coordinates_physics)) with ``K_{k+1/2}`` defined on the
interface ``k+1/2`` between layer ``k`` and ``k+1``

```math
\frac{\partial X_k}{\partial t} = \frac{g}{\Delta p_k} \left[ a_{k+1/2} (X_{k+1} - X_k) - a_{k-1/2} (X_k - X_{k-1}) \right],
\qquad a_{k+1/2} = \frac{\rho_{k+1/2} K_{k+1/2}}{z_k - z_{k+1}}
```

with ``\Delta p_k`` the pressure thickness of layer ``k``, the conductance ``a`` [kg/m²/s],
the density ``\rho_{k+1/2} = p_{k+1/2} / (R_d T_{v, k+1/2})`` from the pressure and the
[Virtual temperature](@ref) on the interface, and the distance between layers
``z_k - z_{k+1} = (\Phi_k - \Phi_{k+1})/g`` from the [Geopotential](@ref) ``\Phi``.
There is no flux through the top of the atmosphere, ``a_{1/2} = 0``, and no diffusive flux
through the surface, ``a_{N+1/2} = 0``, as the surface fluxes are computed separately.
This formulation works for any [Vertical coordinates](@ref vertical_coordinates_physics) and
conserves the mass-weighted column integral ``\sum_k X_k \Delta p_k`` exactly.

Temperature is diffused through the dry static energy ``s = c_p T + \Phi``, holding the geopotential
fixed, so that ``X = T + \Phi/c_p`` and the temperature tendency is the tendency of ``X``.
An adiabatic column with constant dry static energy is therefore not changed by the diffusion.

## Implicit time integration

The diffusivity near the surface, where layers can be less than 100 m thick, can easily
exceed the stability limit ``\Delta z^2 / 2K`` of an explicit time step. The diffusion is
therefore integrated implicitly for every column, with ``K`` evaluated
at the previous time step [KalnayKanamitsu1988](@citep). Writing the
diffusion operator above as ``L``, the tendency ``y = (X^+ - X^-)/\Delta t`` from the previous
value ``X^-`` to the new one ``X^+`` is the flux divergence of a state ``\hat{X}``

```math
y = L \hat{X} + S, \qquad \hat{X} = \gamma X^+ + (1 - \gamma) X^- = X^- + \gamma \Delta t \, y
```

with ``\Delta t`` the time step of the parameterizations (``2\Delta t`` for leapfrog),
``S`` the tendencies of the parameterizations called before the vertical diffusion (see below)
and the implicitness factor ``\gamma`` (option `implicitness`): ``\gamma = 0`` would be an explicit
step, ``\gamma = 1/2`` Crank-Nicolson, ``\gamma = 1`` is backward Euler and ``\gamma > 1`` is called
over-implicit. The tendency ``y`` is the solution of the tridiagonal system

```math
(1 - \gamma \Delta t L) y = L X^- + S
```

which is solved with the Thomas algorithm, a forward elimination and a back substitution over the
layers. This is the backward Euler solve with the time step ``\gamma \Delta t``, while the model
still steps with ``\Delta t``. The system is diagonally dominant so that the solution is
unconditionally stable for a given diffusivity and conserves the column integral, and for
``\gamma \geq 1`` the new value ``X^+ = (1 - 1/\gamma) X^- + \hat{X}/\gamma`` lies between ``X^-`` and
``\hat{X}``, so that the diffusion alone does not create new extrema.

The default is ``\gamma = 1.5``, the value of the IFS, see section 3.5 "Solution of the EDMF equations",
equations (3.65) and (3.66), of [IFSCycle49r1](@citet): "In the model, ``\alpha = 1.5``, to avoid non-linear
instability from the K-coefficients", where ``\alpha`` is our ``\gamma``. The reason is that the diffusivity is
not constant but that of the previous time step and grows with the diffused state itself, e.g.
``K \propto u_* \propto |\mathbf{u}|`` in the boundary layer. Backward Euler with such an explicit
coefficient is only stable for a weak dependence [KalnayKanamitsu1988](@citep): for a damping rate
``\propto X^P`` and a large diffusivity a perturbation of a forced equilibrium is multiplied by
``(\gamma - 1 - P)/\gamma`` every step, which is ``-P`` for backward Euler, an oscillation between two
states that grows for ``P > 1``, and stable for ``P < 2\gamma - 1``. With ``\gamma = 1`` this model at
T512 with 24 layers developed such oscillations (also known as fibrillations) in single columns with
60 to 75 m/s wind and diffusivities of several hundred m²/s in the lowermost layers, with a period of
four time steps (two steps of leapfrog), which grew until the simulation failed. For ``\gamma > 1`` the tendencies
``S`` enter ``\hat{X}`` with the weight ``\gamma`` too, which is why the IFS uses ``\gamma = 1`` for
tracers, where the combination with a sink can otherwise create negative values (section 3.7 of
[IFSCycle49r1](@citet)). Here all diffused variables use the same factor. Backward Euler is

```julia
vertical_diffusion = BulkRichardsonDiffusion(spectral_grid, implicitness = 1)
```

The source ``S`` are the tendencies of the parameterizations called before the vertical diffusion,
by default the surface fluxes into the lowermost layer. They are applied first and the result is
diffused, so that the surface fluxes are spread over the boundary layer within the time step instead of
heating, moistening or decelerating only the lowermost layer. This is a split version of the
coupled implicit solve of surface fluxes and vertical diffusion in [Frierson2006](@citet).

## Diffusivity

The diffusivity on every interface is the maximum of three parts, capped at ``K_{max}`` (default 1000 m²/s)

```math
K = \min(\max(K_{profile}, K_{local}, K_{background}), K_{max})
```

and is computed separately for momentum ``K_m`` and for heat and humidity ``K_h``, which
differ only in ``K_{local}``. The implicit solver only needs the conductances
``a = \rho K / \Delta z`` [kg/m²/s] from above, so these are stored, not the diffusivities, as
`vertical_conductance_momentum` and `vertical_conductance_heat` in
`simulation.variables.parameterizations`, next to the boundary layer height ``h`` as
`boundary_layer_height`. The value belongs to the interface below each layer and is zero for
the lowermost layer (no flux through the surface). The diffusivity is
``K = a \Delta z / \rho`` with the density and layer distance as defined above.
The boundary layer height is written to the output with
`add!(model, OutputVariables.BoundaryLayerHeight())` (`blh` in m), which is also part of the
output group `OutputVariables.BoundaryLayer()`.

### Boundary layer profile

Following [Frierson2006](@citet) we define the boundary layer height ``h`` with the bulk Richardson number

```math
Ri(z) = \frac{gz \left[ s_v(z) - s_v(z_N) \right]}{c_p T_v(z_N) |v(z)|^2}
```

between the lowermost layer ``N`` at height ``z_N`` and the layer at height ``z`` above the surface,
using the virtual dry static energy ``s_v = c_p T_v + \Phi`` and the wind speed ``|v|``, at least
1 m/s [HoltslagBoville1993](@citep). ``h`` is the height where ``Ri(z)`` first reaches the critical Richardson
number ``Ri_c = 1`` (option `critical_Richardson`) going upwards from the lowermost layer with
``Ri(z_N) = 0``, interpolated linearly between the layers. ``Ri_c = 1`` is the value of
[Frierson2006](@citet) and at the upper end of the literature values of 0.25–1 [HoltslagBoville1993](@citep).
For a given stratification the height grows like ``\sqrt{Ri_c}`` times the wind speed: a much larger value such
as ``Ri_c = 10`` gives tropical boundary layers of 1–2 km, but boundary layers that reach the tropopause under
the mid-latitude jets. At T512 with 24 layers that mixing through the troposphere warmed the lowermost
layer by 4–6 K, reversed the sensible heat flux over the oceans and reduced the precipitation by 30 %
within half a year. Within the boundary layer

```math
K_{profile}(z) = \begin{cases}
    K_b(z) \quad &\text{for} \quad z \leq f_b h \\
    K_b(f_b h) \frac{z}{f_b h} \left[ 1 - \frac{z - f_b h}{(1 - f_b)h} \right]^2
        \quad &\text{for} \quad f_b h < z \leq h \\
\end{cases}
```

and zero above. ``f_b = 0.1`` is the surface layer fraction of the boundary layer, above which
``K`` decreases to zero at ``z = h``. ``K_b(z)`` is

```math
K_b(z) = \kappa u_\star z \left[ 1 + \frac{Ri_a}{Ri_{cut}} \frac{\log(z/z_0)}{1 - Ri_a/Ri_{cut}}\right]^{-1}
```

with the von Kármán constant ``\kappa = 0.4`` and the friction velocity ``u_\star = \sqrt{C} V_s``
from the drag coefficient for momentum ``C`` of the boundary layer and the surface wind speed ``V_s``,
see [Bulk Richardson-based drag coefficient](@ref). The stability correction in brackets uses the
bulk Richardson number ``Ri_a`` between the surface and the lowermost layer (computed as for the drag,
from the soil and sea surface temperatures, weighted with the land fraction)
clamped to ``0 \leq Ri_a \leq Ri_{cut}``, so that it is 1 for unstable and neutral conditions and 0
for ``Ri_a \geq Ri_{cut}``. Frierson et al. (2006) use ``Ri_c`` here. We use a separate cutoff
``Ri_{cut} = 10`` (option `Richardson_cutoff`), as for the drag, so that the profile only vanishes
where buoyancy dominates the shear by an order of magnitude. ``Ri_c`` keeps its conventional meaning as
the threshold that defines the boundary layer height. ``z_0`` is the roughness length for momentum.

### Local diffusivity

Above the boundary layer, and within it where the profile is weaker,
the diffusivity follows the local first-order closure of [Louis1979](@citet)

```math
K_{local} = l^2 S f(Ri), \qquad Ri = \frac{N^2}{S^2}, \qquad
S = \left| \frac{\partial \mathbf{v}}{\partial z} \right|, \qquad
N^2 = \frac{g}{c_p T_v} \frac{\partial s_v}{\partial z}
```

with the gradient Richardson number ``Ri`` on the interfaces from the wind shear ``S`` and the
squared buoyancy frequency ``N^2``. The mixing length is ``1/l = 1/\kappa z + 1/\lambda``
with the asymptotic mixing length ``\lambda = 150`` m (option `mixing_length`) below 1 km
that decays towards ``\lambda = 30`` m (option `free_mixing_length`) in the free atmosphere as
``\lambda(z) = 30 + 120 \exp(1 - z/1000)`` above 1 km [HoltslagBoville1993](@citep).
For stable stratification ``Ri > 0`` the long-tail stability functions of
[LouisTiedtkeGeleyn1982](@citet)

```math
f_m(Ri) = \frac{1}{1 + 2bRi/\sqrt{1 + dRi}}, \qquad f_h(Ri) = \frac{1}{1 + 3bRi\sqrt{1 + dRi}}, \qquad b = d = 5
```

are used for momentum ``f_m`` and heat ``f_h``. They decay slowly with ``Ri`` so that some mixing
remains also in strong inversions, which avoids that the surface decouples from the atmosphere
in polar night [Sandu2013](@citep). For unstable stratification ``Ri < 0`` (``N^2 < 0``, superadiabatic)
[HoltslagBoville1993](@citet) use ``f = \sqrt{1 - 18 Ri}`` for both, i.e.
``K_{local} = l^2 \sqrt{S^2 - 18N^2}``, which also mixes superadiabatic layers without wind shear.
With ``\lambda = 30`` m in the free atmosphere, however, a superadiabatic layer there is only
mixed on a time scale of about a day. We therefore split this into a shear part and a buoyancy part

```math
K_{local} = l^2 S + l_c^2 \left( \sqrt{S^2 - 18N^2} - S \right) \quad \text{for} \quad N^2 < 0
```

where the mixing length ``l_c`` of the buoyancy part uses an asymptotic mixing length of at least
150 m (option `convective_mixing_length`) at all heights. This is continuous at ``N^2 = 0``
and mixes superadiabatic layers within a few time steps also in the free troposphere.
`convective_mixing_length = 0` recovers the formulation of [HoltslagBoville1993](@citet).

### Background diffusivity

A background diffusivity ``K_{background} = K_0 \exp(-z/H)`` with ``K_0 = 0.1`` m²/s
(option `background_diffusivity`) and ``H = 1000`` m (option `background_height`) keeps a
minimum of mixing near the surface for any stratification and wind.
