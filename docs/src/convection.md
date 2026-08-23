# Convection

Convection is the atmospheric process of rising motion because of positively
buoyant air parcels compared to its surroundings. In hydrostatic models
like the primitive equation model in SpeedyWeather.jl convection has to
be parameterized as the vertical velocity is not a prognostic variable
that depends on vertical stability but rather diagnosed to satisfy
horizontal divergence. Convection can be shallow and non-precipitating
denoting that buoyant air masses can rise but do not reach saturation
until they reach a level of zero buoyancy. But convection can also be
_deep_ denoting that saturation has been reached during ascent whereby
the latent heat release from condensation provides additional energy
for further ascent. Deep convection is therefore usually also
precipitating as the condensed humidity forms cloud droplets that
eventually fall down as convective precipitation.
See also [Large-scale condensation](@ref) in comparison.

## Convection implementations

Currently implemented are
```@example convection
using InteractiveUtils # hide
using SpeedyWeather
subtypes(SpeedyWeather.AbstractConvection)
```
which are described in the following.

## [Simplified Betts-Miller convection](@id BettsMiller)

We follow the simplification of the Betts-Miller convection scheme
[^Betts1986][^BettsMiller1986] as studied by Frierson, 2007 [^Frierson2007].
The central idea of this scheme is to represent the effect of
convection as an adjustment towards a (pseudo-) moist adiabat
reference profile and its associated humidity profile. Meaning that
conceptually for every vertical column in the atmosphere we

1. Diagnose the vertical temperature and humidity profile of the environment relative to the adiabat up to the level of zero buoyancy.
2. Decide whether convection should take place and whether it is deep (precipitating) or shallow (non-precipitating).
3. Relax temperature and humidity towards (corrected) profiles from 1.

## Reference profiles

The dry adiabat is

```math
T = T_0 \left(\frac{p}{p_0}\right)^\frac{R}{c_p}
```

The temperature ``T`` of an air parcel at pressure ``p`` is determined by the
temperature ``T_0`` it had at pressure ``p_0`` (this can be surface but it does not have to be)
and the gas constant for dry air ``R = 287.04 J/K/kg`` and the heat capacity ``c_p = 1004.64 J/K/kg``.
The pseudo adiabat follows the dry adiabat until saturation is reached
(the lifting condensation level, often abbreviated to LCL), ``q > q^\star``.
Then it follows the pseudoadiabatic lapse rate

```math
\Gamma = -\frac{dT}{dz} = \frac{g}{c_p}\left[
    \frac{  1 + \frac{q^\star L_v  }{(1-q^\star)^2     R_d T_v}}{
            1 + \frac{q^\star L_v^2}{(1-q^\star)^2 c_p R_v T^2}}\right]
```

with gravity ``g``, heat capacity ``c_p``, the saturation specific humidity of the parcel ``q^\star``
(which is its specific humidity given that it has already reached saturation),
latent heat of vaporization ``L_v``, dry gas constant ``R_d``, water vapor
gas constant ``R_v``, and [Virtual temperature](@ref) ``T_v``.
Starting with a temperature ``T`` and humidity ``q = q^\star`` at the lifting condensation level
temperature aloft changes with $dT = -\frac{d\Phi}{c_p}(...)$ between two layers separated
``d\Phi`` in geopotential ``\Phi`` apart. On that new layer, ``q^\star`` is recalculated
as well as the virtual temperature ``T_v = T(1 + \mu q^\star)``. ``\mu`` is derived from
the ratio of dry to vapor gas constants see [Virtual temperature](@ref).
Note that the pseudoadiabatic ascent is independent of the environmental temperature
and humidity and function of temperature and humidity of the parcel only
(although that one starts with surface temperature and humidity from the environment).
Solely the level of zero buoyancy is determined by comparing the parcel's virtual
temperature ``T_v`` to the virtual temperature of the environment ``T_{v,e}`` at that level.
Level of zero buoyancy is reached when ``T_v = T_{v,e}`` but continues for ``T_v > T_{v,e}``
which means that the parcel is still buoyant. Note that the virtual temperature includes
the effect that humidity has on its density.

The (absolute) temperature a lifted parcel has during ascent (following its pseudoadiabat, dry and/
or moist, until reaching the level of zero buoyancy) is then taken as the reference temperature
profile ``T_{ref}`` that the Betts-Miller convective parameterization relaxes towards as a first guess
(with a following adjustment as discussed below). The humidity profile is taken as
``q_{ref} = RH_{SBM}\,q^\star(T_{ref})`` with a parameter ``RH_{SBM}`` (default ``RH_{SBM} = 0.7``)
of the scheme (Simplified Betts-Miller, SBM) that determines a constant relative humidity of
the reference profile, applied to the saturation humidity ``q^\star`` evaluated at the reference
temperature ``T_{ref}`` of that layer.

An illustration of adiabat reference line is below

![image](https://github.com/user-attachments/assets/ea8d8a0b-20c5-47cc-88ff-47112eeeb200)

## Entrainment

A real convective plume mixes with the environmental air around it as it rises, which dilutes
the parcel's temperature and moisture excess, reducing its buoyancy and moisture content and
generally shallowing convection compared to the undiluted ascent described above. The scheme
optionally applies this dilution at every level ``k`` during the ascent, after the dry or moist
adiabatic step:

```math
\begin{aligned}
T_p &\leftarrow (1 - \varepsilon_k) T_p + \varepsilon_k T_{e,k} \\
q_p &\leftarrow (1 - \varepsilon_k) q_p + \varepsilon_k q_{e,k}
\end{aligned}
```

with parcel temperature/humidity ``T_p, q_p``, environmental temperature/humidity at that level
``T_{e,k}, q_{e,k}``, and an entrainment rate ``\varepsilon_k \in [0, 1]`` that can vary by level.
Because ``\varepsilon_k`` is a mixing fraction *per model level* rather than a rate per unit
height, the resulting profile depends on the vertical resolution (`nlayers`): the same
`AbstractEntrainment` profile mixes in more environmental air over a coarser column, since each
level spans a larger part of the atmosphere. In the moist branch of the ascent, mixing happens
after condensation at that level, so the saturation adjustment (and thus the LZB and the
associated ``T_{ref}, q_{ref}`` at that level) already reflects the diluted parcel; the latent
heat implied by the change in ``q_p`` from mixing itself is not separately accounted for.

Available profiles, selected via the `entrainment` option of [`BettsMillerConvection`](@ref) or
[`BettsMillerDryConvection`](@ref):

```@example convection
subtypes(SpeedyWeather.AbstractEntrainment)
```

- `NoEntrainment` (the default): ``\varepsilon_k = 0`` everywhere, reproducing the scheme exactly
  as described above with no dilution.
- `LinearEntrainment`: ``\varepsilon_k`` ramps linearly from `surface_entrainment` at the surface
  (``\sigma = 1``) to zero at `σ_entrainment`, and is zero above that level.
- `ConstantEntrainment`: a single `entrainment_rate` applied at every level.

Since `NoEntrainment` is the default, entrainment does not change any model's behaviour unless
explicitly enabled, e.g. `BettsMillerConvection(spectral_grid; entrainment = LinearEntrainment(spectral_grid))`.

## First-guess relaxation

With the [Reference profiles](@ref) ``T_{ref}, q_{ref}`` obtained, we relax the actual
environmental temperature ``T`` and specific humidity ``q`` in the column

```math
\begin{aligned}
\delta q &= - \frac{q - q_{ref}}{\tau_{SBM}} \\
\delta T &= - \frac{T - T_{ref}}{\tau_{SBM}}
\end{aligned}
```

with the second parameter of the parameterization, the time scale ``\tau_{SBM}``.
Note that because this is a first-guess relaxation, these tendencies are not actually
the resulting tendencies from this scheme. Those will be calculated in [Corrected relaxation](@ref).

Note that above the level of zero buoyancy no relaxation takes place ``\delta T = \delta q = 0``,
or, equivalently ``T = T_{ref}``, ``q = q_{ref}`` there.
Vertically integration from the level of zero buoyancy ``p_{LZB}`` to the surface ``p_0`` (so that
``dp > 0``, as ``p_{LZB} < p_0``) yields

```math
\begin{aligned}
P_q &= - \int_{p_{LZB}}^{p_0} \delta q \frac{dp}{g} \\
P_T &= \int_{p_{LZB}}^{p_0} \frac{c_p}{L_v} \delta T \frac{dp}{g}
\end{aligned}
```

``P_q`` is the precipitation in units of ``kg / m^2 / s`` due to drying (as a consequence of
the humidity tendency) and ``P_T`` is the precipitation in the same units due to warming
(as resulting from temperature tendencies). Note that they are the vertically difference
between current profiles and the references profiles, so if ``P_q > 0`` this would mean
that a convective adjustment to ``q_{ref}`` would release humidity from the column
through condensation, but ``P_q`` can also be negative. Consequently similar for ``P_T``.


## Convective criteria

We now distinguish three cases

1. Deep convection when ``P_T > 0`` and ``P_q > 0``
2. Shallow convection when ``P_T > 0`` and ``P_q \le 0``
3. No convection for ``P_T \le 0``.

Note that to evaluate these cases it is not necessary to divide by ``\tau_{SBM}`` in
the first-guess relaxation, neither are the ``1/g`` and ``\tfrac{c_p}{g L_v}``
necessary to multiply during the vertical integration as all are positive constants.
While this changes the units of ``P_T, P_q`` they can be reused in the following.

This is illustrated in the following

![image](https://github.com/user-attachments/assets/ec5dd0b8-b350-4412-920d-80c4ea9dbb4c)

## Deep convection

Following Frierson, 2007 [^Frierson2007] in order to conserve enthalpy we correct
the reference profile for temperature ``T_{ref} \to T_{ref, 2}`` so that ``P_T = P_q``.

```math
T_{ref, 2} = T_{ref} + \frac{1}{\Delta p c_p} \int_{p_0}^{p_{LZB}} \left[ c_p (T - T_{ref}) + L_v (q - q_{ref}) \right] dp
```

``\Delta p`` is the pressure difference ``p_{LZB} - p_0``.
The terms inside the integral are rearranged compared to Frierson, 2007 to show
that the vertical integral in [First-guess relaxation](@ref) really only has to be computed once.

## Shallow convection

In the following we describe the "qref" scheme from Frierson, 2007 which corrects
reference profiles for both temperature and humidity to guarantee that ``P_q = 0``,
i.e. no precipitation during convection. In that sense, shallow convection is
non-precipitating. Although shallow convection is supposed to be shallow
we do not change the height of the convection and keep using the ``p_{LZB}``
determined during the calculation of the [Reference profiles](@ref).

```math
\begin{aligned}
\Delta q &= \int_{p_{LZB}}^{p_0} (q - q_{ref}) dp \\
Q_{ref}  &= \int_{p_{LZB}}^{p_0} -q_{ref}~dp \\
f_q      &= 1 - \frac{\Delta q}{Q_{ref}} \\
q_{ref, 2} &= f_q q_{ref} \\
\Delta T &= -\frac{1}{\Delta p} \int_{p_{LZB}}^{p_0} (T - T_{ref}) dp \\
T_{ref,2} &= T_{ref} - \Delta T
\end{aligned}
```

This scheme is non-precipitating (``P_q = 0``) as derived below. ``P_q`` follows from the vertical integral of the moisture tendency ``\delta q``, which in this scheme is a relaxation towards ``q_{\text{ref,2}}`` (see [Corrected relaxation](@ref)).

```math
P_q = -\int_{p_{LZB}}^{p_0} \frac{\delta q}{g} \, dp = \int_{p_{LZB}}^{p_0} \frac{q - q_{\text{ref,2}}}{g \tau_{SBM}} dp
```

Inserting from above yields

```math
P_q = \int_{p_{LZB}}^{p_0} \frac{q - \left(1 - \frac{\Delta q}{Q_{\text{ref}}} \right) q_{\text{ref}} }{g \tau_{SBM}} dp = \frac{1}{g \tau_{SBM}} \int_{p_{LZB}}^{p_0} \left( q - q_{\text{ref}} + \frac{\Delta q}{Q_{\text{ref}}} q_{\text{ref}} \right) dp
```

The integral becomes
```math
\Delta q + \frac{\Delta q}{Q_{\text{ref}}} \int_{p_{LZB}}^{p_0}q_{\text{ref}} dp = \Delta q + \frac{\Delta q}{Q_{\text{ref}}} (- Q_{\text{ref}}) = 0
```

So this scheme is indeed non-precipitating, i.e. $P_q = 0$.

## Corrected relaxation

After the reference profiles have been corrected in [Deep convection](@ref)
and [Shallow convection](@ref) we actually calculate tendencies from

```math
\begin{aligned}
\delta q &= - \frac{q - q_{ref, 2}}{\tau_{SBM}} \\
\delta T &= - \frac{T - T_{ref, 2}}{\tau_{SBM}}
\end{aligned}
```

with ``\tau_{SBM} = 4h`` as default.

## Convective precipitation

The convective precipitation ``P`` results then from the vertical integration of the
``\delta q`` tendencies,
similar to [Large-scale precipitation](@ref).

```math
P = -\int_{p_{LZB}}^{p_0} \frac{\Delta t}{g \rho} \delta q\,dp
```

Only layers with ``\delta q < 0`` (drying, i.e. condensation) contribute to this integral;
layers that moisten (``\delta q > 0``) are clamped to zero rather than allowed to offset the
total, since the scheme does not model reevaporation of falling convective rain. The result is
additionally multiplied by a deep-convection indicator so that ``P = 0`` identically for shallow
convection: with the ``\delta q \ge 0``-clamp in place, the qref correction's exact ``P_q = 0``
result (derived in [Shallow convection](@ref)) is not automatically reproduced layer-by-layer, so
the indicator is what actually enforces no precipitation in the shallow case, not merely a
redundant restatement of it.

In the shallow convection case ``P=0`` due to the correction (and the indicator above) even though
in the first guess relaxation ``P<0`` was possible, but for deep convection ``P>0`` by definition.

### Convective snow

If the temperature of the lowermost model layer is below a `freezing_threshold` (default
``273.15\,K``), all of ``P`` falls as snow instead of rain for that column and time step,
controlled by the `snow` option (`true` by default). This is a single check on the layer the
parcel starts its ascent from, unlike [Large-scale condensation](@ref)'s falling-flux melt
cascade through intermediate layers: convective precipitation is generated and deposited within
the same time step rather than fluxed layer by layer, so there is no intermediate layer for it
to fall through and potentially melt in. The phase swap conserves total precipitation mass
(rain + snow is unchanged, only the partitioning between the two is affected), and shallow
convection produces neither since ``P = 0`` there regardless of temperature. With `snow = false`
all convective precipitation is always rain, recovering the previous behaviour.

## Dry convection

In the primitive equation model with humidity the [Betts-Miller convection scheme](@ref BettsMiller)
as described above is defined. Without humidity, a dry version reduces to the
[Shallow convection](@ref) case. The two different shallow convection schemes in
Frierson 2007[^Frierson2007], the "shallower" shallow convection scheme and the "qref"
(as implemented here in [Shallow convection](@ref)) in that case also reduce to
the same formulation. The dry Betts-Miller convection scheme is the default
in the primitive equation model without humidity.


## References

[^Betts1986]: Betts, A. K., 1986: A new convective adjustment scheme. Part I: Observational and theoretical basis. Quart. J. Roy. Meteor. Soc.,112, 677-691. DOI: [10.1002/qj.49711247307](https://doi.org/10.1002/qj.49711247307)

[^BettsMiller1986]: Betts, A. K. and M. J. Miller, 1986: A new convective adjustment scheme. Part II: Single column tests using GATE wave, BOMEX, ATEX and Arctic air-mass data sets. Quart. J. Roy. Meteor. Soc.,112, 693-709. DOI: [10.1002/qj.49711247308](https://doi.org/10.1002/qj.49711247308)

[^Frierson2007]: Frierson, D. M. W., 2007: The Dynamics of Idealized Convection Schemes and Their Effect on the Zonally Averaged Tropical Circulation. J. Atmos. Sci., 64, 1959-1976. DOI:[10.1175/JAS3935.1](https://doi.org/10.1175/JAS3935.1)
