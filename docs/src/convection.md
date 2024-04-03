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

## Simplified Betts-Miller convective adjustment

We follow the simplification of the Betts-Miller convection scheme
[^Betts1986][^BettsMiller1986] as studied by Frierson, 2007 [^Frierson2007].
The central idea of this scheme is to represent the effect of 
convection as an adjustment towards a (pseudo-) moist adiabat 
reference profile and its associated humidity profile. Meaning that
conceptually for every vertical column in the atmosphere we

1. Diagnose the vertical temperature and humidity profile of the environment relative to the surface and calculate the adiabat to the level of zero buoyancy.
2. Decide whether convection should take place and whether it is deep (precipitating) or shallow (non-precipitating).
3. Relax temperature and humidity towards (adjusted) profiles from 1.

## Reference profiles

The dry adiabat is

```math
T = T_0 (\frac{p}{p_0})^\frac{R}{c_p}
```

The temperature ``T`` of an air parcel at pressure ``p`` is determined by the
temperature ``T_0`` it had at pressure ``p_0`` (this can be surface but it does not have to be)
and the gas constant for dry air ``R = 287.04 J/K/kg`` and the heat capacity ``c_p = 1004.64 J/K/kg``.
The pseudo adiabat follows the dry adiabat until saturation is reached
(the lifting condensation level, often abbreviated to LCL), ``q > q^\star``.
Then it follows the pseudoadiabatic lapse rate

```math
\Gamma = -\frac{dT}{dz} = \frac{g}{c_p}\left(
    \frac{  1 + \frac{q^\star L_v  }{(1-q^\star)^2     R_d T_v}}{
            1 + \frac{q^\star L_v^2}{(1-q^\star)^2 c_p R_v T^2}}\right)
```

with gravity ``g``, heat capacity ``c_p``, the saturation specific humidity of the parcel ``q^\star``
(which is its specific humidity given that it has already reached saturation),
latent heat of vaporization ``L_v``, dry gas constant ``R_d``, water vapour
gas constant ``R_v``, and [Virtual temperature](@ref) ``T_v``.
Starting with a temperature ``T`` and humidity ``q = q^\star`` at the lifting condensation level 
temperature aloft changes with $dT = -\frac{d\Phi}{c_p}(...)$ between two layers separated
``d\Phi`` in geopotential ``\Phi`` apart. On that new layer, ``q^\star`` is recalculated
as well as the virtual temperature ``T_v = T(1 + \mu q^\star)``. ``\mu`` is derived from
the ratio of dry to vapour gas constants see [Virtual temperature](@ref).
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
``q_{ref} = RH_{SBM}T_{ref}`` with a parameter ``RH_{SBM}`` (default ``RH_{SBM} = 0.7``)
of the scheme (Simplified Betts-Miller, SBM) that determines a constant relative humidity of
the reference profile.

## First-guess relaxation

With the [Reference profiles](@ref) ``T_{ref}, q_{ref}`` obtained, we relax the actual
environmental temperature ``T`` and specific humidity ``q`` in the column

```math
\begin{aligned}
\delta q &= - \frac{q - q_{ref}}{\tau_{SBM}} \\
\detla T &= - \frac{T - T_{ref}}{\tau_{SBM}}
\end{aligned}
```

with the second parameter of the parameterization, the time scale ``\tau_{SBM}``.
Note that above the level of zero buoyancy no relaxation takes place ``\delta T = \delta q = 0``,
or, equivalently ``T = T_{ref}``, ``q = q_{ref}`` there.
Vertically integration from surface ``p_0`` to level of zero buoyancy in 
pressure coordinates ``p_{LZB}`` yields

```math
\begin{aligned}
P_q &= - \int_{p_0}^p_{LZB} \delta q \frac{dp}{g} \\
P_T &= \int_{p_0}^p_{LZB} \frac{c_p}{L_v} \delta T \frac{dp}{g}
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
2. Shallow convection when ``P_T > 0`` and ``P_q <= 0``
3. No convection for ``P_T <= 0``.

Note that to evaluate these cases it is not necessary to divide by ``\tau_{SBM}`` in
the first-guess relaxation, neither are the ``1/g`` and ``\tfrac{c_p}{g L_v}``
necessary to multiply during the vertical integration as all are positive constants.
While this changes the units of ``P_T, P_q`` they can be reused in the following.

## Deep convection

Following Frierson, 2007 [^Frierson2007] in order to conserve enthalpy we correct
the reference profile for temperature ``T_{ref} \to T_{ref, 2}`` so that ``P_T = P_q``.

```math
T_{ref, 2} = T_{ref} + \frac{1}{\delta p c_p} \int_{p_0}^p_{LZB} c_p (T - T_{ref}) + L_v (q - q_{ref}) dp
```

with the terms inside the integral rearranged compared to Frierson, 2007 to show
that the vertical integral in [First-guess relaxation](@ref) really only has to be computed once.

## Shallow convection

In the following we describe the qref scheme from Frierson, 2007 which corrects
reference profiles for both temperature and humidity to guarantee that ``P_q = 0``,
i.e. no precipitation during convection. In that sense, shallow convection is
non-precipitating.






## References

[^Betts1986]: Betts, A. K., 1986: A new convective adjustment scheme. Part I: Observational and theoretical basis. Quart. J. Roy. Meteor. Soc.,112, 677-691. DOI: [10.1002/qj.49711247307](https://doi.org/10.1002/qj.49711247307)

[^BettsMiller1986]: Betts, A. K. and M. J. Miller, 1986: A new convective adjustment scheme. Part II: Single column tests using GATE wave, BOMEX, ATEX and Arctic air-mass data sets. Quart. J. Roy. Meteor. Soc.,112, 693-709. DOI: [10.1002/qj.49711247308](https://doi.org/10.1002/qj.49711247308)

[^Frierson2007]: Frierson, D. M. W., 2007: The Dynamics of Idealized Convection Schemes and Their Effect on the Zonally Averaged Tropical Circulation. J. Atmos. Sci., 64, 1959-1976. DOI:[10.1175/JAS3935.1](https://doi.org/10.1175/JAS3935.1)