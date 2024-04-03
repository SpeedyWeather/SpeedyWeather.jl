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

### Reference profiles

The dry adiabat is

```math
T = T_s (\frac{p}{p_s})^\frac{R}{c_p}
```

The temperature ``T`` of an air parcel at pressure ``p`` is determined by the
temperature ``T_s`` it had at pressure ``p_s`` and the gas constant for dry air
``R`` and the heat capacity ``c_p``. The pseudo adiabat follows the dry adiabat
until saturation is reached, ``q > q^\star``. Then it follows the pseudoadiabatic
lapse rate

```math
\Gamma = -\frac{dT}{dz} = \frac{g}{c_p}\left(
    \frac{1 + \frac{qL_v}{(1-q)^2R_dT_v}}{1 + \frac{qL^2}{(1-q)^2 c_p R_vT^2}}\right)
```




## References

[^Betts1986]: Betts, A. K., 1986: A new convective adjustment scheme. Part I: Observational and theoretical basis. Quart. J. Roy. Meteor. Soc.,112, 677-691. DOI: [10.1002/qj.49711247307](https://doi.org/10.1002/qj.49711247307)

[^BettsMiller1986]: Betts, A. K. and M. J. Miller, 1986: A new convective adjustment scheme. Part II: Single column tests using GATE wave, BOMEX, ATEX and Arctic air-mass data sets. Quart. J. Roy. Meteor. Soc.,112, 693-709. DOI: [10.1002/qj.49711247308](https://doi.org/10.1002/qj.49711247308)

[^Frierson2007]: Frierson, D. M. W., 2007: The Dynamics of Idealized Convection Schemes and Their Effect on the Zonally Averaged Tropical Circulation. J. Atmos. Sci., 64, 1959-1976. DOI:[10.1175/JAS3935.1](https://doi.org/10.1175/JAS3935.1)