# Saturation vapor pressure

The saturation vapor pressure ``e_s(T)`` is the partial pressure of water vapor in equilibrium
with a flat surface of liquid water (or ice) at temperature ``T``. It defines the saturation
specific humidity

```math
q_s(T, p) = \epsilon \frac{e_s(T)}{p}, \quad \epsilon = \frac{R_d}{R_v}
```

that is used by [Large-scale condensation](@ref), [Convection](@ref), the cloud diagnostics in
[Radiation](@ref) and the [Surface fluxes](@ref) of humidity. Several formulations of ``e_s(T)``
are implemented, which differ in accuracy and in whether they include an ice phase. The
formulation is a property of the atmosphere, stored in `model.atmosphere.saturation`.

## Clausius-Clapeyron and its approximations

The Clausius-Clapeyron (or more precisely Clapeyron) relation is exact for any phase transition

```math
\frac{de_s}{dT} = \frac{L}{T (v_v - v_l)}
```

with latent heat ``L`` and the specific volumes ``v_v, v_l`` of vapor and liquid. Three
approximations lead to the formula commonly called Clausius-Clapeyron:

1. the condensed volume is negligible, ``v_l \ll v_v`` (error ``\sim 10^{-5}``),
2. the vapor is an ideal gas, ``v_v = R_v T / e_s`` (error ``\sim 0.1\%``),
3. the latent heat ``L`` is constant.

The last is by far the largest error: ``L_v`` decreases by about 2.4 kJ/kg per Kelvin.
Integrating with a constant ``L_v`` overestimates ``e_s`` by about 2.4% at 30˚C and 4% at 40˚C.
Other formulations either integrate with ``L`` linear in temperature, or are empirical fits
(Magnus or Tetens type) to reference data of pure water vapor, e.g.
[MurphyKoop2005](@citet). The presence of dry air changes ``e_s`` by only ≈0.4%,
the so-called enhancement factor, which is ignored here.

## Implemented formulations

```@example svp
using InteractiveUtils # hide
using SpeedyWeather
subtypes(SpeedyWeather.AbstractSaturationVaporPressure)
```

| Formulation | Formula | Phases | Notes |
|---|---|---|---|
| `ClausiusClapeyronConstantL` (default) | ``e_0 \exp\left(\frac{L_v}{R_v}\left(\frac{1}{T_0} - \frac{1}{T}\right)\right)`` | liquid | ``L_v, R_v`` from the atmosphere; +2–4% at 30–40˚C |
| `ClausiusClapeyronLinearL` | ``e_0 \left(\frac{T}{T_0}\right)^a \exp\left(b\left(\frac{1}{T_0} - \frac{1}{T}\right)\right)`` | liquid | [Romps2017](@citet), ``L`` linear in ``T``, ≈0.1% accurate |
| `IFSMixedPhase` | ``\alpha e_{s,l} + (1 - \alpha) e_{s,i}``, ``e_{s,x} = e_0 \exp\left(a_x \frac{T - T_0}{T - b_x}\right)`` | liquid, ice, mixed | as in the ECMWF IFS [IFSCycle49r1](@cite), ``\alpha = \left(\frac{T - T_{ice}}{T_0 - T_{ice}}\right)^2`` between -23˚C and 0˚C |
| `MagnusBolton` | ``e_0 \exp\left(a \frac{T_c}{T_c + b}\right)``, ``T_c`` in ˚C | liquid | [Bolton1980](@citet), ≈0.1% accurate between -35˚C and 35˚C |

For Romps (2017) the exponents are ``a = (c_{pv} - c_{vl})/R_v`` and
``b = (E_{0v} - (c_{vv} - c_{vl}) T_0)/R_v`` with ``T_0, e_0`` at the triple point.
All coefficients are fields of the respective formulation and can be changed by keyword.
In Pa for the four formulations in the order above
```@example svp
atmosphere = EarthAtmosphere(Float64)
formulations = (ClausiusClapeyronConstantL, ClausiusClapeyronLinearL, IFSMixedPhase, MagnusBolton)
for T in (233.15, 253.15, 273.15, 293.15, 313.15)
    e = [SpeedyWeather.saturation_vapor_pressure(F(Float64), T, atmosphere) for F in formulations]
    println(round(T - 273.15), "˚C: ", join(round.(e, digits = 2), ", "))
end
```

For reference, [MurphyKoop2005](@citet) give 19.0 (12.8 over ice), 125.6 (103.2 over ice),
611.2, 2339 and 7385 Pa at these temperatures.

Here is an overview of the accuracy of each formulation against [MurphyKoop2005](@citet) quantified as relative error.
The IFS formulation considers a mixed phase between -23˚ and 0˚C which we therefore also apply to the the 
[MurphyKoop2005](@citet) reference and other formulations too to compare on a single axis from -100˚ to +100˚
covering both vapor saturation over ice and liquid water.

<img width="2000" height="1040" alt="image" src="https://github.com/user-attachments/assets/dcfe9a46-0708-484f-bc52-a9e9724624da" />

## User interface

The default is Clausius-Clapeyron with a constant latent heat (as in previous versions).
Choose another formulation through the atmosphere

```@example svp
spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
saturation = IFSMixedPhase(spectral_grid)       # or MagnusBolton, ClausiusClapeyronLinearL, ...
atmosphere = EarthAtmosphere(spectral_grid; saturation)
model = PrimitiveWetModel(spectral_grid; atmosphere)
model.atmosphere.saturation
```

Coefficients are changed by keyword, e.g. `ClausiusClapeyronConstantL(spectral_grid, e₀ = 611.2)`.
The functions

```@example svp
T = 280f0       # temperature [K]
p = 85000f0     # pressure [Pa]
SpeedyWeather.saturation_vapor_pressure(T, atmosphere),         # [Pa]
SpeedyWeather.grad_saturation_vapor_pressure(T, atmosphere),    # deₛ/dT [Pa/K]
SpeedyWeather.saturation_humidity(T, p, atmosphere)             # [kg/kg]
```

evaluate the chosen formulation, which each also provides its derivative via
`SpeedyWeather.grad(saturation, T, atmosphere)`. The derivative is analytical, cheap and
branch-free, so it can be used inside kernels and differentiated through.

The derivative of the saturation humidity follows from it as

```math
\frac{\partial q_s}{\partial T} = \frac{\epsilon}{p} \frac{de_s}{dT}
```

available as `SpeedyWeather.grad_saturation_humidity(T, p, atmosphere)`. The
[Large-scale condensation](@ref) uses it for its implicit correction and for limiting
reevaporation, so that both stay consistent with the chosen formulation.

!!! note "Convection still assumes constant latent heat"
    The moist-adiabatic lapse rate in [Convection](@ref) still uses
    ``\partial q_s/\partial T = q_s L_v/(R_v T^2)`` implicitly, which is exact only for the
    default `ClausiusClapeyronConstantL`.

## Define your own formulation

A new formulation is a subtype of `SpeedyWeather.AbstractSaturationVaporPressure` with
a method for `saturation_vapor_pressure(svp, T, atmosphere)` and `grad(svp, T, atmosphere)`

```@example svp
@kwdef struct ConstantSaturation{NF} <: SpeedyWeather.AbstractSaturationVaporPressure
    e₀::NF = 1000
end

SpeedyWeather.saturation_vapor_pressure(svp::ConstantSaturation, T, atmosphere) = svp.e₀
SpeedyWeather.grad(svp::ConstantSaturation, T, atmosphere) = zero(svp.e₀)

atmosphere = EarthAtmosphere(spectral_grid, saturation = ConstantSaturation{spectral_grid.NF}())
SpeedyWeather.saturation_vapor_pressure(T, atmosphere)
```

The `atmosphere` argument provides constants like `atmosphere.latent_heat_condensation` and
`atmosphere.R_vapor` if a formulation should be consistent with them. A struct with a single
field (as here) has to be created with its number format, e.g. `ConstantSaturation{Float32}()`,
as `ConstantSaturation(spectral_grid)` is ambiguous with the struct's positional constructor. Use `ifelse` rather than
`if` to stay branch-free.
