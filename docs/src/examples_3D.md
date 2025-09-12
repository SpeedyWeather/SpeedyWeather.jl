# Examples 3D

The following showcases several examples of SpeedyWeather.jl simulating
the [Primitive equations](@ref primitive_equation_model) with and without
humidity and with and without physical parameterizations.

See also [Examples 2D](@ref Examples) for examples with the
[Barotropic vorticity equation](@ref barotropic_vorticity_model) and the
[shallow water model](@ref shallow_water_model).

## Jablonowski-Williamson baroclinic wave

```@example jablonowski
using SpeedyWeather
spectral_grid = SpectralGrid(trunc=31, nlayers=8, Grid=FullGaussianGrid, dealiasing=3)

orography = ZonalRidge(spectral_grid)
initial_conditions = InitialConditions(
    vordiv = ZonalWind(),
    temp = JablonowskiTemperature(),
    pres = ConstantPressure())

model = PrimitiveDryModel(spectral_grid; orography, initial_conditions, physics=false)
simulation = initialize!(model)
run!(simulation, period=Day(9))
nothing # hide
```

The Jablonowski-Williamson baroclinic wave test case[^JW06] using the
[Primitive equation model](@ref primitive_equation_model) particularly the dry model,
as we switch off all physics with `physics=false`.
We want to use 8 vertical levels, and a lower resolution of T31 on a
[full Gaussian grid](@ref FullGaussianGrid).
The Jablonowski-Williamson initial conditions are `ZonalWind` for vorticity and divergence
(curl and divergence of ``u, v``), `JablonowskiTemperature` for temperature and
`ZeroInitially` for pressure. The orography is just a `ZonalRidge`.
There is no forcing and the initial conditions are baroclinically unstable which kicks
off a wave propagating eastward. This wave becomes obvious when visualised with

```@example jablonowski
using CairoMakie

vor = simulation.diagnostic_variables.grid.vor_grid[:, end]
heatmap(vor, title="Surface relative vorticity")
save("jablonowski.png", ans) # hide
nothing # hide
```
![Jablonowski plot](jablonowski.png)

## Held-Suarez forcing

```@example heldsuarez
using SpeedyWeather
spectral_grid = SpectralGrid(trunc=31, nlayers=8)

# construct model with only Held-Suarez forcing, no other physics
model = PrimitiveDryModel(
    spectral_grid,

    # Held-Suarez forcing and drag
    temperature_relaxation = HeldSuarez(spectral_grid),
    boundary_layer_drag = LinearDrag(spectral_grid),

    # switch off other physics
    convection = nothing,
    shortwave_radiation = nothing,
    longwave_radiation = nothing,
    vertical_diffusion = nothing,

    # switch off surface fluxes (makes ocean/land/land-sea mask redundant)
    surface_wind = nothing,
    surface_heat_flux = nothing,

    # use Earth's orography
    orography = EarthOrography(spectral_grid)
)

simulation = initialize!(model)
run!(simulation, period=Day(20))
nothing # hide
```

The code above defines the Held-Suarez forcing [^HS94] in terms of temperature relaxation
and a linear drag term that is applied near the planetary boundary but switches off
all other physics in the primitive equation model without humidity.
Switching off the surface wind would also automatically turn off the surface evaporation
(not relevant in the primitive _dry_ model) and sensible heat flux as that one is proportional
to the surface wind (which is zero with `nothing`).
But to also avoid the calculation being run at all we use `nothing` for model components
passed to the model constructor.
Some model components use a `NoSomething` but most can just be set to `nothing`.
`NoSomething`s do not require the spectral grid to be passed on, but as a convention
we allow every model component to have it for construction even if not required.

Visualising surface temperature with

```@example heldsuarez
using CairoMakie

temp = simulation.diagnostic_variables.grid.temp_grid[:, end]
heatmap(temp, title="Surface temperature [K]", colormap=:thermal)

save("heldsuarez.png", ans) # hide
nothing # hide
```
![Held-Suarez](heldsuarez.png)

## Aquaplanet

```@example aquaplanet
using SpeedyWeather

# components
spectral_grid = SpectralGrid(trunc=31, nlayers=8)
ocean = AquaPlanet(spectral_grid, temp_equator=302, temp_poles=273)
land_sea_mask = AquaPlanetMask(spectral_grid)
orography = NoOrography(spectral_grid)

# create model, initialize, run
model = PrimitiveWetModel(spectral_grid; ocean, land_sea_mask, orography)
simulation = initialize!(model)
run!(simulation, period=Day(20))
nothing # hide
```

Here we have defined an aquaplanet simulation by
- creating an `ocean::AquaPlanet`. This will use constant sea surface temperatures that only vary with latitude.
- creating a `land_sea_mask::AquaPlanetMask` this will use a land-sea mask with `false`=ocean everywhere.
- creating an `orography::NoOrography` which will have no orography and zero surface geopotential.

All passed on to the model constructor for a `PrimitiveWetModel`, we have now a model with humidity
and physics parameterization as they are defined by default (typing `model` will give you an overview
of its components). We could have change the `model.land` and `model.vegetation` components too,
but given the land-sea masks masks those contributions to the surface fluxes anyway, this is not
necessary. Note that neither sea surface temperature, land-sea mask
or orography have to agree. It is possible to have an ocean on top of a mountain.
For an ocean grid-cell that is (partially) masked by the land-sea mask, its value will
be (fractionally) ignored in the calculation of surface fluxes (potentially leading
to a zero flux depending on land surface temperatures).

Now with the following we visualize the surface humidity after the 50 days of
simulation. We use 50 days as without mountains it takes longer for the initial conditions to
become unstable. The surface humidity shows small-scale patches in the tropics, which is a result
of the convection scheme, causing updrafts and downdrafts in both humidity and temperature.

```@example aquaplanet
using CairoMakie

humid = simulation.diagnostic_variables.grid.humid_grid[:, end]
heatmap(humid, title="Surface specific humidity [kg/kg]", colormap=:oslo)

save("aquaplanet.png", ans) # hide
nothing # hide
```
![Aquaplanet](aquaplanet.png)

## Aquaplanet without (deep) convection

Now we want to compare the previous simulation to a simulation without
deep convection, called `DryBettsMiller`, because it is the
[Betts-Miller convection](@ref BettsMiller) but with humidity set to zero
in which case the convection is always non-precipitating _shallow_
(because the missing latent heat release from condensation makes it shallower)
convection. In fact, this convection is the default when using
the `PrimitiveDryModel`. Instead of redefining the [Aquaplanet](@ref) setup again,
we simply reuse these components `spectral_grid`, `ocean`, `land_sea_mask`
and `orography` (because `spectral_grid` hasn't changed this is possible).

```@example aquaplanet
# Execute the code from Aquaplanet above first!
convection = DryBettsMiller(spectral_grid, time_scale=Hour(4))

# reuse other model components from before
model = PrimitiveWetModel(spectral_grid; ocean, land_sea_mask, orography, convection)

simulation = initialize!(model)
run!(simulation, period=Day(20))

humid = simulation.diagnostic_variables.grid.humid_grid[:, end]
heatmap(humid, title="No deep convection: Surface specific humidity [kg/kg]", colormap=:oslo)
save("aquaplanet_nodeepconvection.png", ans) # hide
nothing # hide
```

But we also want to compare this to a setup where convection is completely
disabled, i.e. `convection = nothing`

```@example aquaplanet
# Execute the code from Aquaplanet above first!
convection = nothing

# reuse other model components from before
model = PrimitiveWetModel(spectral_grid; ocean, land_sea_mask, orography, convection)

simulation = initialize!(model)
run!(simulation, period=Day(20))

humid = simulation.diagnostic_variables.grid.humid_grid[:, end]
heatmap(humid, title="No convection: Surface specific humidity [kg/kg]", colormap=:oslo)
save("aquaplanet_noconvection.png", ans) # hide
nothing # hide
```

And the comparison looks like

![Aquaplanet, no deep convection](aquaplanet_nodeepconvection.png)
![Aquaplanet, no convection](aquaplanet_noconvection.png)

## Large-scale vs convective precipitation

```@example precipitation
using SpeedyWeather

# components
spectral_grid = SpectralGrid(trunc=31, nlayers=8)
large_scale_condensation = ImplicitCondensation(spectral_grid, snow=false)
convection = SimplifiedBettsMiller(spectral_grid)

# create model, initialize, run
model = PrimitiveWetModel(spectral_grid; large_scale_condensation, convection)
simulation = initialize!(model)
run!(simulation, period=Day(10))
nothing # hide
```

We run the default `PrimitiveWetModel` with `ImplicitCondensation` as large-scale condensation
(see [Implicit large-scale condensation](@ref)) and the `SimplifiedBettsMiller`
for convection (see [Simplified Betts-Miller](@ref BettsMiller)). These schemes
have some additional parameters, we leave them as default for now, but you could
do `ImplicitCondensation(spectral_grid, relative_humidity_threshold = 0.8)` to
let it rain at 80% instead of 100% relative humidity. We now want to analyse
the precipitation that comes from these parameterizations.

Note that the following only considers liquid precipitation for simplicity.
We set `snow=false` in `ImplicitCondensation` and therefore deal with rain only.

```@example precipitation
using CairoMakie

(; rain_large_scale, rain_convection) = simulation.diagnostic_variables.physics
m2mm = 1000     # convert from [m] to [mm]
heatmap(m2mm*rain_large_scale, title="Large-scale precipiation (rain) [mm]: Accumulated over 10 days", colormap=:dense)
save("large-scale_precipitation_acc.png", ans) # hide
nothing # hide
```
![Large-scale precipitation](large-scale_precipitation_acc.png)

Precipitation (rain, both large-scale and convective) are written into the
`simulation.diagnostic_variables.physics` which, however, accumulate all precipitation
during simulation. In the NetCDF output, precipitation rate (in mm/hr) is calculated
from accumulated precipitation as a post-processing step.
More interactively, you can also reset these accumulators and integrate for another 6 hours
to get the precipitation only in that period.

```@example precipitation
# reset accumulators and simulate 6 hours
simulation.diagnostic_variables.physics.rain_large_scale .= 0
simulation.diagnostic_variables.physics.rain_convection .= 0
run!(simulation, period=Hour(6))

# visualise, rain_* arrays are flat copies, no need to read them out again!
m2mm_hr = (1000*Hour(1)/Hour(6))    # convert from [m] to [mm/hr]
heatmap(m2mm_hr*rain_large_scale, title="Large-scale precipiation (rain) [mm/hr]", colormap=:dense)
save("large-scale_precipitation.png", ans) # hide
heatmap(m2mm_hr*rain_convection, title="Convective precipiation (rain) [mm/hr]", colormap=:dense)
save("convective_precipitation.png", ans) # hide
nothing # hide
```
![Large-scale precipitation](large-scale_precipitation.png)
![Convective precipitation](convective_precipitation.png)

As the precipitation fields are accumulated meters over the integration period
we divide by 6 hours to get a precipitation rate ``[m/s]``
but then multiply with 1 hour and 1000 to get the typical precipitation unit of ``[mm/hr]``.

## References

[^JW06]: Jablonowski, C. and Williamson, D.L. (2006), A baroclinic instability test case for atmospheric model dynamical cores. Q.J.R. Meteorol. Soc., 132: 2943-2975. DOI:[10.1256/qj.06.12](https://doi.org/10.1256/qj.06.12)

[^HS94]: Held, I. M. & Suarez, M. J. A Proposal for the Intercomparison of the Dynamical Cores of Atmospheric General Circulation Models. Bulletin of the American Meteorological Society 75, 1825-1830 (1994). DOI:[10.1175/1520-0477(1994)075<1825:APFTIO>2.0.CO;2](https://doi.org/10.1175/1520-0477(1994)075<1825:APFTIO>2.0.CO;2)
