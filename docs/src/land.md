# Land surface model

The land surface in SpeedyWeather is represented through several
components each of which can be changed in a modular way.

- The actual [LandModel](@ref) or [DryLandModel](@ref) describing the equations for soil temperature and soil moisture, and vegetation or rivers.
- [The land-sea mask](@ref)
- The surface [Albedo](@ref)
- The [Orography](@ref)

As these components are largely independent of another,
it is possible to have a white (albedo high) ocean at the
top of Mt Everest, or to paint the Sahara black and moving
it below sea-level. The first one of these bullet points is
discussed below for the others see the respective sections.

## Dry vs wet land

The type hierarchy of theactual land surface model is defined as

```@example land
using SpeedyWeather
using InteractiveUtils # hide
subtypes(SpeedyWeather.AbstractLand)
```
So that a dry land does not have moisture whether the atmosphere
is a `PrimitiveWetModel` (with humidity so it can rain) or a
`PrimitiveDryModel` (without humidity). In contrast a wet land
does have some soil moisture (which can be zero) and in combination
with a `PrimitiveWetModel` will increase with precipitation for
example. You can combine dry and wet land and dry and wet atmosphere
freely.

### LandModel

The default `LandModel` in SpeedyWeather contains
(at the moment other than 2 soil layers are not supported or experimental)

```@example land
spectral_grid = SpectralGrid(trunc=31, nlayers=8)
geometry = LandGeometry(spectral_grid, nlayers=2) # that's also the default, therefore it's optional here
land = LandModel(spectral_grid; geometry)
```

With `LandGeometry` currently used to define the number of soil layers and their layer thickness

```@example land
land.geometry
```

And `land.thermodynamics` used to define thermodynamics-relevant
surface properties as conceptual "constants" used in the
other components of the land surface model `land`.

```@example land
land.thermodynamics
```

To change these you can either mutate the fields or create a new model component `thermodynamics`
passed on to the land model constructor

```@example land
thermodynamics = LandThermodynamics(spectral_grid, heat_conductivity_dry_soil=0.25)
land = LandModel(spectral_grid; thermodynamics)
land.thermodynamics
```

(note that `; thermodynamics` is a Julia shortcut instead of `thermodynamics=thermodynamics`
we often use matching keyword arguments by their name).
Similarly you can change `land.geometry` which however is not mutable.

Generally, a `LandModel` or `DryLandModel` is similarly constructed to a `PrimitiveWetModel`
or `PrimitiveDryModel`. One starts by defining its non-default components,
always passes on the spectral grid as the first argument, and then constructs
the actual land model by passing on those components.
Finally, we use that land model to construct the atmospheric model

```@example land
model = PrimitiveWetModel(spectral_grid; land)
model.land
```

is now the land defined above used when integrating a SpeedyWeather `model`.

In case a non-default number of soil layers is used, the `LandGeometry` also needs to be passed to the `NetCDFOutput` constructor to allocate the correct dimensions of the output variables, when an output is desired:

```julia
output = NetCDFOutput("output.nc", model, land_geometry)
```

### DryLandModel

Alternatively, one can use the `DryLandModel` to explicitly disable any
functionality around soil moisture. By doing so, soil moisture will
always be zero, or treated like zero skipping unnecessary computations.
It is created in the same way (reusing some non-default `thermodynamics` from above)

```@example land
land = DryLandModel(spectral_grid; thermodynamics)
```

but it does not contain soil moisture, vegetation or rivers in contrast to the
`LandModel`.

## Land soil temperature

Currently implemented soil temperatures are

```@example land
subtypes(SpeedyWeather.AbstractLandTemperature)
```

You can use them by passing them on to a
`LandModel`/`DryLandModel` model constructor

```@example land
temperature = LandBucketTemperature(spectral_grid)
land = LandModel(spectral_grid; temperature)
```

such that

```@example land
land.temperature
```

is the `LandBucketTemperature` just defined. Similarly

```@example land
land = DryLandModel(spectral_grid; temperature)
```

if you do not want the land to hold any moisture, vegetation or rivers.

### LandBucketTemperature

`LandBucketTemperature` is a prognostic model of the soil temperature
in the land surface model, interacting two-way with the surface air
temperature. It can warm up through radiation and other surface heat fluxes,
retain thermal energy and release this back to the atmosphere either
in the form of longwave radiative fluxes or sensible heat fluxes
(latent heat fluxes depend on soil moisture, see [Surface fluxes](@ref)).
It is a bucket model such that interaction between neighbouring
grid cells ("buckets") of the land surface only interact through the
atmosphere with another, there are no direct horizontal fluxes
between cells. In the sense of soil moisture, you can fill a bucket
from above with rainfall, it may leak/drain at the bottom but buckets
are laterally isolated from another. A similar concept applies to the
heat fluxes of the `LandBucketTemperature`.

The `LandBucketTemperature` here follows MITgcm's 2-layer model, as defined
[here](https://mitgcm.readthedocs.io/en/latest/phys_pkgs/land.html).
As this is a 2-layer model, `SpectralGrid(nlayers_soil=2)` is required.
The equations are

```math
\begin{aligned}
\Delta z_1 C_1 \frac{\mathrm{d}T_1}{\mathrm{d}t} &= F - \lambda\frac{T_1 - T_2}{(\Delta z_1 + \Delta z_2)/2} \\
\Delta z_2 C_2 \frac{\mathrm{d}T_2}{\mathrm{d}t} &= \lambda\frac{T_1 - T_2}{(\Delta z_1 + \Delta z_2)/2}
\end{aligned}
```

for two layers of thicknesses ``\Delta z_1 = 0.1~m`` (top) and ``\Delta z_2 = 4.0~m`` (layer below)
and respective temperatures ``T_1, T_2``. The total surface downward heat flux
`F` (in ``W/m^2``) forces the surface layer 1, and diffusion scales with the heat conductivity
``\lambda = 0.42 W/m/K`` in the opposite direction of the heat gradient between the layers.
For the parameter choices here it is typical that the surface layer is dominated by the daily cycle,
but the layer below by the seasonal cycle.

The heat capacities ``C_1, C_2`` are diagnosed from the heat capacity of water
``C_w = 4.2 \times 10^6 J/m^3/K`` and dry soil ``C_s = 1.13 \times 10^6 J/m^3/K``
given the soil moistures ``W_1, W_2`` (ratio of available water to field capacity)
of the respective layers.

```math
\begin{aligned}
C_1 &= C_w W_1 \gamma + C_s \\
C_2 &= C_w W_2 \gamma + C_s
\end{aligned}
```

with ``\gamma = 0.24`` being the field capacity per meter soil.

## Land soil moisture

Currently implemented soil moistures are

```@example land
subtypes(SpeedyWeather.AbstractSoilMoisture)
```

You can use them by passing them on to a
`LandModel` (not the `DryLandModel` though which does not have moisture)
model constructor

```@example land
soil_moisture = LandBucketMoisture(spectral_grid)
land = LandModel(spectral_grid; soil_moisture)
```

### LandBucketMoisture

`LandBucketMoisture` defines the prognostic equation for
soil moisture in the land surface model. It is a bucket model
in the sense that every grid cell can fill like a bucket with
rainfall from above dry out by evaporation, can drain into layers below
or into a river runoff. But there is generally no lateral
transport, only through the atmosphere.

The `LandBucketMoisture` here follows MITgcm's 2-layer model, as defined
[here](https://mitgcm.readthedocs.io/en/latest/phys_pkgs/land.html).
As this is a 2-layer model, `SpectralGrid(nlayers_soil=2)` is required.
The equations are

```math
\begin{aligned}
\frac{\mathrm{d}W_1}{\mathrm{d}t} &= \frac{P - E - R}{f_1} + \frac{W_2 - W_1}{\tau} \\
\frac{\mathrm{d}W_2}{\mathrm{d}t} &= -\frac{f_1}{f_2}\frac{W_2 - W_1}{\tau} \\
\end{aligned}
```

for soil moistures ``W_1, W_2`` in the respective layers (1 top, 2 below)
defined as ratio of available water to field capacity, ``f_i = \gamma \Delta z_i``
with ``\gamma = 0.24`` the field capacity per meter soil and
``\Delta z_1 = 0.1~m`` the top layer thickness by default, and
``\Delta z_2 = 4.0~m`` the layer below. The top layer is forced by precipitation
``P`` minus evaporation ``E`` minus river runoff ``R``. The second term is a
diffusion term of soil moisture between the two layers, acting on a time scale
of ``\tau = 2~``days.

At the moment (and generally if not coupled to an ocean model) the river runoff
lets water disappear. ``W_1, W_2`` are bounded by ``[0, 1]`` so that if
more precipitation ``P`` (or in combination with negative evaporation ``E``,
meaning condensation) occurs than the land can hold we compute the excess water
as ``\delta W_1 = W_1 - 1``, set the soil moisture ``W_1 = 1`` back to the
maximum and add a fraction ``p = 0.5`` of that excess water to the layer below
``W_2 = W_2 + p \delta W_1 \tfrac{f_1}{f_2}`` the other half of that excess
water is put into the river runoff.

## Snow model

SpeedyWeather has a simple snow model that allows for snow to lie on land
and impact albedo and isolate air-land fluxes. Snow is created from
the precipitation schemes, can melt on the ground depending on soil temperature
where it is then passed to the soil moisture. We also include a
runoff/relaxation term to prevent snow piling up without bounds as
soil temperatures below freezing would never remove such snow.
In reality this is where one needs a ice sheet model to convert the
snow to ice and simulate ice flow like in a glacier or in the
Greenland and Antarctic ice sheets.

### LandSnowModel

`SnowModel` stores a single snow bucket with depth ``S`` in units of equivalent liquid water height
(`prognostic_variables.land.snow_depth`) and solves the following equation

```math
\frac{dS}{dt} = P - M - R
```
with precipitation ``P``, melting ``M`` and runoff ``R``.
Snow accumulates from the column-integrated large-scale
(currently not from convection) snow precipitation rate ``P``, `snow_rate` (``kg/m²/s``),
and melts once the top soil layer exceeds the melt threshold ``T_{melt}`` (default ``275~K``)
through the term ``M`` and the runoff is implemented as a weak relaxation term
back to 0 with a multi-year time scale.

The available melt energy in the top
layer of thickness ``z₁`` uses the dry-soil heat capacity ``cₛ``:

```math
E_{avail} = cₛ \, \max(T_1 - T_{melt}, 0) \, \frac{z₁}{Δt}.
```

This is the maximum melt rate

```math
\text{melt}_{rate_{max}} = \frac{E_{avail}}{ρ_w Lᵢ},
```

with ``ρ_w`` water density and ``Lᵢ`` latent heat of fusion. But note that
this formulation allows to melt more snow than there actually is, so we need
to cap the amount to not end up with negative snow. We implement this constrain
not for terms individually but for the sum of all terms, see below.

A slow runoff/relaxation term prevents perennial snow packs from growing without bound,

```math
\text{runoff}_{rate_{max}} = \frac{S}{\tau_{runoff}},
```

controlled by `runoff_time_scale` (default ``4`` years). Snowfall, melt and runoff form a raw
tendency ``\text{d}S_{max}`` that is further limited so we never remove more snow than is present
(including what falls during the current step). The actual removal is reported as

```math
\text{snow\_melt\_rate} = \text{melt}_{rate_{max}} + \text{runoff}_{rate_{max}} + \text{excess},
```

in ``kg/m²/s``, where the `excess` term (negative for melting/runoff trying to remove more snow than there is)
only appears when the naive tendency would overdraw the bucket.
`snow_melt_rate` therefore combines true melt with the runoff leak and is zero over ocean points.
Snow depth is clipped to zero and stored as equivalent liquid water height, not physical snow thickness.

The snow budget links into other surface schemes:

- `snow_melt_rate` provides a latent heat sink to the land temperature budget.
- The same flux (converted to ``m/s`` by dividing by ``ρ_w``) feeds the soil moisture tendency alongside rain.
- Snow depth drives the snow-albedo calculation and insulates surface heat/humidity fluxes (see below).

Snow cover over land is diagnosed from snow depth using either a linear ramp
``σₛ = \min(S / \text{snow\_depth\_scale}, 1)`` or the default saturating form

```math
σₛ = \frac{S}{S + \text{snow\_depth\_scale}},
```

set via the `snow_cover` keyword of `LandSnowAlbedo`. Snow then adds an albedo to the land albedo
(diagnosed from vegetation cover)

```math
albedo_{land} = albedo_{land} + σₛ albedo_{snow}
```

so snow depth from the snow bucket immediately brightens land grid cells.
The total albedo is higher over already brighter areas (low vegetation cover)
and lower over darker areas. This somewhat reflects that in forests the
snow cover is broken up and snow lies in between trees.
`DefaultAlbedo` uses `LandSnowAlbedo`; there is currently no time-evolving snow albedo.

## Albedo

Albedo is the surface reflectivity to downward solar shortwave radiation.
A value of 1 indicates that all of the radiative flux is reflected at
the Earth's surface and sent back up through the atmospheric column.
In contrast, a value of 0 means no reflection and all of that radiative
flux is absorbed, typically heating ocean or land surface.
The following albedo's are currently implemented

```@example land
subtypes(SpeedyWeather.AbstractAlbedo)
```
Albedo is generally a 2D global diagnostic variable for ocean and land separately.
You can keep it constant or diagnose it from other fields on every time step:
`OceanSeaIceAlbedo` mixes in sea ice concentration and `LandSnowAlbedo` uses the
snow depth bucket described above. Conceptually albedo is not a static boundary
condition but recomputed each step so transient surface changes (e.g. snowfall)
brighten the surface without losing the underlying bare albedo. If you want to set
albedo manually with `set!` then use `ManualAlbedo` which has its own albedo 2D field
that is copied into the diagnostic variables on every time step.
See example below.
`AlbedoClimatology` does the same but `ManualAlbedo` does not need
to read an albedo from file at initialization.

`Albedo` itself is a container for separate albedo's for ocean and land as
averaging those in grid cells which are partially land, partially ocean will
yield inaccurate results. Think 10% land having a lower heat capacity than
land but being treated with an albedo that comes from 90% ocean.
Not very realistic. The default albedo can be created with

```@example land
albedo = DefaultAlbedo(spectral_grid)
```

and inspected with

```@example land
albedo.ocean
```

and `albedo.land`. You can mix those albedos too, they are internally
two independent albedos that are applied to fluxes separately, e.g.

```@example land
albedo = OceanLandAlbedo(GlobalConstantAlbedo(spectral_grid), AlbedoClimatology(spectral_grid))
```

constructs an albedo that is a global constant (default 0.3) for the ocean
but the `AlbedoClimatology` read from file used for the land.
The first argument for `Albedo` is used for ocean the second for land
but you can use keywords too.

Alternatively you can also drop the separation into ocean/land albedo
(e.g. idealised simulations or aqua planet, rocky planet). And just
use

```@example land
albedo = GlobalConstantAlbedo(spectral_grid)
```

and this definition of the albedo will be used for both ocean and land fluxes.
In all cases you can then pass on the albedo to the model constructor, e.g.

```@example land
albedo = OceanLandAlbedo(GlobalConstantAlbedo(spectral_grid, albedo=0.1), ManualAlbedo(spectral_grid))
set!(albedo.land, (λ, φ) -> 0.2 + 0.3*abs(φ)/90)

model = PrimitiveWetModel(spectral_grid; albedo)
model.albedo
```

The albedo in the `model` is now the one defined just in the lines above,
using a globally constant albedo of 0.1 for the ocean but a higher albedo
over land which also increases to 0.5 towards the poles.

You can always output the land-sea mask weighted albedo with
`add!(model, SpeedyWeather.AlbedoOutput())` or inspect it as follows

```@example land
simulation = initialize!(model)
run!(simulation, steps=1)   # run for a step to "diagnose" albedo = ocean/land weighted

using CairoMakie
(; albedo) = simulation.diagnostic_variables.physics
heatmap(albedo, title="Custom albedo, separately defined for ocean/land")
save("ocean_land_albedo.png", ans) # hide
nothing # hide
```
![Ocean-land albedo](ocean_land_albedo.png)
