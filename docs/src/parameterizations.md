# Parameterizations

The following is an overview of how our parameterizations from a software engineering
perspective are internally defined and how a new parameterization can be accordingly
implemented. For the mathematical formulation and the physics they represent see

- [Vertical diffusion](@ref)
- [Convection](@ref)
- [Large-scale condensation](@ref)
- [Radiation](@ref)
- [Surface fluxes](@ref) 

We generally recommend reading [Extending SpeedyWeather](@ref) first,
which explains the logic of how to extend many of the components in
SpeedyWeather. The same logic applies here and we will not iterate on
many of the details, but want to highlight which abstract supertype
new parameterizations have to subtype respectively and which functions
and signatures they have to extend.

In general, every parameterization "class" (e.g. convection) is just a
*conceptual* class for clarity. You can define a custom convection
parameterization that acts as a longwave radiation and vice versa.
This also means that if you want to implement a parameterization
that does not fit into any of the "classes" described here you
can still implement it under any name and any class. From a software
engineering perspective they are all the same except that they
are executed in the order as outlined in [Pass on to model construction](@ref).
That's also why below we write for every parameterization
"expected to write into `some.array_name`" as this would correspond
conceptually to this class, but no hard requirement exists that a
parameterization actually does that.

We start by highlighting some general do's and don'ts for
parameterization before listing specifics for individual parameterizations.

## Use `ColumnVariables` work arrays

When defining a new (mutable) parameterization with (mutable) fields
do make sure that is constant during the model integration. While you
can and are encouraged to use the `initialize!` function to precompute
arrays (e.g. something that depends on latitude using `model.geometry.latd`)
these should not be used as work arrays on every time step of the
model integration. The reason is that the parameterization are executed
in a parallel loop over all grid points and a mutating parameterization object
would create a race condition with undefined behaviour.

Instead, `column::ColumnVariables` has several work arrays that you can
reuse `column.a` and `.b`, `.c`, `.d`. Depending on the number of threads
there will be several `column` objects to avoid the race condition if
several threads would compute the parameterizations for several columns
in parallel. An example is [the Simplified Betts-Miller](@ref BettsMiller)
convection scheme which needs to compute reference profiles which should
not live inside the `model.convection` object (as there's always only one of those).
Instead this parameterization does the following inside
`convection!(column::ColumnVariables, ...)`

```julia
# use work arrays for temp_ref_profile, humid_ref_profile
temp_ref_profile = column.a
humid_ref_profile = column.b
```

These work arrays have an unknown state so you should overwrite every entry
and you also should not use them to retain information after that parameterization
has been executed.

## Accumulate do not overwrite

Every parameterization either computes tendencies directly or indirectly via
fluxes (upward or downward, see [Fluxes to tendencies](@ref)). Both of these are
arrays in which *every* parameterization writes into, meaning they should be
*accumulated not overwritten*. Otherwise any parameterization that executed
beforehand is effectively disabled. Hence, do

```julia
column.temp_tend[k] += something_you_calculated
```

not `column.temp_tend[k] = something_you_calculated` which would overwrite
any previous tendency. The tendencies are reset to zero for every grid point
at the beginning of the evaluation of the parameterization for that grid point,
meaning you can do `tend += a` even for the first parameterization that writes into
a given tendency as this translates to `tend = 0 + a`.

## Pass on to model construction

After defining a (custom) parameterization it is recommended to also define
a generator function that takes in the `SpectralGrid` object
(see [How to run SpeedyWeather.jl](@ref)) as first (positional) argument,
all other arguments can then be passed on as keyword arguments with defaults
defined. Creating the default convection parameterization for example would be
```@example parameterization
using SpeedyWeather
spectral_grid = SpectralGrid(trunc=31, nlev=8)
convection = SimplifiedBettsMiller(spectral_grid, time_scale=Hour(4))
```
Further keyword arguments can be added or omitted all together (using the default
setup), only the `spectral_grid` is required. We have chosen here the name
of that parameterization to be `convection` but you could choose any other name too.
However, with this choice one can conveniently pass on via matching
the `convection` field inside `PrimitiveWetModel`, see
[Keyword Arguments](https://docs.julialang.org/en/v1/manual/functions/#Keyword-Arguments).

```@example parameterization
model = PrimitiveWetModel(;spectal_grid, convection)
```
otherwise we would need to write

```@example parameterization
my_convection = SimplifiedBettsMiller(spectral_grid)
model = PrimitiveWetModel(;spectal_grid, convection=my_convection)
```
The following is an overview of what the parameterization fields inside the
model are called. See also [Tree structure](@ref), and therein [PrimitiveDryModel](@ref)
and [PrimitiveWetModel](@ref)

- `model.boundary_layer_drag`
- `model.temperature_relaxation`
- `model.vertical_diffusion`
- `model.convection`
- `model.large_scale_condensation` (`PrimitiveWetModel` only)
- `model.shortwave_radiation`
- `model.longwave_radiation`
- `model.surface_thermodynamics`
- `model.surface_wind`
- `model.surface_heat_flux`
- `model.surface_evaporation` (`PrimitiveWetModel` only)

Note that the parameterizations are executed in the order of the list
above. That way, for example, radiation can depend on calculations in
large-scale condensation but not vice versa (only at the next time step).

## Custom boundary layer drag

A boundary layer drag can serve two purposes: (1) Actually define a tendency
to the momentum equations that acts as a drag term, or (2) calculate the
drag coefficient ``C`` in `column.boundary_layer_drag` that is used
in the [Surface fluxes](@ref).

- subtype `CustomDrag <: AbstractBoundaryLayer`
- define `initialize!(::CustomDrag, ::PrimitiveEquation)`
- define `boundary_layer_drag!(::ColumnVariables, ::CustomDrag, ::PrimitiveEquation)`
- expected to accumulate (`+=`) into `column.u_tend` and `column.v_tend`
- or calculate `column.boundary_layer_drag` to be used in surface fluxes

## Custom temperature relaxation

By default, there is no temperature relaxation in the primitive equation
models (i.e. `temperature_relaxation = NoTemperatureRelaxation()`).
This parameterization exists for the [Held-Suarez forcing](@ref).

- subtype `CustomTemperatureRelaxation <: AbstractTemperatureRelaxation`
- define `initialize!(::CustomTemperatureRelaxation, ::PrimitiveEquation)`
- define `temperature_relaxation!(::ColumnVariables, ::CustomTemperatureRelaxation, ::PrimitiveEquation)`
- expected to accumulate (`+=`) into `column.temp_tend`

## Custom vertical diffusion

While vertical diffusion may be applied to temperature (usually via some form of static energy
to account for adiabatic diffusion), humidity and/or momentum, they are grouped together.
You can define a vertical diffusion for only one or several of these variables where you then
can internally call functions like `diffuse_temperature!(...)` for each variable.
For vertical diffusion

- subtype `CustomVerticalDiffusion <: AbstractVerticalDiffusion`
- define `initialize!(::CustomVerticalDiffusion, ::PrimitiveEquation)`
- define `vertical_diffusion!(::ColumnVariables, ::CustomVerticalDiffusion, ::PrimitiveEquation)`
- expected to accumulate (`+=`) into `column.temp_tend`, and similarly for `humid`, `u`, and/or `v`
- or using fluxes like `column.flux_temp_upward`

## Custom convection

- subtype `CustomConvection <: AbstractConvection`
- define `initialize!(::CustomConvection, ::PrimitiveEquation)`
- define `convection!(::ColumnVariables, ::CustomConvection, ::PrimitiveEquation)`
- expected to accumulate (`+=`) into `column.temp_tend` and `column.humid_tend`
- or using fluxes, `.flux_temp_upward` or similarly for `humid` or `downward`

Note that we define convection here for a model of type `PrimitiveEquation`, i.e.
both dry and moist convection. If your `CustomConvection` only makes sense for
one of them use `::PrimitiveDry` or `::PrimitiveWet` instead.

## Custom large-scale condensation

- subtype `CustomCondensation <: AbstractCondensation`
- define `initialize!(::CustomCondensation, ::PrimitiveWet)`
- define `condensation!(::ColumnVariables, ::CustomCondensation, ::PrimitiveWet)`
- expected to accumulate (`+=`) into `column.humid_tend` and `column.temp_tend`
- or using fluxes, `.flux_humid_downward` or similarly for `temp` or `upward`

## Custom radiation

`AbstractRadiation` has two subtypes, `AbstractShortwave` and `AbstractLongwave`
representing two (from a software engineering perspective) independent
parameterizations that are called one after another (short then long).
For shortwave

- subtype `CustomShortwave <: AbstractShortwave`
- define `initialize!(::CustomShortwave, ::PrimitiveEquation)`
- define `shortwave_radiation!(::ColumnVariables, ::CustomShortwave, ::PrimitiveEquation)`
- expected to accumulate (`+=`) into `column.flux_temp_upward`, `.flux_temp_downward`
- or directly into the tendency `.temp_tend`

For longwave this is similar but using `<: AbstractLongwave` and `longwave_radiation!`.

## Custom surface fluxes

[Surface fluxes](@ref) are the most complicated to customize as they depend on
the [Ocean](@ref) and Land model, [The land-sea mask](@ref), and by default the
[## Bulk Richardson-based drag coefficient](@ref), see [Custom boundary layer drag](@ref).
The computation of the surface fluxes is split into three components that
are called one after another

1. Surface thermodynamics to calculate the surface values of lowermost layer variables
- subtype `CustomSurfaceThermodynamics <: AbstractSurfaceThermodynamics`
- define `initialize!(::CustomSurfaceThermodynamics, ::PrimitiveEquation)`
- define `surface_thermodynamics!(::ColumnVariables, ::CustomSurfaceThermodynamics, ::PrimitiveEquation)`
- expected to set `column.surface_temp`, `.surface_humid`, `.surface_air_density`
2. Surface wind to calculate wind stress (momentum flux) as well as surface wind used for
3. Surface (sensible) heat flux
4. Surface evaporation


