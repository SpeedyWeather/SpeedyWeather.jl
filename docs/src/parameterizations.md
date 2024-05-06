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

Note that the parameterizations are executed in the order of the list
above. That way, for example, radiation can depend on calculations in
large-scale condensation but not vice versa (only at the next time step).
We start by highlighting some general do's and dont's for
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

## Custom surface fluxes

[Surface fluxes](@ref) are the most complicated to customize as they depend on
the [Ocean](@ref) and Land model and [The land-sea mask](@ref). We give a brief overview,
