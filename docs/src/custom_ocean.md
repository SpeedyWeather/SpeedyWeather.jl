# Ocean

The ocean in SpeedyWeather.jl is defined with two horizontal fields in the
prognostic variables which has a field `ocean`, i.e. `simulation.variables.prognostic.ocean`.

- `ocean.sea_surface_temperature` with units of Kelvin [K].
- `ocean.sea_ice_concentration` with units of area fraction [1].

Both are two-dimensional grids using the same grid type and resolution as
the dynamical core. So both sea surface temperature and sea ice concentration
are globally defined but their relative contribution to surface fluxes is
weighted by [The land-sea mask](@ref). Sea surface temperature must be
finite everywhere, including at (partially or fully) land-masked grid cells:
the land-sea mask only weights the resulting flux multiplicatively, it does
not guard against propagating a `NaN` through it, so `NaN` sea surface
temperatures are not supported. Ocean models therefore fill grid cells
without "real" ocean with a fallback temperature (e.g. `SlabOcean`'s
`land_temperature` option) rather than leaving them `NaN`. The sea ice
concentration is simply set to zero everywhere where there is no sea ice.

Note that neither sea surface temperature, land-sea mask
or orography have to agree. It is possible to have an ocean on top of a mountain.
For an ocean grid-cell that is (partially) masked by the land-sea mask, its value will
still be (fractionally) weighted into the calculation of surface fluxes
(potentially leading to a small but nonzero flux depending on land surface
temperatures), so a physically sensible, finite value should be provided
everywhere, even where the land-sea mask is 1 (fully land).

# Custom ocean model

Now the ocean model is expected to change `ocean.sea_surface_temperature`
and/or `ocean.sea_ice_concentration` on a given time step.
A new ocean model has to be defined as

```julia
struct CustomOceanModel <: AbstractOcean
    # fields, coefficients, whatever is constant, 
end
```

and can have parameters like `CustomOceanModel{T}` and any fields. 
`CustomOceanModel` then needs to extend the following functions

```julia
function initialize!(
    ocean_model::CustomOceanModel,
    model::PrimitiveEquation)
    # your code here to initialize the ocean model itself
    # you can use other fields from model, e.g. model.geometry
end

function initialize!(
    vars::Variables,
    ocean_model::CustomOceanModel,
    model::PrimitiveEquation,
)
    # your code here to initialize the prognostic variables for the ocean
    # namely, ocean.sea_surface_temperature, ocean.sea_ice_concentration, e.g.
    # vars.prognostic.ocean.sea_surface_temperature .= 300      # 300K everywhere
end
```

Note that the first is only to initialize the `CustomOceanModel` not the
prognostic variables. For example `SeasonalOceanClimatology <: AbstractOcean`
loads in climatological sea surface temperatures for every time month in the
first `initialize!` but only writes them (given `time`) into the prognostic
variables in the second `initialize!`. They are internally therefore also
called in that order. Note that the function signatures should not be changed
except to define a new method for `CustomOceanModel` or whichever name you chose.

Then you have to extend the `timestep!` function which has a signature like
```julia
function timestep!(
    vars::Variables,
    ocean_model::CustomOceanModel,
    model::PrimitiveEquation,
)
    # your code here to change the vars.prognostic.ocean.sea_surface_temperature
end
```
which is called on every time step before the land and before the parameterization
and therefore also before the dynamics. You can schedule the execution with
[Schedules](@ref) or you can use the `ocean.time` time to determine when last
the ocean time step was executed and whether it should be executed now, e.g.
`(time - ocean.time) < ocean_model.Δt && return nothing` would not execute
unless the period of the `ocean_model.Δt` time step has passed. Note that
the `ocean.sea_surface_temperature` or `.sea_ice_concentration` are unchanged
if the ocean time step is not executed, meaning that the sea surface temperatures
for example can lag behind the dynamical core for some days essentially assuming
constant temperatures throughout that period. Any ocean
model with constant temperatures and sea ice should just `return nothing`.