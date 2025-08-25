# Ocean

The ocean in SpeedyWeather.jl is defined with two horizontal fields in the
prognostic variables which has a field `ocean`, i.e. `simulation.prognostic_variables.ocean`.

- `ocean.sea_surface_temperature` with units of Kelvin [K].
- `ocean.sea_ice_concentration` with units of area fraction [1].

Both are two-dimensional grids using the same grid type and resolution as
the dynamical core. So both sea surface temperature and sea ice concentration
are globally defined but their mask is defined with [The land-sea mask](@ref).
However, one should still set grid cells where the sea surface temperature
is not defined to `NaN` in which case any fluxes are zero. This is important
when a fractional land-sea mask does not align with the sea surface
temperatures to not produce unphysical fluxes. The sea ice concentration
is simply set to zero everywhere where there is no sea ice.

Note that neither sea surface temperature, land-sea mask
or orography have to agree. It is possible to have an ocean on top of a mountain.
For an ocean grid-cell that is (partially) masked by the land-sea mask, its value will
be (fractionally) ignored in the calculation of surface fluxes (potentially leading
to a zero flux depending on land surface temperatures). For an ocean grid cell
that is `NaN` but not masked by the land-sea mask, its value is always ignored.

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
    ocean::PrognosticVariablesOcean,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    ocean_model::CustomOceanModel,
    model::PrimitiveEquation,
)
    # your code here to initialize the prognostic variables for the ocean
    # namely, ocean.sea_surface_temperature, ocean.sea_ice_concentration, e.g.
    # ocean.sea_surface_temperature .= 300      # 300K everywhere
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
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    ocean_model::CustomOceanModel,
    model::PrimitiveEquation,
)
    # your code here to change the progn.ocean.sea_surface_temperature
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