# Land-sea mask

The following describes how a custom land-sea mask can be defined.
SpeedyWeather uses a _fractional_ land-sea mask, i.e. for every grid-point

- 1 indicates land
- 0 indicates ocean
- (0,1) indicates a grid-cell partially covered by ocean and land

The land-sea mask determines solely how to weight the surface fluxes
coming from land or from the ocean. For the sensible heat fluxes this uses
land and sea surface temperatures and weights the respective fluxes
proportional to the fractional mask. Similar for evaporation.
You can therefore define an ocean on top of a mountain, or a land without
heat fluxes when the land-surface temperature is not defined, i.e. `NaN`.
Let ``F_L, F_S`` be the fluxes coming from land and sea, respectively.
Then the land-sea mask ``a \in [0,1]`` weights them as follows for the
total flux ``F``

```math
F = aF_L + (1-a)F_S
```

but ``F=F_L`` if the sea flux is NaN (because the ocean temperature is not defined)
and ``F=F_S`` if the land flux is NaN (because the land temperature or soil moisture
is not defined, for sensible heat fluxes or evaporation), and ``F=0`` if both fluxes
are NaN.

Setting the land-sea mask to ocean therefore will disable any fluxes that
may come from land, and vice versa. However, with an ocean-everywhere land-sea mask
you must also define sea surface temperatures everywhere, otherwise the fluxes
in those regions will be zero.

## Manual land-sea mask

You can create the default land-sea mask as follows

```@example landseamask
using SpeedyWeather
spectral_grid = SpectralGrid(trunc=31, nlev=8)
land_sea_mask = LandSeaMask(spectral_grid)
```

which will automatically interpolate the land-sea mask onto grid and resolution
as defined in `spectral_grid` at initialization. The actual mask is in
`land_sea_mask.land_sea_mask` and you can visualise it with

```@example landseamask
model = PrimitiveWetModel(;spectral_grid, land_sea_mask)
simulation = initialize!(model)     # triggers also initialization of model.land_sea_mask
plot(land_sea_mask.mask)
```

Now before you run a simulation you could manually change the land-sea mask by

```@example landseamask
# unpack, this is a flat copy, changing it will also change the mask inside model
(; mask) = land_sea_mask

# ocean everywhere, or
mask .= 0    

# random land-sea mask, or
for i in eachindex(mask)
    mask[i] = rand()     
end

# ocean only between 10˚S and 10˚N
for (j, ring) in enumerate(RingGrids.eachring(mask))
    for ij in ring
        mask[ij] = abs(model.geometry.latd[j]) > 10 ? 1 : 0
    end
end
```

And now you can run the simulation as usual with `run!(simulation)`. Most useful
for the generation of custom land-sea masks in this manual way is probably the
`model.geometry` component which has all sorts of coordinates like `latd`
(latitudes in degrees on rings) or `latds, londs` (latitude, longitude in degrees
for every grid point).

## Earth's land-sea mask

The Earth's [`LandSeaMask`](@ref) has itself the option to load another
land-sea mask from file, but you also have to specify the grid that mask
from files comes on. It will then attempt to read it via `NCDatasets`
and interpolate onto the model grid.

## AquaPlanetMask

Predefined is also the [`AquaPlanetMask`](@ref) which can be created as
```@example landseamask
land_sea_mask = AquaPlanetMask(spectral_grid)
```
and is equivalent to using Earth's [`LandSeaMask`](@ref) but setting
the entire mask to zero afterwards `land_sea_mask.mask .= 0`.

## Custom land-sea mask

You can define a custom land-sea mask type as outlined in [`AbstractLandSeaMask`](@ref).
A custom land-sea mask have to be defined as a new type (`struct` or `mutable struct`)

```julia
CustomMask{NF<:AbstractFloat, Grid<:AbstractGrid{NF}} <: AbstractLandSeaMask{NF, Grid}
```

and needs to have at least a field called `mask::Grid` that uses a `Grid` as defined
by the spectral grid object, so of correct size and with the number format `NF`.
All `AbstractLandSeaMask` have a convenient generator function to be used like
`mask = CustomMask(spectral_grid, option=argument)`, but you may add your own or customize by
defining `CustomMask(args...)` which should return an instance of type `CustomMask{NF, Grid}`
with parameters matching the spectral grid. Then the initialize function has to be extended for
that new mask

```julia
initialize!(mask::CustomMask, model::PrimitiveEquation)
```

which generally is used to tweak the mask.mask grid as you like, using
any other options you have included in `CustomMask` as fields or anything else (preferrably read-only,
because this is only to initialize the land-sea mask, nothing else) from `model`. You can
for example read something from file, set some values manually, or use coordinates from `model.geometry`.

## Time-dependent land-sea mask

It is possible to define an [intrusive callback](@ref intrusive_callbacks) to change the
land-sea mask during integration. The grid in `model.land_sea_mask.mask`
is mutable, meaning you can change the values of grid points in-place but not replace
the entire mask or change its size. If that mask is changed, this will be reflected
in all relevant model components. For example, we can define a callback that
floods the entire planet at the beginning of the 21st century as

```@example landseamask
Base.@kwdef struct MilleniumFlood <: SpeedyWeather.AbstractCallback
    schedule::Schedule = Schedule(DateTime(2000,1,1))
end

# initialize the schedule
function SpeedyWeather.initialize!(
    callback::MilleniumFlood,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::ModelSetup,
)
    initialize!(callback.schedule, progn.clock)
end

function SpeedyWeather.callback!(
    callback::MilleniumFlood,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::ModelSetup,
)
    # escape immediately if not scheduled yet
    isscheduled(callback.schedule, progn.clock) || return nothing

    # otherwise set the entire land-sea mask to ocean
    model.land_sea_mask.mask .= 0
    @info "Everything flooded on $(progn.clock.time)"
end

# nothing needs to be done when finishing
SpeedyWeather.finish!(::MilleniumFlood, args...) = nothing
```

Note that the flooding will take place only at the start of the 21st century,
last indefinitely, but not if the model integration period does not cover that
exact event, see [Schedules](@ref). Initializing a model a few days earlier
would then have this `MilleniumFlood` take place

```@example landseamask
land_sea_mask = LandSeaMask(spectral_grid)      # start with Earth's land-sea mask
model = PrimitiveWetModel(;spectral_grid, land_sea_mask)
add!(model, MilleniumFlood())   # or MilleniumFlood(::DateTime) for any non-default date
model.feedback.verbose = false # hide

simulation = initialize!(model, time=DateTime(1999,12,29))
run!(simulation, period=Day(5))
plot(model.land_sea_mask.mask)
```

And the land-sea mask has succesfully been set to ocean everywhere at the start
of the 21st century.
