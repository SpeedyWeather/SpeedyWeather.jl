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
are executed in the order as defined in [`get_parametrizations](@ref).
That's also why below we write for every parameterization
"expected to write into `some.array_name`" as this would correspond
conceptually to this class, but no hard requirement exists that a
parameterization actually does that.

We start by highlighting some general do's and don'ts for
parameterization before listing specifics for individual parameterizations.

!!! info "Parameterizations for PrimitiveEquation models only"
    The parameterizations described here can only be used for the primitive
    equation models `PrimitiveDryModel` and `PrimitiveWetModel` as the
    parameterizations are defined to act on a vertical column.
    For the 2D models `BarotropicModel` and `ShallowWaterModel` additional
    terms have to be defined as a custom forcing or drag, see [Extending SpeedyWeather](@ref).

## Define your own parameterizations

When defining a new paramerization it is required to subtype `AbstractParameterization` 
and implement the [`variables`](@ref), [`initialize!`](@ref), and [`parameterization!`](@ref) 
that define its behaviour. We'll first introduce the general idea, before 
giving a concrete example.

When defining a new parameterization with (mutable) fields
do make sure that it is constant during the model integration. While you
can and are encouraged to use the `initialize!` function to precompute
arrays (e.g. something that depends on latitude using `model.geometry.latd`)
these should not be used as work arrays on every time step of the
model integration. The reason is that the parameterization are executed
in a parallel loop over all grid points and a mutating parameterization object
would create a race condition with undefined behaviour.

Instead, you can define new prognostic and diagnostic variables for the 
parameterization to write into with the `variables` function: 

```julia 
variables(::MyParameterization) = (DiagnosticVariable(name=:flux_variable, dims=Grid2D(), units="W/m^2", desc="custom flux variable"), PrognosticVariable(name=:my_prognostic_variable, dims=Spectral2D(), units="K/s", desc="custom prognostic variable"))
```

In this example we allocate a new diagnostic variable `flux_variable` and a new prognostic variable `my_prognostic_variable` for our parameterization. The flux variable is defined as a two-dimensional
variable on our grid, and the prognostic variable is defined as a spectral variable. Three-dimensional variables are also possible by using `Grid3D` and `Spectral3D` as `dims`.

These variables are then passed to the `parameterization!` function inside of the regular `PrognosticVariables` and `DiagnosticVariables` objects. Additionally, `DiagnosticVariables` has several work 
arrays that you canreuse `diagn.grid.a` and `.b`, `.c`, `.d`. These work arrays have 
an unknown state so you should overwrite every entry and you also should not use them 
to retain information after that parameterization has been executed.

## Define the parameterization! function 

The actual parameterization computation is defined in the [`parameterization!`](@ref) function. 
This function takes in the prognostic and diagnostic variables as well as the model 
object and should compute the tendencies and fluxes that are then accumulated into 
the respective arrays. It is computed within a KernelAbstraction.jl kernel and therefore
has to be defined with GPU support in mind. This means e.g. no dynamic dispatches, only scalar 
indexing and ideally no allocations. For more details see [KernelAbstraction.jl](https://github.com/JuliaClimate/KernelAbstractions.jl) and our example parameterzations that we implemented. The signature of the function is 

```@docs 
parameterization!(ij, diagn::DiagnosticVariables, progn::PrognosticVariables, parameterization::AbstractParameterization, model_parameters)
```

## Accumulate do not overwrite

Every parameterization either computes tendencies directly or indirectly via
fluxes (upward or downward, see [Fluxes to tendencies](@ref)). Both of these are
arrays in which *every* parameterization writes into, meaning they should be
*accumulated not overwritten*. Otherwise any parameterization that executed
beforehand is effectively disabled. Hence, do

```julia
diagn.tendendies.temp_tend[k] += something_you_calculated
```

not `diagn.tendendies.temp_tend[k] = something_you_calculated` which would overwrite
any previous tendency. 

## Define generator function

After defining a (custom) parameterization it is recommended to also define
a generator function that takes in the `SpectralGrid` object
(see [How to run SpeedyWeather.jl](@ref)) as first (positional) argument,
all other arguments can then be passed on as keyword arguments with defaults
defined. Creating the default convection parameterization for example would be
```@example parameterization
using SpeedyWeather
spectral_grid = SpectralGrid(trunc=31, nlayers=8)
convection = SimplifiedBettsMiller(spectral_grid, time_scale=Hour(4))
```
Further keyword arguments can be added or omitted all together (using the default
setup), only the `spectral_grid` is required. 

## Use your parameterization

In principle, there are two different ways to use a new parameterization within SpeedyWeather.jl:

1. Define a new parameterization that replaces an existing one
2. Define a completely new parameterization and pass it on to the model constructor additional to existing ones

### Replace an existing parameterization

If you want to implement e.g. a new convection scheme, you probably want it to replace 
the already existing default convection scheme. Continuing the example from above, when 
defining a new scheme with the name matching those already present inside `PrimitiveWetModel` (see
[Keyword Arguments](https://docs.julialang.org/en/v1/manual/functions/#Keyword-Arguments)), 
we can simply pass it to the model constructor:

```@example parameterization
model = PrimitiveWetModel(spectral_grid; convection)
nothing # hide
```
otherwise we would need to write

```@example parameterization
my_convection = SimplifiedBettsMiller(spectral_grid)
model = PrimitiveWetModel(spectral_grid, convection=my_convection)
nothing # hide
```
The following is an overview of what the parameterization fields inside the
model are called. See also [Tree structure](@ref), and therein [PrimitiveDryModel](@ref)
and [PrimitiveWetModel](@ref)

- `model.vertical_diffusion`
- `model.convection`
- `model.large_scale_condensation` (`PrimitiveWetModel` only)
- `model.albedo`
- `model.optical_depth`
- `model.shortwave_radiation`
- `model.longwave_radiation`
- `model.boundary_layer_drag`
- `model.surface_condition`
- `model.surface_momentum_flux`
- `model.surface_heat_flux`
- `model.surface_humidity_flux` (`PrimitiveWetModel` only)
- `model.stochastic_physics`

Note that the parameterizations are executed in the order of the list
above. That way, for example, radiation can depend on calculations in
large-scale condensation but not vice versa (only at the next time step). 
In principle, it is possible to change this order by overwriting the default
[`get_parameterizations`](@ref) function of the model. 

### Customizing existing parameterizations

All of our existing parameterizations follow the same interface as defined before: 
A parameterization is expected to implement the following functions:

1. Define a generator function `MyParameterization(spectral_grid; kwargs...)`
2. A `initialize!(::MyParameterization, ::PrimitiveEquation)` function 
3. A `parameterization!(ij, diagn::DiagnosticVariables, progn::PrognosticVariables, ::MyParameterization, ::PrimitiveEquation)` function 
4. A `variables(::MyParameterization)` function

Our existing parameterizations also define further abstract subtypes of `AbstractParameterization` 
that can also be used to used to reuse some of the functionality of existing parameterizations. 
These include among others: `AbstractBoundaryLayers`, `AbstractLargeScaleCondensation`, `AbstractTemperatureRelaxation`, `AbstractVerticalDiffusion`,`AbstractConvection`, `AbstractAlbedo`, `AbstractOpticalDepth`, `AbstractShortwaveRadiation`, `AbstractLongwaveRadiation`, `AbstractBoundaryLayerDrag`, `AbstractSurfaceCondition`, `AbstractSurfaceMomentumFlux`, `AbstractSurfaceHeatFlux`, `AbstractSurfaceHumidityFlux`, `AbstractStochasticPhysics`, `AbstractSurfaceWind`. If you extend on these parameterizations, you can inherit from them and only need to implement the functions that are not already implemented by the parent parameterization. It's best to check the source code of the parent parameterization to see what functions are already implemented.

### Registering a new parameterization

If you want to implement a new parameterization that doesn't replace an existing one, you can register it by overwriting the `get_parameterizations(model)` function for the model your using. 
The default implementation is 

```julia 
function SpeedyWeather.get_parameterizations(model::PrimitiveWetModel)
    return (# diffusion
            vertical_diffusion = model.vertical_diffusion,
            
            # hydrological cycle
            convection = model.convection,
            large_scale_condensation = model.large_scale_condensation,
            
            # radiation
            albedo = model.albedo,
            optical_depth = model.optical_depth,
            shortwave_radiation = model.shortwave_radiation,
            longwave_radiation = model.longwave_radiation,
            
            # surface fluxes
            boundary_layer_drag = model.boundary_layer_drag,
            surface_condition = model.surface_condition,
            surface_momentum_flux = model.surface_momentum_flux,
            surface_heat_flux = model.surface_heat_flux,
            surface_humidity_flux = model.surface_humidity_flux,
            
            # stochastic physics
            stochastic_physics = model.stochastic_physics,
    )
end
```

You can add your own parameterizations to this function. The order of the parametrizations in this tuple is important as it defines the order in which they are called. 

## Example: Albedo 

Let's implement a very simple albedo parameterization as an example how to define a new parameterization. All this parametrization does is to set the albedo to a constant value over land and to linearly interpolate between the albedo of the ocean and the seaice depending on the current sea ice concentration. This albedo is written into the diagnostic variable `diagn.physics.albedo` and used later by subsequent parameterizations. 

Let's get started, 