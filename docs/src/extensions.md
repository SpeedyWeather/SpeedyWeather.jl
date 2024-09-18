# Extending SpeedyWeather

Generally, SpeedyWeather is built in a very modular, extensible way.
While that sounds fantastic in general, it does not save you from understanding
its [modular logic](@ref logic) before you can extend SpeedyWeather.jl easily yourself.
We highly recommend you to read the following sections if you would like to
extend SpeedyWeather in some way, but it also gives you a good understanding
of how we build SpeedyWeather in the first place. Because in the end there
is no difference between internally or externally defined model components.
Having said that, there is a question of the [Scope of variables](@ref)
meaning that some functions or types require its module to be explicitly named,
like `SpeedyWeather.some_function` instead of just `some_function`. But that's
really it.

Before and especially after reading this section you are welcome to
[raise an issue](https://github.com/SpeedyWeather/SpeedyWeather.jl/issues)
about whatever you would like to do with SpeedyWeather. We are happy to help.

## [SpeedyWeather's modular logic](@id logic)

Almost every component in SpeedyWeather is implemented in three steps:

1. Define a new type,
2. define its initialization,
3. extend a function that defines what it does.

To be a bit more explicit (Julia always encourages you to think more abstractly,
which can be difficult to get started...) we will use the example of defining
a new forcing for the `Barotropic` or `ShallowWater` models. But the concept
is the same whether you want to define a forcing, a drag, or a new parameterization
for the primitive equations, etc. For the latter, see [Parameterizations](@ref)
In general, you can define a new component also just in a notebook or in the Julia REPL,
you do not have to branch off from the repository and write directly into it.
However, note the [Scope of variables](@ref) if you define a component externally.

To define a new forcing type, at the most basic level you would do
```@example extending
using SpeedyWeather

struct MyForcing{NF} <: SpeedyWeather.AbstractForcing
    # define some parameters and work arrays here
    a::NF
    v::Vector{NF}
end
```
In Julia this introduces a new (so-called compound) type that is a subtype of `AbstractForcing`,
we have a bunch of these abstract super types defined (see [Abstract model components](@ref))
and you want to piggy-back on them because of multiple-dispatch. This new type could also be 
a `mutable struct`, could have keywords defined with `@kwdef` and can
also be parametric with respect to the number format `NF` or grid, but let's skip
those details for now. Conceptually you include into the type any parameters
(example the float `a` here) that you may need and especially those that you want to change
(ideally not work arrays, see discussion in [Use `ColumnVariables` work arrays](@ref)).
This type will get a user-facing interface so that one can quickly create a new
forcing but with altered parameters. Generally you should also include any
kind of precomputed arrays (here a vector `v`). For example, you want to
apply your forcing only in certain parts of the globe? Then you probably want
to define a mask here that somehow includes the information of your region.
For a more concrete example see [Custom forcing and drag](@ref).

To define the new type's initialization, at the most basic level you need
to extend the `initialize!` function for this new type. A dummy example:
```@example extending
function initialize!(forcing::MyForcing, model::AbstractModel)
    # fill in/change any fields of your new forcing here
    forcing.v[1] = 1
    # you can use information from other model components too
    forcing.v[2] = model.planet.gravity
end
```
This function is called _once_ during model initialisation (which is in fact
just the initialisation of all its components, like the forcing here)
and it allows you to precompute values or arrays also based on
parameters of other model components. Like in this example, we want to
use the gravity that is defined in `model.planet`. If you need a value
for gravity in your forcing you could add a gravity field therein, but
then if you change `planet.gravity` this change would not propagate
into your forcing! Another example would be to use `model.geometry.coslat`
if you need to use the cosine of latitude for some precomputation,
which, however, depends on the resolution and so should not be hardcoded
into your forcing.

As the last step we have to extend the `forcing!` function which is the function
that is called on _every_ step of the time integration. This new method
for `forcing!` needs to have the following function signature
```@example extending
function forcing!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    forcing::MyForcing,
    model::AbstractModel,
    lf::Integer,
)
    # whatever the forcing is supposed to do, in the end you want
    # to write into the tendency fields
    diagn.tendencies.u_tend_grid = forcing.a
    diagn.tendencies.v_tend_grid = forcing.a
    diagn.tendencies.vor_tend = forcing.a
end
```
`DiagnosticVariables` is the type of the first argument, because it contains
the tendencies you will want to change, so this is supposed to be read and write.
The other arguments should be treated read-only. You can make use of anything else
in `model`, but often we unpack the model in a function barrier (which can help with
type inference and therefore performance). But let's skip that detail for now.
Generally, try to precompute what you can in
`initialize!`. For the forcing you will need to force the velocities `u, v` in
grid-point space or the vorticity `vor`, divergence `div` in spectral space.
This is not a constrain in most applications we came across, but in case it
is in yours please reach out.

## Scope of variables

The above (conceptual!) examples leave out some of the details, particularly around the
scope of variables when you want to define a new forcing interactively inside a notebook
or the REPL (which is actually the recommended way!!). To respect the scope of
variables, a bunch of functions will need their module to be explicit specified.
In general, you should be familiar with Julia's 
[scope of variables](https://docs.julialang.org/en/v1/manual/variables-and-scoping/)
logic.

The `initialize!` function is a function *inside* the SpeedyWeather module,
as we want to define a new method for it *outside* that can be called *inside*
we actually need to write
```@example extending
function SpeedyWeather.initialize!(forcing::MyForcing, model::SpeedyWeather.AbstractModel)
    # how to initialize it
end
```
And similar for `SpeedyWeather.forcing!`.

You also probably want to make use of functions that are already defined inside
SpeedyWeather or its submodules `SpeedyTransforms`, or `RingGrids`. If something
does not seem to be defined, although you can see it in the documentation or
directly in the code, you probably need to specify its module too! Alternatively,
note that you can also always do `import SpeedWeather: AbstractModel` to bring
a given variable into global scope which removes the necessity to write
`SpeedyWeather.AbstractModel`.

## Abstract model components

You may wonder which abstract model components there are, you can always check this
with
```@example extending
using InteractiveUtils # hide
subtypes(SpeedyWeather.AbstractModelComponent)
```

we illustrated the modular logic here using `AbstractForcing` and
`AbstractDrag` is very similar. However, other model components also
largely follow [SpeedyWeather's modular logic](@ref logic) as
for example outlined in [Defining a callback](@ref) or 
[Defining a new orography type](@ref). If you do not find much
documentation about a new custom type where you would like
extend SpeedyWeather's functionality it is probably because
we have not experimented much with this either. But that
does not mean it is not possible. Just reach out by
creating an issue in this case.

Similarly, `AbstractParameterization` has several
subtypes that define conceptual classes of parameterizations, namely
```@example extending
subtypes(SpeedyWeather.AbstractParameterization)
```

but these are discussed in more detail in [Parameterizations](@ref).
For a more concrete example of how to define a new forcing for the 2D models,
see [Custom forcing and drag](@ref).