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
for the primitive equations, etc. In general, you can define a new component
also just in a notebook or in the Julia REPL, you do not have to branch off
from the repository and write directly into it. However, note the [Scope of variables](@ref)
if you define a component externally.

To define a new forcing type, at the most basic level you would do
```julia
struct MyForcing <: AbstractForcing
    # define some parameters and work arrays here
    a::Float64
    v::Vector{Float64}
end
```
In Julia this introduces a new (so-called compound) type that is a subtype of `AbstractForcing`,
we have a bunch of these abstract super types defined and you want to
piggy-back on them because of multiple-dispatch. This new type could also be 
a `mutable struct`, could have keywords defined with `Base.@kwdef` and can
also be parametric with respect to the number format or grid, but let's skip
those details for now. Conceptually you include into the type any parameters
(example the float `a` here) that you may need and especially those that you want to change.
This type will get a user-facing interface so that one can quickly create a new
forcing but with altered parameters. Generally you should also include any
kind of precomputed or work arrays (here a vector `v`). For example, you want to
apply your forcing only in certain parts of the globe? Then you probably want
to define a mask here that somehow includes the information of your region.
For a more concrete example see [Add custom forcing](@ref).

To define the new type's initialization, at the most basic level you need
to extend the `initialize!` function for this new type
```julia
function initialize!(forcing::MyForcing,model::ModelSetup)
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
```julia
function forcing!(  diagn::DiagnosticVariablesLayer,
                    forcing::MyForcing,
                    time::DateTime,
                    model::ModelSetup)
    # whatever the forcing is supposed to do, in the end you want
    # to write into the tendency fields
    diagn.tendencies.u_tend_grid[1] = forcing.a
    diagn.tendencies.v_tend_grid[1] = forcing.a
    diagn.tendencies.vor_tend[1] = forcing.a
end
```
`DiagnosticVariablesLayer` is the type of the first argument, because it contains
the tendencies you will want to change, so this is supposed to be read and write.
The other arguments should be treated read-only. You make use of the `time`
or anything else in the `model`, but the latter likely comes with a performance
penalty which is why we often unpack the model in a function barrier. But let's
skip that detail for now. Generally, try to precompute what you can in
`initialize!`. For the forcing you will need to force the velocities `u,v` in
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

Defining a new type for example as a subtype of AbstractForcing, you will actually
want to do
```julia
struct MyForcing <: SpeedyWeather.AbstractForcing
```
Because `AbstractForcing` is defined inside SpeedyWeather, but we do not export
all functions/types in order to not flood your global scope to avoid
naming conflicts. Rule of thumb: If you get an error hinting at something is
not defined, make sure you are in the right scope!

Then the `initialize!` function is a function *inside* the SpeedyWeather module,
as we want to define a new method for it *outside* that can be called *inside*
we actually need to write
```julia
function SpeedyWeather.initialize!(forcing::MyForcing,model::SpeedyWeather.ModelSetup)
```
And similar for `SpeedyWeather.forcing!`.

You also probably want to make use of functions that are already defined inside
SpeedyWeather or its submodules `SpeedyTransforms`, or `RingGrids`. If something
does not seem to be defined, although you can see it in the documentation or
directly in the code, you probably need to specify its module too!

## Add custom forcing

The following example is a bit more concrete than the previous conceptual example,
but we try to add a few more details that are important, or you at least should
be aware of it. In this example we want to add a `StochasticStirring` forcing
as defined in [Vallis et al., 2004](https://doi.org/10.1175/1520-0469(2004)061%3C0264:AMASDM%3E2.0.CO;2)

```math
\begin{aligned}
\frac{\partial \zeta}{\partial t} &+ \nabla \cdot (\mathbf{u}(\zeta + f)) =
S - r\zeta - \nu\nabla^{4}\zeta \\ 
S_{l,m}^i &= A(1-\exp(-2\tfrac{\Delta t}{\tau}))Q^i_{l,m} + \exp(-\tfrac{dt}{\tau})S_{l,m}^{i-1} \\
\end{aligned}
```

So there is a term `S` that is supposed to force the vorticity equation in the
[Barotropic vorticity model]. However, this term is also stochastically
evolving in time, meaning we have to store the previous time steps, `i-1`,
in spectral space, because that's where the forcing is defined: for degree
`l` and order `m` of the spherical harmonics. `A` is a real amplitude.
`\Delta t` the time step of the model, `\tau` the decorrelation time scale
of the stochastic process. `Q` is for every spherical harmonic a complex random uniform
number in `[-1,1]` in both its real and imaginary components. 
So we actually define our `StochasticStirring` forcing as follows and
will explain the details in second
```julia
Base.@kwdef struct StochasticStirring{NF} <: SpeedyWeather.AbstractForcing{NF}
    
    # DIMENSIONS from SpectralGrid
    "Spectral resolution as max degree of spherical harmonics"
    trunc::Int
    
    "Number of latitude rings, used for latitudinal mask"
    nlat::Int

    
    # OPTIONS
    "Decorrelation time scale τ [days]"
    decorrelation_time::Float64 = 2

    "Stirring strength A [1/s²]"
    strength::Float64 = 1e-11

    "Stirring latitude [˚N]"
    latitude::Float64 = 45

    "Stirring width [˚]"
    width::Float64 = 24

    
    # TO BE INITIALISED
    "Stochastic stirring term S"
    S::LowerTriangularMatrix{Complex{NF}} = zeros(LowerTriangularMatrix{Complex{NF}},trunc+2,trunc+1)
    
    "a = A*sqrt(1 - exp(-2dt/τ)), the noise factor times the stirring strength [1/s²]"
    a::Base.RefValue{NF} = Ref(zero(NF))
        
    "b = exp(-dt/τ), the auto-regressive factor [1]"
    b::Base.RefValue{NF} = Ref(zero(NF))
        
    "Latitudinal mask, confined to mid-latitude storm track by default [1]"
    lat_mask::Vector{NF} = zeros(NF,nlat)
end
```


## Add custom drag