# Differentiability and Adjoint Model

SpeedyWeather.jl is written with differentiability in mind. This means that our model is differentiable by automatic differentiation (AD). If you are interested in machine learning (ML), this means that you can integrate our model directly into your ML models without the need to first train your neural networks offline. For atmospheric modellers this means that you get an adjoint model for free which is always generated automatically, so that we don't need to maintain it separately. This allows you to calibrate SpeedyWeather.jl in a fully automatic and data-driven way.

!!! warning "Work in progress"
    The differentiability of SpeedyWeather.jl is still work in progress and some parts of this documentation might be not be always updated to the latest state. We will extend this documentation over time. Don't hesitate to contact us via GitHub issues or mail when you have questions or want to collaborate.

For the differentiability of our model we rely on [Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl). If you've used Enzyme before, just go ahead and try to differentiate the model! It should work. We have checked the correctness of the gradients extensively against a finite differences differentiation with [FiniteDifferences.jl](https://github.com/JuliaDiff/FiniteDifferences.jl/). In the following we present a simple example how we can take the gradient of a single timestep of the [primitive equation model](@ref primitive_equation_model) with respect to one of the model parameter.

!!! warning "Enzyme with Julia 1.12"
    Currently there are still some issues with Enzyme in Julia 1.11/1.12, we recommend to use Julia 1.10 for the following. 
    * Julia 1.10: All differentiation functionality tested and working
    * Julia 1.11: Basic end-to-end differentation functionality tested and working
    * Julia 1.12: Only individual functions and components and `BarotropicModel` differentiate correctly, no end-to-end differentiation of `PrimitiveWetModel` possible at the moment

    With Julia 1.12 and 1.13, use [Reactant](@ref "Differentiating with Reactant") instead.

## Differentiating through a single timestep

First we initialize the model as usual (see [Time stepping](@ref time_stepping) for
background on the time integration scheme used internally):

```julia
using SpeedyWeather, Enzyme

spectral_grid = SpectralGrid(truncation=24, nlayers=3)
model = PrimitiveWetModel(; spectral_grid)
simulation = initialize!(model)
initialize!(simulation)
run!(simulation, period=Day(10)) # spin-up the model a bit
```

Then, we get all [variables](@ref "Variables") we need from our `simulation`

```julia
(; variables, model) = simulation
(; Δt, Δt_millisec) = model.time_stepping
dt = 2Δt
```

Next, we will prepare to use Enzyme. Enzyme saves the gradient information in a shadow of the original input. For the inputs this shadow is initialized zero, whereas for the output the shadow is used as the seed of the AD. In other words, as we are doing reverse-mode AD, the shadow of the output is the value that is backpropageted by the reverse-mode AD. Ok, let's initialize everything:

```julia
dvars = make_zero(variables) # shadow for the variables
dvars.prognostic.vorticity .= 1    # seed the reverse AD
dmodel = make_zero(model) # here, we'll accumulate all parameter derivatives
```

Then, we can already do the differentiation with Enzyme

```julia
autodiff(Reverse, SpeedyWeather.timestep!, Const, Duplicated(variables, dvars), Const(dt), Duplicated(model, dmodel))
```

The derivatives are accumulated in the `dmodel` shadow. So, if we e.g. want to know the derivative with respect to the gravity constant, we just have to inspect:

```julia
dmodel.planet.gravity
```

## Differentiating longer trajectories and checkpointing 

In a very similar fashion as for the [single timestep](@ref "Differentiating through a single timestep"), we can also differentiate longer trajectories. For this, we need to use checkpointing to avoid storing all intermediate states in memory. This will also at the same time keep the compile time of the gradient still manageable. For a full example on how to do this for a sensitivity analysis, see the [sensitivity example](https://github.com/SpeedyWeather/SpeedyWeather.jl/tree/main/SpeedyWeather/test/differentiability/sensitivity_examples/).

## Differentiating with Reactant

!!! warning "Work in progress"
    Differentiating SpeedyWeather.jl with Reactant is still work in progress and experimental.

Instead of differentiating the Julia code of the model with Enzyme directly ("vanilla Enzyme",
as above), you can also run the model with [Reactant.jl](https://github.com/EnzymeAD/Reactant.jl)
(see [Reactant Support](@ref "Reactant Support")) and differentiate the program that Reactant
compiles. Enzyme then differentiates Reactant's intermediate representation (MLIR) instead of the
Julia code. In contrast to vanilla Enzyme, this also works with Julia 1.12 and Julia 1.13.

This needs a model created with the `ReactantDevice`:

```julia
using SpeedyWeather, Reactant, Enzyme
using CUDA      # needed for Reactant to compile KernelAbstractions kernels, also on CPUs

spectral_grid = SpectralGrid(architecture = SpeedyWeather.ReactantDevice(), nlayers = 8)
spectral_transform = MatrixSpectralTransform(spectral_grid)
model = PrimitiveWetModel(spectral_grid; spectral_transform)
simulation = initialize!(model)
run!(simulation, period = Day(10))      # spin-up the model a bit
initialize!(simulation, steps = 1)      # set the clock for one time step

(; variables, model) = simulation
```

As a loss we take the squared norm of the vorticity after one time step

```julia
function loss(vars, model)
    SpeedyWeather.time_step!(vars, model.time_stepping, model)
    return sum(abs2, vars.prognostic.vorticity.data)
end
```

and differentiate it with Enzyme in reverse mode, compiled with Reactant's `@jit` (or `@compile`):

```julia
dvars = make_zero(variables)
compile_options = SpeedyWeather.reactant_compile_options(; raise_first = true)
@jit compile_options = compile_options autodiff(Reverse, loss, Active, Duplicated(variables, dvars), Const(model))
```

`dvars.prognostic.vorticity` then holds the derivative of the loss with respect to the vorticity
before the time step. `SpeedyWeather.reactant_compile_options()` are the options SpeedyWeather compiles
all its programs with (`raise = true`: the KernelAbstractions kernels are raised to StableHLO as the
remainder of the model code, and the optimization patterns in `SpeedyWeather.REACTANT_EXCLUDED_PASSES`
excluded, see [GPU and Architectures](@ref)); the differentiation additionally needs `raise_first = true`.

To differentiate through several time steps, use `SpeedyWeather.time_steps!`, which runs `n`
time steps in one traced loop (so the compiled program contains the time step once, not `n`
times). How the loop is differentiated is set with its keyword arguments (as Reactant's `@trace`):
`checkpointing = false` (default) stores every time step for the reverse pass,
`checkpointing = Reactant.Periodic(k)` stores `k` checkpoints and recomputes the time steps between
them, and `checkpointing = Reactant.Binomial(budget)` uses binomial checkpointing with at most `budget`
checkpoints; `mincut = true` (default) minimizes what is stored per time step

```julia
function loss(vars, model, n, checkpointing)
    SpeedyWeather.time_steps!(Simulation(vars, model), n; checkpointing)
    return sum(abs2, vars.prognostic.vorticity.data)
end

n = 24          # time steps
checkpointing = Reactant.Periodic(4)
initialize!(simulation, steps = n)
dvars = make_zero(variables)
compile_options = SpeedyWeather.reactant_compile_options(; raise_first = true)
@jit compile_options = compile_options autodiff(Reverse, loss, Active,
    Duplicated(variables, dvars), Const(model), Const(n), Const(checkpointing))
```

!!! note "Loops of few time steps"
    Reactant unrolls traced loops of at most 4 iterations, so the checkpointing options only take
    effect for `n > 4`. `mincut = false` currently doesn't finish compiling the gradient of the loop
    (it is untested). Differentiating the loop needs Reactant v0.2.294 or newer.

!!! note "Precision on GPUs"
    On NVIDIA GPUs with TF32 support (e.g. A40, A100, H100) XLA runs Float32 matrix products in TF32
    by default, which changes the results (and gradients) considerably. For gradients in Float32
    precision see [Precision of matrix products on GPUs](@ref).

The gradients with Reactant (of one time step, and of several time steps with each of the
checkpointing options above) are compared to vanilla Enzyme on the CPU for the same state:
in our CI for the `BarotropicModel` (`SpeedyWeather/test/reactant/differentation.jl`) and,
as this is too slow for the CI, for the `PrimitiveWetModel` with the script
`SpeedyWeather/test/reactant/differentiation_primitive_wet.jl` (the CI only checks that the gradient of
several time steps of the `PrimitiveWetModel` compiles and is finite). Note that this comparison needs Julia 1.10 or 1.11 because of vanilla Enzyme, the differentiation with Reactant itself doesn't.

## Correctness 

The correctness of gradients is tested in our CI in `SpeedyWeather/test/differentiability` by comparing the results of the reverse mode to analytic results for function where it's available, against the forward mode and against finite difference when no analytic results are possible. These differentiability tests also feature many scripts using the differentiabilty of SpeedyWeather and might be a good resource if you want to apply it, too!

## Parameter handling

SpeedyWeather also provides automated parameter handling for all [models](@ref "Models") and subcomponents via an extension of [ModelParameters.jl](https://github.com/rafaqz/ModelParameters.jl). Parameters can be automatically collected via the `parameters` method:

```julia
spectral_grid = SpectralGrid(truncation=24, nlayers=1) 
model = BarotropicModel(; spectral_grid)
params = parameters(model)

# output (truncated)
Parameters:
┌───────┬──────────────────────────┬──────────┬────────────┬──────────────────────────────┬─────────────────────────────────────┬──────────────────────────────────────────────────────────────────────────────────────────┐
│   idx │                fieldname │      val │  component │               component_type │                              bounds │                                                                                     desc │
│ Int64 │                   Symbol │  Float32 │     Symbol │                     DataType │ IntervalSets.TypedEndpointsInterval │                                                                                   String │
├───────┼──────────────────────────┼──────────┼────────────┼──────────────────────────────┼─────────────────────────────────────┼──────────────────────────────────────────────────────────────────────────────────────────┤
│     1 │                 rotation │  7.29e-5 │     planet │               Earth{Float32} │       -Inf .. Inf (open) (RealLine) │                                            angular frequency of Earth's rotation [rad/s] │
│     2 │                  gravity │     9.81 │     planet │               Earth{Float32} │ 0.0 .. Inf (closed-open) (HalfLine) │                                                       gravitational acceleration [m/s^2] │
│     3 │               axial_tilt │     23.4 │     planet │               Earth{Float32} │                           -90 .. 90 │                                                angle [˚] rotation axis tilt wrt to orbit │
│     4 │           solar_constant │   1365.0 │     planet │               Earth{Float32} │ 0.0 .. Inf (closed-open) (HalfLine) │                                    Total solar irradiance at the distance of 1 AU [W/m²] │
│     5 │         mol_mass_dry_air │  28.9649 │ atmosphere │     EarthAtmosphere{Float32} │        0.0 .. Inf (open) (HalfLine) │                                                            molar mass of dry air [g/mol] │
│     6 │           mol_mass_vapor │  18.0153 │ atmosphere │     EarthAtmosphere{Float32} │        0.0 .. Inf (open) (HalfLine) │                                                        molar mass of water vapor [g/mol] │
│     7 │            heat_capacity │   1004.0 │ atmosphere │     EarthAtmosphere{Float32} │ 0.0 .. Inf (closed-open) (HalfLine) │                                           specific heat at constant pressure cₚ [J/K/kg] │
│     8 │                  R_vapor │  461.524 │ atmosphere │     EarthAtmosphere{Float32} │        0.0 .. Inf (open) (HalfLine) │                                           specific gas constant for water vapor [J/kg/K] │
│     9 │                mol_ratio │  0.62197 │ atmosphere │     EarthAtmosphere{Float32} │        0.0 .. Inf (open) (HalfLine) │                        Ratio of gas constants: dry air / water vapor, often called ε [1] │
│    10 │              μ_virt_temp │ 0.607794 │ atmosphere │     EarthAtmosphere{Float32} │        0.0 .. Inf (open) (HalfLine) │ Virtual temperature Tᵥ calculation, Tᵥ = T(1 + μ*q), humidity q, absolute tempereature T │
│    11 │                        κ │ 0.285911 │ atmosphere │     EarthAtmosphere{Float32} │       -Inf .. Inf (open) (RealLine) │                                         = R_dry/cₚ, gas const for air over heat capacity │
│    12 │            water_density │   1000.0 │ atmosphere │     EarthAtmosphere{Float32} │        0.0 .. Inf (open) (HalfLine) │                                                                    water density [kg/m³] │
│    13 │ latent_heat_condensation │  2.501e6 │ atmosphere │     EarthAtmosphere{Float32} │ 0.0 .. Inf (closed-open) (HalfLine) │                                                       latent heat of condensation [J/kg] │
│    14 │  latent_heat_sublimation │  2.801e6 │ atmosphere │     EarthAtmosphere{Float32} │ 0.0 .. Inf (closed-open) (HalfLine) │                                                        latent heat of sublimation [J/kg] │
│    15 │       reference_pressure │ 100000.0 │ atmosphere │     EarthAtmosphere{Float32} │        0.0 .. Inf (open) (HalfLine) │                                                          surface reference pressure [Pa] │
│    16 │    reference_temperature │    288.0 │ atmosphere │     EarthAtmosphere{Float32} │ 0.0 .. Inf (closed-open) (HalfLine) │                                                        surface reference temperature [K] │
│    17 │         moist_lapse_rate │    0.005 │ atmosphere │     EarthAtmosphere{Float32} │       -Inf .. Inf (open) (RealLine) │                                   reference moist-adiabatic temperature lapse rate [K/m] │
│    18 │           dry_lapse_rate │   0.0098 │ atmosphere │     EarthAtmosphere{Float32} │       -Inf .. Inf (open) (RealLine) │                                     reference dry-adiabatic temperature lapse rate [K/m] │
│    19 │          layer_thickness │   8500.0 │ atmosphere │     EarthAtmosphere{Float32} │        0.0 .. Inf (open) (HalfLine) │                                          layer thickness for the shallow water model [m] │
│    20 │                 strength │  3.0e-12 │    forcing │      KolmogorovFlow{Float32} │       -Inf .. Inf (open) (RealLine) │                                                      [OPTION] Strength of forcing [1/s²] │
│    21 │               wavenumber │      8.0 │    forcing │      KolmogorovFlow{Float32} │        0.0 .. Inf (open) (HalfLine) │                    [OPTION] Wavenumber of forcing in meridional direction (pole to pole) │
│    22 │                        c │   1.0e-7 │       drag │ LinearVorticityDrag{Float32} │ 0.0 .. Inf (closed-open) (HalfLine) │                                                          [OPTION] drag coefficient [1/s] │
└───────┴──────────────────────────┴──────────┴────────────┴──────────────────────────────┴─────────────────────────────────────┴──────────────────────────────────────────────────────────────────────────────────────────┘
```

The returned `ParameterTable` object implements the `Model` interface from ModelParmaeters.jl which allows you to interact with the parameter metadata in tabular form. For example, we could extract the values of the parameters with `params[:,:val]` or the bounds with `params[:,:bounds]`. Subsets of parameters can also be extracted by indexing `params` with one or more `String` variable names (or prefxes), e.g:

```julia
param_subset = params[["planet.gravity", "atmosphere.heat_capacity"]]

# output (truncated)

Parameters:
┌───────┬───────────────┬─────────┬────────────┬──────────────────────────┬───────────────────────────────────────┬────────────────────────────────────────────────┐
│   idx │     fieldname │     val │  component │           component_type │                                bounds │                                           desc │
│ Int64 │        Symbol │ Float32 │     Symbol │                 DataType │ DomainSets.HalfLine{Float64, :closed} │                                         String │
├───────┼───────────────┼─────────┼────────────┼──────────────────────────┼───────────────────────────────────────┼────────────────────────────────────────────────┤
│     1 │       gravity │    9.81 │     planet │           Earth{Float32} │   0.0 .. Inf (closed-open) (HalfLine) │             gravitational acceleration [m/s^2] │
│     2 │ heat_capacity │  1004.0 │ atmosphere │ EarthAtmosphere{Float32} │   0.0 .. Inf (closed-open) (HalfLine) │ specific heat at constant pressure cₚ [J/K/kg] │
└───────┴───────────────┴─────────┴────────────┴──────────────────────────┴───────────────────────────────────────┴────────────────────────────────────────────────┘
```

## Vectorizing parameters

Many sensitivity analysis, optimization, or uncertainty quantification algorithms require the parameters to be supplied as one or more vectors of values. `ParameterTable` provides a dispatch for `Base.vec` that flattens the model parameters into a [ComponentVector](https://docs.sciml.ai/ComponentArrays/stable/quickstart/):

```julia
param_vec = vec(params)

# output

ComponentVector{Float32}(planet = (rotation = 7.29f-5, gravity = 9.81f0, axial_tilt = 23.4f0, solar_constant = 1365.0f0), atmosphere = (mol_mass_dry_air = 28.9649f0, mol_mass_vapor = 18.0153f0, heat_capacity = 1004.0f0, R_vapor = 461.52438f0, mol_ratio = 0.62197006f0, μ_virt_temp = 0.60779446f0, κ = 0.2859107f0, water_density = 1000.0f0, latent_heat_condensation = 2.501f6, latent_heat_sublimation = 2.801f6, reference_pressure = 100000.0f0, reference_temperature = 288.0f0, moist_lapse_rate = 0.005f0, dry_lapse_rate = 0.0098f0, layer_thickness = 8500.0f0), forcing = (strength = 3.0f-12, wavenumber = 8.0f0), drag = (c = 1.0f-7))
```

`ComponentVector`s behave like normal `Array`s but additionally allow you to access the components following the original nested structure in the model, e.g. `param_vec.planet.solar_constant` will extract the solar constant parameter from the `Earth` component.

We can use the resulting parameter vector to calculate sensitivities over a single time step:

```julia
initialize!(simulation)
run!(simulation, period=Day(10))
(; Δt) = simulation.model.time_stepping
ps = parameters(model)
pvec = vec(ps)
dp = zero(pvec)
dvars = make_zero(variables) # shadow for the variables
dvars.prognostic.vorticity .= 1    # seed the reverse AD

function timestep_with_new_params!(vars, dt, model, p)
    new_model = SpeedyWeather.reconstruct(model, p)
    SpeedyWeather.timestep!(vars, dt, new_model)
    return nothing
end

autodiff(Reverse, timestep_with_new_params!, Const, Duplicated(variables, dvars), Const(dt), Duplicated(model, make_zero(model)))
```

Note, however, that a full sensitivity analysis over long integration periods is computationally much more demanding, and is something that we are currently working on.

Stay tuned!
