# Differentiability and Adjoint Model

SpeedyWeather.jl is written with differentiability in mind. This means that our model is differentiable by automatic differentiation (AD). If you are interested in machine learning (ML), this means that you can integrate our model directly into your ML models without the need to first train your neural networks offline. For atmospheric modellers this means that you get an adjoint model for free which is always generated automatically, so that we don't need to maintain it separately. This allows you to calibrate SpeedyWeather.jl in a fully automatic and data-driven way.

!!! warn Work in progress
    The differentiability of SpeedyWeather.jl is still work in progress and some parts of this documentation might be not be always updated to the latest state. We will extend this documentation over time. Don't hesitate to contact us via GitHub issues or mail when you have questions or want to colloborate.

For the differentiability of our model we rely on [Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl). If you've used Enzyme before, just go ahead and try to differentiate the model! It should work. We have checked the correctness of the gradients extensively against a finite differences differentiation with [FiniteDifferences.jl](https://github.com/JuliaDiff/FiniteDifferences.jl/). In the following we present a simple example how we can take the gradient of a single timestep of the primitive equation model with respect to one of the model parameter. 

!!! warn Enzyme with Julia 1.11
    Currently there are still some issues with Enzyme in Julia 1.11, we recommend to use Julia 1.10 for the following

## Differentiating through a single timestep

First we initialize the model as usual: 

```julia
using SpeedyWeather, Enzyme 

spectral_grid = SpectralGrid(trunc=23, nlayers=3)           
model = PrimitiveWetModel(; spectral_grid) 
simulation = initialize!(model)  
initialize!(simulation)
run!(simulation, period=Day(10)) # spin-up the model a bit
```

Then, we get all variables we need from our `simulation`

```julia
(; prognostic_variables, diagnostic_variables, model) = simulation
(; Δt, Δt_millisec) = model.time_stepping
dt = 2Δt

progn = prognostic_variables
diagn = diagnostic_variables
```

Next, we will prepare to use Enzyme. Enzyme saves the gradient information in a shadow of the original input. For the inputs this shadow is initialized zero, whereas for the output the shadow is used as the seed of the AD. In other words, as we are doing reverse-mode AD, the shadow of the output is the value that is backpropageted by the reverse-mode AD. Ok, let's initialize everything: 

```julia
dprogn = one(progn) # shadow for the progn values 
ddiagn = make_zero(diagn) # shadow for the diagn values 
dmodel = make_zero(model) # here, we'll accumulate all parameter derivatives 
```

Then, we can already do the differentiation with Enzyme

```julia
autodiff(Reverse, SpeedyWeather.timestep!, Const, Duplicated(progn, dprogn), Duplicated(diagn, ddiagn), Const(dt), Duplicated(model, dmodel))
```

The derivitaves are accumulated in the `dmodel` shadow. So, if we e.g. want to know the derivative with respect to the gravity constant, we just have to inspect: 

```julia 
dmodel.planet.gravity 
```
## Parameter handling

SpeedyWeather also provides automated parameter handling for all models and subcomponents via an extension of [ModelParameters.jl](https://github.com/rafaqz/ModelParameters.jl). Parameters can be automatically collected via the `parameters` method:

```julia
spectral_grid = SpectralGrid(trunc=23, nlayers=1) 
model = Barotropic(; spectral_grid)
params = parameters(model)

# output (truncated)
Parameters:
┌───────┬──────────────────────────┬──────────┬────────────┬──────────────────────────────┬─────────────────────────────────────┬──────────────────────────────────────────────────────────────────────────────────────────┐       
│   idx │                fieldname │      val │  component │               componentttype │                              bounds │                                                                                     desc │       
│ Int64 │                   Symbol │  Float32 │     Symbol │                     DataType │ IntervalSets.TypedEndpointsInterval │                                                                                   String │       
├───────┼──────────────────────────┼──────────┼────────────┼──────────────────────────────┼─────────────────────────────────────┼──────────────────────────────────────────────────────────────────────────────────────────┤       
│     1 │                 rotation │  7.29e-5 │     planet │               Earth{Float32} │       -Inf .. Inf (open) (RealLine) │                                            angular frequency of Earth's rotation [rad/s] │       
│     2 │                  gravity │     9.81 │     planet │               Earth{Float32} │ 0.0 .. Inf (closed-open) (HalfLine) │                                                       gravitational acceleration [m/s^2] │       
│     3 │               axial_tilt │     23.4 │     planet │               Earth{Float32} │                           -90 .. 90 │                                                angle [˚] rotation axis tilt wrt to orbit │       
│     4 │           solar_constant │   1365.0 │     planet │               Earth{Float32} │ 0.0 .. Inf (closed-open) (HalfLine) │                                    Total solar irradiance at the distance of 1 AU [W/m²] │       
│     5 │         mol_mass_dry_air │  28.9649 │ atmosphere │     EarthAtmosphere{Float32} │        0.0 .. Inf (open) (HalfLine) │                                                            molar mass of dry air [g/mol] │       
│     6 │          mol_mass_vapour │  18.0153 │ atmosphere │     EarthAtmosphere{Float32} │        0.0 .. Inf (open) (HalfLine) │                                                       molar mass of water vapour [g/mol] │       
│     7 │            heat_capacity │   1004.0 │ atmosphere │     EarthAtmosphere{Float32} │ 0.0 .. Inf (closed-open) (HalfLine) │                                           specific heat at constant pressure cₚ [J/K/kg] │       
│     8 │                 R_vapour │  461.524 │ atmosphere │     EarthAtmosphere{Float32} │        0.0 .. Inf (open) (HalfLine) │                                          specific gas constant for water vapour [J/kg/K] │       
│     9 │                mol_ratio │  0.62197 │ atmosphere │     EarthAtmosphere{Float32} │        0.0 .. Inf (open) (HalfLine) │                       Ratio of gas constants: dry air / water vapour, often called ε [1] │       
│    10 │              μ_virt_temp │ 0.607794 │ atmosphere │     EarthAtmosphere{Float32} │        0.0 .. Inf (open) (HalfLine) │ Virtual temperature Tᵥ calculation, Tᵥ = T(1 + μ*q), humidity q, absolute tempereature T │       
│    11 │                        κ │ 0.285911 │ atmosphere │     EarthAtmosphere{Float32} │       -Inf .. Inf (open) (RealLine) │                                         = R_dry/cₚ, gas const for air over heat capacity │       
│    12 │            water_density │   1000.0 │ atmosphere │     EarthAtmosphere{Float32} │        0.0 .. Inf (open) (HalfLine) │                                                                    water density [kg/m³] │       
│    13 │ latent_heat_condensation │  2.501e6 │ atmosphere │     EarthAtmosphere{Float32} │ 0.0 .. Inf (closed-open) (HalfLine) │                                                       latent heat of condensation [J/kg] │       
│    14 │  latent_heat_sublimation │  2.801e6 │ atmosphere │     EarthAtmosphere{Float32} │ 0.0 .. Inf (closed-open) (HalfLine) │                                                        latent heat of sublimation [J/kg] │       
│    15 │                 pres_ref │ 100000.0 │ atmosphere │     EarthAtmosphere{Float32} │        0.0 .. Inf (open) (HalfLine) │                                                          surface reference pressure [Pa] │       
│    16 │                 temp_ref │    288.0 │ atmosphere │     EarthAtmosphere{Float32} │ 0.0 .. Inf (closed-open) (HalfLine) │                                                        surface reference temperature [K] │       
│    17 │         moist_lapse_rate │    0.005 │ atmosphere │     EarthAtmosphere{Float32} │       -Inf .. Inf (open) (RealLine) │                                   reference moist-adiabatic temperature lapse rate [K/m] │       
│    18 │           dry_lapse_rate │   0.0098 │ atmosphere │     EarthAtmosphere{Float32} │       -Inf .. Inf (open) (RealLine) │                                     reference dry-adiabatic temperature lapse rate [K/m] │       
│    19 │          layer_thickness │   8500.0 │ atmosphere │     EarthAtmosphere{Float32} │        0.0 .. Inf (open) (HalfLine) │                                          layer thickness for the shallow water model [m] │       
│    20 │                 strength │  3.0e-12 │    forcing │      KolmogorovFlow{Float32} │       -Inf .. Inf (open) (RealLine) │                                                      [OPTION] Strength of forcing [1/s²] │       
│    21 │               wavenumber │      8.0 │    forcing │      KolmogorovFlow{Float32} │        0.0 .. Inf (open) (HalfLine) │                    [OPTION] Wavenumber of forcing in meridional direction (pole to pole) │       
│    22 │                        c │   1.0e-7 │       drag │ LinearVorticityDrag{Float32} │ 0.0 .. Inf (closed-open) (HalfLine) │                                                          [OPTION] drag coefficient [1/s] │       
└───────┴──────────────────────────┴──────────┴────────────┴──────────────────────────────┴─────────────────────────────────────┴──────────────────────────────────────────────────────────────────────────────────────────┘
```

The returned `SpeedyParams` object implements the `Model` interface from ModelParmaeters.jl which allows you to interact with the parameter metadata in tablar form. For example, we could extract the values of the parameters with `params[:,:val]` or the bounds with `params[:,:bounds]`. Subsets of parameters can also be extracted by indexing `params` with one or more `String` variable names (or prefxes), e.g:

```julia
param_subset = params[["planet.gravity", "atmosphere.heat_capacity"]]

# output (truncated)

Parameters:
┌───────┬───────────────┬─────────┬────────────┬──────────────────────────┬───────────────────────────────────────┬────────────────────────────────────────────────┐
│   idx │     fieldname │     val │  component │           componentttype │                                bounds │                                           desc │
│ Int64 │        Symbol │ Float32 │     Symbol │                 DataType │ DomainSets.HalfLine{Float64, :closed} │                                         String │
├───────┼───────────────┼─────────┼────────────┼──────────────────────────┼───────────────────────────────────────┼────────────────────────────────────────────────┤
│     1 │       gravity │    9.81 │     planet │           Earth{Float32} │   0.0 .. Inf (closed-open) (HalfLine) │             gravitational acceleration [m/s^2] │
│     2 │ heat_capacity │  1004.0 │ atmosphere │ EarthAtmosphere{Float32} │   0.0 .. Inf (closed-open) (HalfLine) │ specific heat at constant pressure cₚ [J/K/kg] │
└───────┴───────────────┴─────────┴────────────┴──────────────────────────┴───────────────────────────────────────┴────────────────────────────────────────────────┘
```

## Vectorizing parameters 

Many sensitivity analysis, optimization, or uncertainty quantification algorithms require the parameters to be supplied as one or more vectors of values. `SpeedyParams` provides a dispatch for `Base.vec` that flattens the model parameters into a [ComponentVector](https://docs.sciml.ai/ComponentArrays/stable/quickstart/):

```julia
param_vec = vec(params)

# output

ComponentVector{Float32}(planet = (rotation = 7.29f-5, gravity = 9.81f0, axial_tilt = 23.4f0, solar_constant = 1365.0f0), atmosphere = (mol_mass_dry_air = 28.9649f0, mol_mass_vapour = 18.0153f0, heat_capacity = 1004.0f0, R_vapour = 461.52438f0, mol_ratio = 0.62197006f0, μ_virt_temp = 0.60779446f0, κ = 0.2859107f0, water_density = 1000.0f0, latent_heat_condensation = 2.501f6, latent_heat_sublimation = 2.801f6, pres_ref = 100000.0f0, temp_ref = 288.0f0, moist_lapse_rate = 0.005f0, dry_lapse_rate = 0.0098f0, layer_thickness = 8500.0f0), forcing = (strength = 3.0f-12, wavenumber = 8.0f0), drag = (c = 1.0f-7))
```

`ComponentVector`s behave like normal `Array`s but additionally allow you to access the components following the original nested structure in the model, e.g. `param_vec.planet.solar_constant` will extract the solar constant parameter from the `Earth` component.

We can use the resulting parameter vector to calculate sensitivities over a single time step:

```julia
initialize!(simulation)
run!(simulation, period=Day(10))
(; Δt, Δt_sec) = simulation.model.time_stepping
ps = parameters(model)
pvec = vec(ps)
dp = zero(pvec)
dprogn = one(progn) # shadow for the prognostic variabels
ddiagn = make_zero(diagn) # shadow for the diagnostic variables

function timestep_with_new_params!(progn, diagn, dt, model, p)
    new_model = SpeedyWeather.reconstruct(model, p)
    SpeedyWeather.timestep!(progn, diagn, dt, new_model)
    return nothing
end

autodiff(Reverse, timestep_with_new_params!, Const, Duplicated(progn, dprogn), Duplicated(diagn, ddiagn), Const(dt), Duplicated(model, make_zero(model)))
```

Note, however, that a full sensitivity analysis over long integration periods is computationally much more demanding, and is something that we are currently working on.

Stay tuned!
