# Differentiability and Adjoint Model

SpeedyWeather.jl is written with differentiability in mind. This means that our model is differentiable by automatic differentiation (AD). If you are interested in machine learning (ML), this means that you can integrate our model directly into your ML models without the need to train your ANNs seperately first. For atmospheric modellers this means that you get an adjoint model for free which is always generated automatically, so that we don't need to maintain it seperatly. So, you can calibrate SpeedyWeather.jl in a fully automatic, data-driven way. 

!!! warn Work in progress
    The differentiability of SpeedyWeather.jl is still work in progress and some parts of this documentation might be not be always updated to the latest state. We will extend this documentation over time. Don't hesitate to contact us via GitHub issues or mail when you have questions or want to colloborate.

For the differentiability of our model we rely on [Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl). If you've used Enzyme before, just go ahead and try to differentiate the model! It should work. We have checked the correctness of the gradients extensively against a finite differences differentiation with [FiniteDifferences.jl](https://github.com/JuliaDiff/FiniteDifferences.jl/). In the following we present a simple example how we can take the gradient of a single timestep of the primitive equation model with respect to one of the model parameter. 

!!! warn Enzyme with Julia 1.11
    Currently there are still some issues with Enzyme in Julia 1.11, we recommend to use Julia 1.10 for the following

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

Doing a full sensitivity analysis through a long integration is computationally much more demanding, and is something that we are currently working on. 
