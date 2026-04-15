# Other output

Next to the [NetCDF output](@ref) there's other ways to output data from
the model, these are described below.

## JLD2 Output 

As an alternative to the [NetCDF output](@ref), it is also possible to directly
output the `Variables` (or one subgroup of it) to a JLD2 file.
This might be interesting if you are really interested in the model internals,
or also for some machine learning tasks. However, this option doesn't feature
any possibilites to regrid or select variables, and it comes with the usual limitations
of serialized JLD2 data: SpeedyWeather.jl always has to be in the scope when loading
the data and the saved files might only load properly with exactly the same version
of SpeedyWeather.jl and Julia as used when saving the data.
Its usage is similar to the NetCDF output above:

```@example output2
using SpeedyWeather
spectral_grid = SpectralGrid(trunc=31, nlayers=1)
output = JLD2Output(output_dt=Hour(1))
model = ShallowWaterModel(spectral_grid, output=output)
model.output
```

With all options shown below 

```@example output2
@doc JLD2Output
```

## Array Output (to RAM)

It's also possible to output the `Variables` (or a subgroup of it) directly into an array that is kept in the memory. This might be useful e.g. for low-resolution simulations you want to work with and visualize quickly, or also for simulations within ML training loops. `ArrayOutput` follows the exact same logic as `JLD2Output` except that the actual output is kept in an array `output`. So, save e.g. all prognostic variables to memory, you can run: 

```@example output3
using SpeedyWeather
spectral_grid = SpectralGrid(trunc=31, nlayers=1)
output = ArrayOutput(output_dt=Hour(1), groups=(:prognostic,))
model = ShallowWaterModel(spectral_grid, output=output)
model.output
```

After a succesfull `run!` the result is stored in `output.output`. 

## Parameter summary

With `output=true` as an argument in the `run!(simulation)` call, the [NetCDFOutput](@ref) by default also
writes a parameter summary into `parameters.txt` in the same folder. This is implemented as
a [Callbacks](@ref) (`<: SpeedyWeather.AbstractCallback`) and can be added independent of
`NetCDFOutput` too. After `output=true` this callback is found here as `ParametersTxt`

```@example output2
# run an example simulation with output
simulation = initialize!(model)
run!(simulation, period=Hour(4), output=true)

# now model.output will have added additional "output" writers
simulation.model.callbacks
```

but it's activity is by default tied to activity of the `NetCDFOutput` with
you can control with `write_only_with_output`.
Creating such a callback independently

```@example output2
parameters_txt = ParametersTxt(write_only_with_output=false)
```

we can add it with a random or manual key as

```@example output2
add!(model, parameters_txt)             # random key
add!(model, :my_key => parameters_txt)  # manual key
```

But note that callbacks are overwritten with identical keys, otherwise treated independently.
Meaning now we have the preexisting `:parameters_txt` callback and then the two we just added,
writing their `parameters.txt` files one after another, overwriting that same file two times.

## Progress txt

Similarly to `ParametersTxt`, a callback `ProgressTxt` is by default added with `output=true`.
They can be created independently too

```@example output2
progress_txt = ProgressTxt(write_only_with_output=false, path="myfolder", filename="letsgo.txt")
```

and added like

```@example output2
add!(model, :progress_txt => progress_txt)
```

## Restart files for variables

`NetCDFOutput` also by default writes a restart file, containing the `simulation.variables.prognostic`
that can be read back in with the `StartFromFile` initial conditions. Implemented as a callback
`WriteVariablesRestartFile` can also be created independently of `NetCDFOutput`, e.g.

```@example output2
restart_file = WriteVariablesRestartFile(write_only_with_output=false, path="folder1", filename="restart.jld2")
```

and added like

```@example output2
add!(model, :restart_file => restart_file)
```

By default `path=""` will use the folder determined by `NetCDFOutput` but otherwise you can
also provide your own. Note that `WriteVariablesRestartFile` will only write the prognostic variables to file.
This is such that you can simulate a spin up and then change model parameters as you like,
to write out specific model components and store them in a file see
[Model component restart file](@ref).

## Model component restart file

If you modified a model component (say by applying a custom orography) you can save this to file too.

```@example output2
orography = ManualOrography(spectral_grid)  # an orography that is untouched at initialize!
model = ShallowWaterModel(spectral_grid; orography)
set!(model, orography=123)                  # set as you like
orography_for_restart = WriteModelComponentFile(component=model.orography, filename="my_orography.jld2")
add!(model, :my_orography => orography_for_restart)
```

Once the simulation ran you can then load this model component from file and use it to construct
a new model with it, or use its information in some other form, e.g. by writing its arrays to
other arrays:

```@example output2
simulation = initialize!(model)
run!(simulation, steps=0, output=true)

# either pass on that same callback (which will read out the path) or provide path directly
my_orography = SpeedyWeather.load_model_component(orography_for_restart)

# construct a new model with that orography
new_model = ShallowWaterModel(spectral_grid, orography = my_orography)
all(new_model.orography.orography .== 123)  # check that the new orography is indeed as customized
```

We do not really want to encourage it, but `WriteModelComponentFile` can be hijacked to write out the entire `model`.
While this easily saves everything of `model` into one file, it always writes many large precomputed arrays to file
whereas they could just be recomputed when constructing a new model. For example, the Legendre polynomials in
`model.spectral_transform` can easily be GBs at higher resolution, see also 
[Precomputed polynomials and allocated memory](@ref). Nevertheless, you can write out the entire model with

```@example output2
add!(model, :model_writer => WriteModelComponentFile(component=model, filename="model.jld2"))
```

Note that this callback contains a `model` that also contains this callback.
This self recursion is not particularly problematic as `model` is just a lazy reference.
However, when you do load in this `model` from file and use it again, note that it again
contains this callback which would write out its model again. You can `delete!` the callback
though, see [Adding a callback](@ref).