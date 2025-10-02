# Other output

Next to the [NetCDF output](@ref) there's other ways to output data from
the model, these are described below.

## JLD2 Output 

As an alternative to the [NetCDF output](@ref), it is also possible to directly
output the `PrognosticVariables` and `DiagnosticVariables` to a JLD2 file.
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

## Restart file

`NetCDFOutput` also by default writes a restart file, containing the `simulation.prognostic_variables`
that can be read back in with teh `StartFromFile` initial conditions. Implemented as a callback
`RestartFile` can also be created independently of `NetCDFOutput`, e.g.

```@example output2
restart_file = RestartFile(write_only_with_output=false, path="folder1", filename="restart.jld2")
```

and added like

```@example output2
add!(model, :restart_file => restart_file)
```

By default `path=""` will use the folder determined by `NetCDFOutput` but otherwise you can
also provide your own.