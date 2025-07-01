# NetCDF output

SpeedyWeather.jl uses NetCDF to output the data of a simulation.
The following describes the details of this and how to change the way in which the NetCDF output is written.
There are many options to this available.

## Creating `NetCDFOutput`

```@example netcdf
using SpeedyWeather
spectral_grid = SpectralGrid()
output = NetCDFOutput(spectral_grid)
```

With `NetCDFOutput(::SpectralGrid, ...)` one creates a `NetCDFOutput` writer with several options,
which are explained in the following. By default, the `NetCDFOutput` is created when constructing
the model, i.e.

```@example netcdf
model = ShallowWaterModel(spectral_grid)
model.output
```

The output writer is a component of every Model, i.e. `BarotropicModel`, `ShallowWaterModel`, `PrimitiveDryModel`
and `PrimitiveWetModel`, and they only differ in their default `output.variables` (e.g. the primitive
models would by default output temperature which does not exist in the 2D models `BarotropicModel` or `ShallowWaterModel`).
But any `NetCDFOutput` can be passed onto the model constructor with the `output` keyword argument.

```@example netcdf
output = NetCDFOutput(spectral_grid, Barotropic)
model = ShallowWaterModel(spectral_grid, output=output)
nothing # hide
```

Here, we created `NetCDFOutput` for the model class `Barotropic` (2nd positional argument, outputting only vorticity and velocity)
but use it in the `ShallowWaterModel`. By default the `NetCDFOutput` is set to inactive, i.e.
`output.active` is `false`. It is only turned on (and initialized) with `run!(simulation, output=true)`.
So you may change the `NetCDFOutput` as you like but only calling `run!(simulation)` will not
trigger it as `output=false` is the default here.

## Output frequency

If we want to increase the frequency of the output we can choose `output_dt` (default `=Hour(6)`) like so
```@example netcdf
output = NetCDFOutput(spectral_grid, ShallowWater, output_dt=Hour(1))
model = ShallowWaterModel(spectral_grid, output=output)
model.output
```
which will now output every hour. It is important to pass on the new output writer `output` to the
model constructor, otherwise it will not be part of your model and the default is used instead.
Note that the choice of `output_dt` can affect the actual time step that is used for the model
integration, which is explained in the following.
Example, we run the model at a resolution of T42 and the time step is going to be
```@example netcdf
spectral_grid = SpectralGrid(trunc=42, nlayers=1)
time_stepping = Leapfrog(spectral_grid)
time_stepping.Δt_sec
```
seconds. Depending on the output frequency (we chose `output_dt = Hour(1)` above)
this will be slightly adjusted during model initialization:
```@example netcdf
output = NetCDFOutput(spectral_grid, ShallowWater, output_dt=Hour(1))
model = ShallowWaterModel(spectral_grid; time_stepping, output)
simulation = initialize!(model)
model.time_stepping.Δt_sec
```
The shorter the output time step the more the model time step needs to be adjusted
to match the desired output time step exactly. This is important so that for daily output at
noon this does not slowly shift towards night over years of model integration.
One can always disable this adjustment with
```@example netcdf
time_stepping = Leapfrog(spectral_grid, adjust_with_output=false)
time_stepping.Δt_sec
```
and a little info will be printed to explain that even though you wanted
`output_dt = Hour(1)` you will not actually get this upon initialization:
```@example netcdf
model = ShallowWaterModel(spectral_grid; time_stepping, output)
simulation = initialize!(model)
```

The time axis of the NetCDF output will now look like
```@example netcdf
using NCDatasets
run!(simulation, period=Day(1), output=true)
id = model.output.id
ds = NCDataset("run_$id/output.nc")
ds["time"][:]
```
which is a bit ugly, that's why `adjust_with_output=true` is the default. In that case we would have
```@example netcdf
time_stepping = Leapfrog(spectral_grid, adjust_with_output=true)
output = NetCDFOutput(spectral_grid, ShallowWater, output_dt=Hour(1))
model = ShallowWaterModel(spectral_grid; time_stepping, output)
simulation = initialize!(model)
run!(simulation, period=Day(1), output=true)
id = model.output.id
ds = NCDataset("run_$id/output.nc")
ds["time"][:]
```
very neatly hourly output in the NetCDF file!

## Output grid

Say we want to run the model at a given horizontal resolution but want to output on another resolution,
the `NetCDFOutput` takes as argument `output_grid::AbstractFullGrid`, any instance of a full grid
can be provided here.
So for example `output_grid=FullClenshawGrid(48)` would interpolate onto a
regular 192x95 longitude-latitude grid of 1.875˚ resolution, regardless the grid and resolution used
for the model integration.
```@example netcdf
my_output_writer = NetCDFOutput(spectral_grid, output_grid=FullClenshawGrid(48))
```
Note that by default the output is on the corresponding full type of the grid type used in the dynamical core
so that interpolation only happens at most in the zonal direction as they share the location of the
latitude rings. You can check this by
```@example netcdf
RingGrids.full_grid_type(OctahedralGaussianGrid)
```
So the corresponding full grid of an `OctahedralGaussianGrid` is the `FullGaussianGrid` and the same resolution
`nlat_half` is chosen by default in the output writer (which you can change though as shown above).
Overview of the corresponding full grids

| Grid | Corresponding full grid |
| ---  | ----------------------- |
| FullGaussianGrid | FullGaussianGrid |
| FullClenshawGrid | FullClenshawGrid |
| OctahadralGaussianGrid | FullGaussianGrid |
| OctaminimalGaussianGrid | FullGaussianGrid |
| OctahedralClenshawGrid | FullClenshawGrid |
| HEALPixGrid | FullHEALPixGrid |
| OctaHEALPixGrid | FullOctaHEALPixGrid |

The grids `FullHEALPixGrid`, `FullOctaHEALPixGrid` share the same latitude rings as their reduced grids,
but have always as many longitude points as there are around the equator. These grids are not
tested in the dynamical core (but you may use them experimentally) and mostly designed for output purposes.

## Output variables

One can easily add or remove variables from being output with the `NetCDFOut` writer. The following
variables are predefined (note they are not exported so you have to prefix `SpeedyWeather.`)

```@example netcdf
using InteractiveUtils # hide
subtypes(SpeedyWeather.AbstractOutputVariable)
```

"Defined" here means that every such type contains information about a variables (long) name,
its units, dimensions, any missing values and compression options. For `HumidityOutput` for example
we have

```@example netcdf
SpeedyWeather.HumidityOutput()
```

You can choose name and unit as you like, e.g. `SpeedyWeather.HumidityOutput(unit = "1")` or change
the compression options, e.g. `SpeedyWeather.HumidityOutput(keepbits = 5)` but more customisation
is discussed in [Customizing netCDF output](@ref).

We can add new output variables with `add!` 

```@example netcdf
output = NetCDFOutput(spectral_grid)            # default variables
add!(output, SpeedyWeather.DivergenceOutput())  # output also divergence
output
```

If you didn't create a `NetCDFOutput` separately, you can also apply this directly to `model`,
either `add!(model, SpeedyWeather.DivergenceOutput())` or `add!(model.output, args...)`,
which technically also just forwards to `add!(model.output.variables, args...)`.
`output.variables` is a dictionary were the variable names (as `Symbol`s) are used as keys,
so `output.variables[:div]` just returns the `SpeedyWeather.DivergenceOutput()` we have
just created using `:div` as key. With those keys one can also `delete!` a variable
from netCDF output

```@example netcdf
delete!(output, :div)
```

If you change the `name` of an output variable, i.e. `SpeedyWeather.DivergenceOutput(name="divergence")`
the key would change accordingly to `:divergence`.

## Grouped variables

For convenience we have defined several output groups, for example `SpeedyWeather.PrecipitationOutput()`
definces both accumulated large-scale and convective precipitation as well as their rates and the
cloud top height. Groups are

```julia
SpeedyWeather.BoundaryOutput()
SpeedyWeather.PrecipitationOutput()
SpeedyWeather.RadiationOutput()
SpeedyWeather.SurfacesFluxesOutput()
SpeedyWeather.AllOutputVariabesl()
```

each of them has to be splatted by appending `...`, e.g.

```@example netcdf
add!(model, SpeedyWeather.SurfaceFluxesOutput()...)
```

## Output path and identification

That's easy by passing on `path="/my/favourite/path/"` and the folder `run_*` with `*` the identification
of the run (that's the `id` keyword, which can be manually set but is also automatically determined as a
number counting up depending on which folders already exist) will be created within.
```julia
julia> path = pwd()
"/Users/milan"
julia> my_output_writer = NetCDFOutput(spectral_grid, path=path)
```
This folder must already exist. If you want to give your run a name/identification you can pass on `id`
```julia
julia> my_output_writer = NetCDFOutput(spectral_grid, id="diffusion_test");
```
which will be used instead of a 4 digit number like 0001, 0002 which is automatically determined if
`id` is not provided. You will see the id of the run in the progress bar
```julia
Weather is speedy: run diffusion_test 100%|███████████████████████| Time: 0:00:12 (19.20 years/day)
```
and the run folder, here `run_diffusion_test`, is also named accordingly
```julia
shell> ls
...
run_diffusion_test
...
```

## Further options

Further options are described in the `NetCDFOutput` docstring, (also accessible via `julia>?NetCDFOutput` for example).
Note that some fields are actual options, but others are derived from the options you provided or are
arrays/objects the output writer needs, but shouldn't be passed on by the user.
The actual options are declared as `[OPTION]` in the following

```@example netcdf
@doc NetCDFOutput
```

## Visualizing output

The saved NetCDF files can be visualized with a wide range of tools, both in Julia, but also in other languages. In order to get a quick view into a NetCDF file, you can use command line tools like `ncview`. For actual visualizations in Julia, it's easy to use [NCDatasets.jl](https://github.com/JuliaGeo/NCDatasets.jl) for accessing the data and [GeoMakie.jl](https://github.com/JuliaGeo/GeoMakie.jl) for plotting it. For a standard animation we already provide the `animate` function within SpeedyWeather.jl's GeoMakie extension that makes it easy to animate a variable from a NetCDF output file or a `Simulation` object, as seen below:

```@example netcdf
using GeoMakie
@doc SpeedyWeather.animate
```

# JLD2 Output 

As an alternative to the NetCDF output, it is also possible to directly output the `PrognosticVariables` and `DiagnosticVariables` to a JLD2 file. This might be interesting if you are really interested in the model internals, or also for some machine learning tasks. However, this option doesn't feature any possibilites to regrid or select variables, and it comes with the usual limitations of serialized JLD2 data: SpeedyWeather.jl always has to be in the scope when loading the data and the saved files might only load properly with exactly the same version of SpeedyWeather.jl and Julia as used when saving the data. Its usage is similar to the NetCDF output above:

```@example netcdf 
spectral_grid = SpectralGrid()
output = JLD2Output(output_dt=Hour(1))
model = ShallowWaterModel(spectral_grid, output=output)
model.output
```

With all options shown below 

```@example netcdf 
@doc JLD2Output
```