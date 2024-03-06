# NetCDF output

SpeedyWeather.jl uses NetCDF to output the data of a simulation.
The following describes the details of this and how to change the way in which the NetCDF output is written.
There are many options to this available.

## Accessing the NetCDF output writer

The output writer is a component of every Model, i.e. `BarotropicModel`, `ShallowWaterModel`, `PrimitiveDryModel` and `PrimitiveWetModel`, hence a non-default output writer can be passed on as a keyword argument to the model constructor

```@example netcdf
using SpeedyWeather
spectral_grid = SpectralGrid()
output = OutputWriter(spectral_grid, ShallowWater)
model = ShallowWaterModel(; spectral_grid, output=output)
nothing # hide
```

So after we have defined the grid through the `SpectralGrid` object we can use and change
the implemented `OutputWriter` by passing on additional arguments.
The `spectral_grid` has to be the first argument then the model type
(`Barotropic`, `ShallowWater`, `PrimitiveDry`, or `PrimitiveWet`)
which helps the output writer to make default choices on which variables to output.
Then we can also pass on further keyword arguments. So let's start with an example.

## Example 1: NetCDF output every hour

If we want to increase the frequency of the output we can choose `output_dt` (default `=Hour(6)`) like so
```@example netcdf
output = OutputWriter(spectral_grid, ShallowWater, output_dt=Hour(1))
model = ShallowWaterModel(; spectral_grid, output=output)
nothing # hide
```
which will now output every hour. It is important to pass on the new output writer `output` to the
model constructor, otherwise it will not be part of your model and the default is used instead.
Note that the choice of `output_dt` can affect the actual time step that is used for the model
integration, which is explained in the following.
Example, we run the model at a resolution of T42 and the time step is going to be
```@example netcdf
spectral_grid = SpectralGrid(trunc=42, nlev=1)
time_stepping = Leapfrog(spectral_grid)
time_stepping.Δt_sec
```
seconds. Depending on the output frequency (we chose `output_dt = Hour(1)` above)
this will be slightly adjusted during model initialization:
```@example netcdf
output = OutputWriter(spectral_grid, ShallowWater, output_dt=Hour(1))
model = ShallowWaterModel(; spectral_grid, time_stepping, output)
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
model = ShallowWaterModel(; spectral_grid, time_stepping, output)
simulation = initialize!(model)
```

The time axis of the NetCDF output will now look like
```@example netcdf
using NCDatasets
model.feedback.verbose = false # hide
run!(simulation, period=Day(1), output=true)
id = model.output.id
ds = NCDataset("run_$id/output.nc")
ds["time"][:]
```
which is a bit ugly, that's why `adjust_with_output=true` is the default. In that case we would have
```@example netcdf
time_stepping = Leapfrog(spectral_grid, adjust_with_output=true)
output = OutputWriter(spectral_grid, ShallowWater, output_dt=Hour(1))
model = ShallowWaterModel(; spectral_grid, time_stepping, output)
simulation = initialize!(model)
run!(simulation, period=Day(1), output=true)
id = model.output.id
ds = NCDataset("run_$id/output.nc")
ds["time"][:]
```
very neatly hourly output in the NetCDF file!

## Example 2: Output onto a higher/lower resolution grid

Say we want to run the model at a given horizontal resolution but want to output on another resolution,
the `OutputWriter` takes as argument `output_Grid<:AbstractFullGrid` and `nlat_half::Int`.
So for example `output_Grid=FullClenshawGrid` and `nlat_half=48` will always interpolate onto a
regular 192x95 longitude-latitude grid of 1.875˚ resolution, regardless the grid and resolution used
for the model integration.
```julia
my_output_writer = OutputWriter(spectral_grid, ShallowWater, output_Grid=FullClenshawGrid, nlat_half=48)
```
Note that by default the output is on the corresponding full of the grid used in the dynamical core
so that interpolation only happens at most in the zonal direction as they share the location of the
latitude rings. You can check this by
```@example netcdf
RingGrids.full_grid(OctahedralGaussianGrid)
```
So the corresponding full grid of an `OctahedralGaussianGrid` is the `FullGaussiangrid` and the same resolution
`nlat_half` is chosen by default in the output writer (which you can change though as shown above).
Overview of the corresponding full grids

| Grid | Corresponding full grid |
| ---  | ----------------------- |
| FullGaussianGrid | FullGaussianGrid |
| FullClenshawGrid | FullClenshawGrid |
| OctahadralGaussianGrid | FullGaussianGrid |
| OctahedralClensawhGrid | FullClenshawGrid |
| HEALPixGrid | FullHEALPixGrid |
| OctaHEALPixGrid | FullOctaHEALPixGrid |

The grids `FullHEALPixGrid`, `FullOctaHEALPixGrid` share the same latitude rings as their reduced grids,
but have always as many longitude points as they are at most around the equator. These grids are not
tested in the dynamical core (but you may use them experimentally) and mostly designed for output purposes.

## Example 3: Changing the output path or identification

That's easy by passing on `path="/my/favourite/path/"` and the folder `run_*` with `*` the identification
of the run (that's the `id` keyword, which can be manually set but is also automatically determined as a
number counting up depending on which folders already exist) will be created within.
```julia
julia> path = pwd()
"/Users/milan"
julia> my_output_writer = OutputWriter(spectral_grid, PrimitiveDry, path=path)
```
This folder must already exist. If you want to give your run a name/identification you can pass on `id`
```julia
julia> my_output_writer = OutputWriter(spectral_grid, PrimitiveDry, id="diffusion_test");
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

Further options are described in the `OutputWriter` docstring, (also accessible via `julia>?OutputWriter` for example).
Note that some fields are actual options, but others are derived from the options you provided or are
arrays/objects the output writer needs, but shouldn't be passed on by the user.
The actual options are declared as `[OPTION]` in the following

```@example netcdf
@doc OutputWriter
```
