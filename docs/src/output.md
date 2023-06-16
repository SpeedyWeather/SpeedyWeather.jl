# NetCDF output

SpeedyWeather.jl uses NetCDF to output the data of a simulation.
The following describes the details of this and how to change the way in which the NetCDF output is written.
There are many options to this available.

## Accessing the NetCDF output writer

The output writer is a component of every Model, i.e. `BarotropicModel`, `ShallowWaterModel`, `PrimitiveDryModel` and `PrimitiveWetModel`, hence a non-default output writer can be passed on as a keyword argument to the model constructor

```julia
julia> using SpeedyWeather
julia> spectral_grid = SpectralGrid()
julia> my_output_writer = OutputWriter(spectral_grid, PrimitiveDry)
julia> model = PrimitiveDryModel(;spectral_grid, output=my_output_writer)
```

So after we have defined the grid through the `SpectralGrid` object we can use and change
the implemented `OutputWriter` by passing on the following arguments
```julia
julia> my_output_writer = OutputWriter(spectral_grid, PrimitiveDry, kwargs...)
```
the `spectral_grid` has to be the first argument then the model type
(`Barotropic`, `ShallowWater`, `PrimitiveDry`, `PrimitiveWet`)
which helps the output writer to make default choices on which variables to output. However, we can
also pass on further keyword arguments. So let's start with an example.

## Example 1: NetCDF output every hour

If we want to increase the frequency of the output we can choose `output_dt` (default `=6` in hours) like so
```julia
julia> my_output_writer = OutputWriter(spectral_grid, PrimitiveDry, output_dt=1)
julia> model = PrimitiveDryModel(;spectral_grid, output=my_output_writer)
```
which will now output every hour. It is important to pass on the new output writer `my_output_writer` to the
model constructor, otherwise it will not be part of your model and the default is used instead.
Note that `output_dt` has to be understood as the minimum frequency or maximum output time step.
Example, we run the model at a resolution of T85 and the time step is going to be 670s
```julia
julia> spectral_grid = SpectralGrid(trunc=85)
julia> time_stepper = Leapfrog(spectral_grid)
Leapfrog{Float32}:
...
 Δt_sec::Int64 = 670
...
```
This means that after 32 time steps 5h 57min and 20s will have passed where output will happen as
the next time step would be >6h. The time axis of the NetCDF output will look like
```julia
julia> using NCDatasets
julia> ds = NCDataset("run_0001/output.nc");
julia> ds["time"][:]
5-element Vector{Dates.DateTime}:
 2000-01-01T00:00:00
 2000-01-01T05:57:20
 2000-01-01T11:54:40
 2000-01-01T17:52:00
 2000-01-01T23:49:20

julia> diff(ds["time"][:])
4-element Vector{Dates.Millisecond}:
 21440000 milliseconds
 21440000 milliseconds
 21440000 milliseconds
 21440000 milliseconds
```
This is so that we don't interpolate in time during output to hit exactly every 6 hours, but at
the same time have a constant spacing in time between output time steps.

## Example 2: Output onto a higher/lower resolution grid

Say we want to run the model at a given horizontal resolution but want to output on another resolution,
the `OutputWriter` takes as argument `output_Grid<:AbstractFullGrid` and `nlat_half::Int`.
So for example `output_Grid=FullClenshawGrid` and `nlat_half=48` will always interpolate onto a
regular 192x95 longitude-latitude grid of 1.875˚ resolution, regardless the grid and resolution used
for the model integration.
```julia
julia> my_output_writer = OutputWriter(spectral_grid, PrimitiveDry, output_Grid=FullClenshawGrid, nlat_half=48)
```
Note that by default the output is on the corresponding full of the grid used in the dynamical core
so that interpolation only happens at most in the zonal direction as they share the location of the
latitude rings. You can check this by
```julia
julia> RingGrids.full_grid(OctahedralGaussianGrid)
FullGaussianGrid
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
```
julia> my_output_writer = OutputWriter(spectral_grid,PrimitiveDry,id="diffusion_test");
```
which will be used instead of a 4 digit number like 0001, 0002 which is automatically determined if
`id` is not provided. You will see the id of the run in the progress bar
```
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

```@docs
OutputWriter
```