# How to run SpeedyWeather.jl

Creating a SpeedyWeather.jl simulation and running it consists conceptually of 4 steps

1. Create a `SpectralGrid` which defines the grid and spectral resolution
2. Create a model
3. Initialize a model to obtain a `Simulation`.
4. Run the simulation.

In the following we will describe these steps in more detail,
but let's start with some examples first.

## Example 1: 2D turbulence on a non-rotating sphere

We want to use the barotropic model to simulate some free-decaying 2D turbulence
on the sphere without rotation. We start by defining the `SpectralGrid` object.
To have a resolution of about 100km, we choose a spectral resolution of
T127 (see [Available horizontal resolutions](@ref)) and `nlev=1` vertical levels.
The `SpectralGrid` object will provide us with some more information
```julia
julia> spectral_grid = SpectralGrid(trunc=127,nlev=1)
SpectralGrid:
 Spectral:   T127 LowerTriangularMatrix{Complex{Float32}}, radius = 6.371e6 m
 Grid:       40320-element, 192-ring OctahedralGaussianGrid{Float32} (quadratic)
 Resolution: 112km (average)
 Vertical:   1-level SigmaCoordinates
```
We could have specified further options, but let's ignore that for now.
Next step we create a planet that's like Earth but not rotating
```julia
julia> still_earth = Earth(rotation=0)
Main.SpeedyWeather.Earth
  rotation: Float64 0.0
  gravity: Float64 9.81
  daily_cycle: Bool true
  length_of_day: Float64 24.0
  seasonal_cycle: Bool true
  length_of_year: Float64 365.25
  equinox: Dates.DateTime
  axial_tilt: Float64 23.4
```
There are other options to create a planet but they are irreleveant for the
barotropic vorticity equations. We also want to specify the initial conditions,
randomly distributed vorticity is already defined
```julia
julia> initial_conditions = StartWithRandomVorticity()
StartWithRandomVorticity
  power_law: Float64 -3.0
  amplitude: Float64 1.0e-5
```
By default, the power of vorticity is spectrally distributed with ``k^{-3}``, ``k`` being the
horizontal wavenumber, and the amplitude is ``10^{-5}\text{ s}^{-1}``.

Now we want to construct a `BarotropicModel`
with these
```julia
julia> model = BarotropicModel(;spectral_grid, initial_conditions, planet=still_earth);
```
The `model` contains all the parameters, but isn't initialized yet, which we can do
with and then run it.
```julia
julia> simulation = initialize!(model);
julia> run!(simulation,n_days=30)
```
The `run!` command will always return the prognostic variables, which, by default, are 
plotted for surface relative vorticity with a unicode plot. The resolution of the plot
is not necessarily representative but it lets us have a quick look at the result

![Barotropic vorticity unicode plot](https://raw.githubusercontent.com/SpeedyWeather/SpeedyWeather.jl/main/docs/img/barotropic_vorticity.jpg)

Woohoo! I can see turbulence! You could pick up where this simulation stopped by simply
doing `run!(simulation,n_days=50)` again. We didn't store any output, which
you can do by `run!(simulation,output=true)`, which will switch on NetCDF output
with default settings. More options on output in [NetCDF output](@ref).

## Example 2: Shallow water with mountains

As a second example, let's investigate the Galewsky et al.[^1] test case for the shallow
water equations with and without mountains. As the shallow water system has also only
one level, we can reuse the `SpectralGrid` from Example 1.
```julia
julia> spectral_grid = SpectralGrid(trunc=127,nlev=1)
```
Now as a first simulation, we want to disable any orography, so we create a `NoOrography`
```julia
julia> orography = NoOrography(spectral_grid)
NoOrography{Float32, OctahedralGaussianGrid{Float32}}
```
Although the orography is zero, you have to pass on `spectral_grid` so that it can
still initialize zero-arrays of the right size and element type. Awesome.
This time the initial conditions should be set the the Galewsky et al.[^1] zonal
jet, which is already defined as
```julia
julia> initial_conditions = ZonalJet()
ZonalJet
  latitude: Float64 45.0
  width: Float64 19.28571428571429
  umax: Float64 80.0
  perturb_lat: Float64 45.0
  perturb_lon: Float64 270.0
  perturb_xwidth: Float64 19.098593171027442
  perturb_ywidth: Float64 3.819718634205488
  perturb_height: Float64 120.0
```
The jet sits at 45˚N with a maximum velocity of 80m/s and a perturbation as described in their paper.
Now we construct a model, but this time a `ShallowWaterModel`
```
julia> model = ShallowWaterModel(;spectral_grid, orography, initial_conditions);
julia> simulation = initialize!(model);
```
![Galewsky jet unicode plot](https://raw.githubusercontent.com/SpeedyWeather/SpeedyWeather.jl/main/docs/img/galewsky.jpg)

Oh yeah. That looks like the wobbly jet in their paper. Let's run it again for another 6 days
but this time also store [NetCDF output](@ref).
```
julia> run!(simulation,n_days=6,output=true)
Weather is speedy: run 0002 100%|███████████████████████| Time: 0:00:12 (115.37 years/day)
```
The progress bar tells us that the simulation run got the identification "0002", meaning that
data is stored in the folder `/run_0002`, so let's plot that data properly (and not just using UnicodePlots).
```julia
julia> using PyPlot, NCDatasets
julia> ds = NCDataset("run_0002/output.nc");
julia> ds["vor"]
vor (384 × 192 × 1 × 25)
  Datatype:    Float32
  Dimensions:  lon × lat × lev × time
  Attributes:
   units                = 1/s
   missing_value        = NaN
   long_name            = relative vorticity
   _FillValue           = NaN
```
Vorticity `vor` is stored as a 384x192x1x25 array, we may want to look at the first time step,
which is the end of the previous simulation (time=6days) which we didn't store output for.
```julia
julia> vor = ds["vor"][:,:,1,1];
julia> lat = ds["lat"][:];
julia> lon = ds["lon"][:];
julia> pcolormesh(lon,lat,vor')
```
Which looks like

![Galewsky jet pyplot](https://raw.githubusercontent.com/SpeedyWeather/SpeedyWeather.jl/main/docs/img/galewsky_nc_6days.png)

You see that the unicode plot heavily coarse-grains the simulation, well it's unicode after all!
And now the last time step, that means time=12days is
```julia
julia> vor = ds["vor"][:,:,1,25];
julia> pcolormesh(lon,lat,vor')
```

![Galewsky jet pyplot](https://raw.githubusercontent.com/SpeedyWeather/SpeedyWeather.jl/main/docs/img/galewsky_nc_12days.png)

The jet broke up into many small eddies, but the turbulence is still confined to the northern hemisphere, cool!
How this may change when we add mountains (we had `NoOrography` above!), say Earth's orography, you may ask?
Let's try it out! We create an `EarthOrography` struct like so

```julia
julia> orography = EarthOrography(spectral_grid)
EarthOrography{Float32, OctahedralGaussianGrid{Float32}}:
 path::String = SpeedyWeather.jl/input_data
 file::String = orography_F512.nc
 scale::Float64 = 1.0
 smoothing::Bool = true
 smoothing_power::Float64 = 1.0
 smoothing_strength::Float64 = 0.1
 smoothing_truncation::Int64 = 85
```

It will read the orography from file as shown, and there are some smoothing options too, but let's not change them.
Same as before, create a model, intialize into a simulation, run. This time directly for 12 days so that we can
compare with the last plot
```julia
julia> model = ShallowWaterModel(;spectral_grid, orography, initial_conditions);
julia> simulation = initialize!(model);
julia> run!(simulation,n_days=12,output=true)
Weather is speedy: run 0003 100%|███████████████████████| Time: 0:00:35 (79.16 years/day)
```
This time the run got the id "0003", but otherwise we do as before.

![Galewsky jet pyplot](https://raw.githubusercontent.com/SpeedyWeather/SpeedyWeather.jl/main/docs/img/galewsky_nc_12days_mountains.png)

Interesting! The initial conditions have zero velocity in the southern hemisphere, but still, one can see
some imprint of the orography on vorticity. You can spot the coastline of Antarctica; the Andes and
Greenland are somewhat visible too. Mountains also completely changed the flow after 12 days,
probably not surprising!

## SpectralGrid

The life of every SpeedyWeather.jl simulation starts with a `SpectralGrid` object.
We have seen some examples above, now let's look into the details

```@docs
SpectralGrid
```

## References

[^1] Galewsky, Scott, Polvani, 2004. *An initial-value problem for testing numerical models of the global shallow-water equations*, Tellus A.
DOI: [10.3402/tellusa.v56i5.14436](https://doi.org/10.3402/tellusa.v56i5.14436)