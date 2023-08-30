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
T127 (see [Available resolutions](@ref)) and `nlev=1` vertical levels.
The `SpectralGrid` object will provide us with some more information
```@example howtorun
using SpeedyWeather

spectral_grid = SpectralGrid(trunc=127,nlev=1)
```
We could have specified further options, but let's ignore that for now.
Next step we create a planet that's like Earth but not rotating
```@example howtorun
still_earth = Earth(rotation=0)
```
There are other options to create a planet but they are irrelevant for the
barotropic vorticity equations. We also want to specify the initial conditions,
randomly distributed vorticity is already defined
```@example howtorun
initial_conditions = StartWithRandomVorticity()
```
By default, the power of vorticity is spectrally distributed with ``k^{-3}``, ``k`` being the
horizontal wavenumber, and the amplitude is ``10^{-5}\text{s}^{-1}``.

Now we want to construct a `BarotropicModel`
with these
```@example howtorun
model = BarotropicModel(;spectral_grid, initial_conditions, planet=still_earth);
nothing # hide
```
The `model` contains all the parameters, but isn't initialized yet, which we can do
with and then run it. The `run!` command will always return the prognostic variables, which, by default, are 
plotted for surface relative vorticity with a unicode plot. The resolution of the plot
is not necessarily representative but it lets us have a quick look at the result
```@example howtorun
simulation = initialize!(model)
run!(simulation,n_days=30)
```

Woohoo! I can see turbulence! You could pick up where this simulation stopped by simply
doing `run!(simulation,n_days=50)` again. We didn't store any output, which
you can do by `run!(simulation,output=true)`, which will switch on NetCDF output
with default settings. More options on output in [NetCDF output](@ref).

## Example 2: Shallow water with mountains

As a second example, let's investigate the Galewsky et al.[^1] test case for the shallow
water equations with and without mountains. As the shallow water system has also only
one level, we can reuse the `SpectralGrid` from Example 1.
```@example howtorun
spectral_grid = SpectralGrid(trunc=127,nlev=1)
```
Now as a first simulation, we want to disable any orography, so we create a `NoOrography`
```@example howtorun
orography = NoOrography(spectral_grid)
```
Although the orography is zero, you have to pass on `spectral_grid` so that it can
still initialize zero-arrays of the right size and element type. Awesome.
This time the initial conditions should be set the the Galewsky et al.[^1] zonal
jet, which is already defined as
```@example howtorun
initial_conditions = ZonalJet()
```
The jet sits at 45˚N with a maximum velocity of 80m/s and a perturbation as described in their paper.
Now we construct a model, but this time a `ShallowWaterModel`
```@example howtorun
model = ShallowWaterModel(;spectral_grid, orography, initial_conditions)
simulation = initialize!(model);
run!(simulation,n_days=6)
```
Oh yeah. That looks like the wobbly jet in their paper. Let's run it again for another 6 days
but this time also store [NetCDF output](@ref).
```@example howtorun
run!(simulation,n_days=6,output=true)
```
The progress bar tells us that the simulation run got the identification "0001"
(which just counts up, so yours might be higher), meaning that
data is stored in the folder `/run_0001`, so let's plot that data properly (and not just using UnicodePlots).
```@example howtorun
using PyPlot, NCDatasets
ds = NCDataset("run_0001/output.nc");
ds["vor"]
```
Vorticity `vor` is stored as a 384x192x1x25 array, we may want to look at the first time step,
which is the end of the previous simulation (time=6days) which we didn't store output for.
```@example howtorun
vor = ds["vor"][:,:,1,1];
lat = ds["lat"][:];
lon = ds["lon"][:];
pcolormesh(lon,lat,vor')
nothing # hide
```
Which looks like

![Galewsky jet pyplot](https://raw.githubusercontent.com/SpeedyWeather/SpeedyWeather.jl/main/docs/img/galewsky_nc_6days.png)

You see that the unicode plot heavily coarse-grains the simulation, well it's unicode after all!
And now the last time step, that means time = 12days is
```@example howtorun
vor = ds["vor"][:,:,1,25];
pcolormesh(lon,lat,vor')
```

![Galewsky jet pyplot](https://raw.githubusercontent.com/SpeedyWeather/SpeedyWeather.jl/main/docs/img/galewsky_nc_12days.png)

The jet broke up into many small eddies, but the turbulence is still confined to the northern hemisphere, cool!
How this may change when we add mountains (we had `NoOrography` above!), say Earth's orography, you may ask?
Let's try it out! We create an `EarthOrography` struct like so

```@example howtorun
orography = EarthOrography(spectral_grid)
```

It will read the orography from file as shown, and there are some smoothing options too, but let's not change them.
Same as before, create a model, intialize into a simulation, run. This time directly for 12 days so that we can
compare with the last plot
```@example howtorun
model = ShallowWaterModel(;spectral_grid, orography, initial_conditions);
simulation = initialize!(model);
run!(simulation,n_days=12,output=true)
```
This time the run got a new run id, `0002` in our case, but otherwise we do as before.

![Galewsky jet pyplot](https://raw.githubusercontent.com/SpeedyWeather/SpeedyWeather.jl/main/docs/img/galewsky_nc_12days_mountains.png)

Interesting! The initial conditions have zero velocity in the southern hemisphere, but still, one can see
some imprint of the orography on vorticity. You can spot the coastline of Antarctica; the Andes and
Greenland are somewhat visible too. Mountains also completely changed the flow after 12 days,
probably not surprising!

## SpectralGrid

The life of every SpeedyWeather.jl simulation starts with a `SpectralGrid` object.
We have seen some examples above, now let's look into the details

```@example howtorun
?SpectralGrid
```

## References

[^1] Galewsky, Scott, Polvani, 2004. *An initial-value problem for testing numerical models of the global shallow-water equations*, Tellus A.
DOI: [10.3402/tellusa.v56i5.14436](https://doi.org/10.3402/tellusa.v56i5.14436)