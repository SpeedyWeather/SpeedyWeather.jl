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
To have a resolution of about 200km, we choose a spectral resolution of
T63 (see [Available horizontal resolutions](@ref)) and `nlev=1` vertical levels.
The `SpectralGrid` object will provide us with some more information
```@example howtorun
using SpeedyWeather
spectral_grid = SpectralGrid(trunc=63,nlev=1)
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
run!(simulation,n_days=20)
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
spectral_grid = SpectralGrid(trunc=85,nlev=1)
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
ioff() # hide
ds = NCDataset("run_0001/output.nc");
ds["vor"]
```
Vorticity `vor` is stored as a lon x lat x vert x time array, we may want to look at the first time step,
which is the end of the previous simulation (time=6days) which we didn't store output for.
```@example howtorun
vor = ds["vor"][:,:,1,1];
lat = ds["lat"][:];
lon = ds["lon"][:];
pcolormesh(lon,lat,vor')
xlabel("longitude")
ylabel("latitude")
tight_layout() # hide
nothing # hide
savefig("galewsky1.png") # hide
```
Which looks like

![Galewsky jet pyplot](galewsky1.png)

You see that the unicode plot heavily coarse-grains the simulation, well it's unicode after all!
And now the last time step, that means time = 12days is
```@example howtorun
vor = ds["vor"][:,:,1,25];
pcolormesh(lon,lat,vor')
xlabel("longitude")
ylabel("latitude")
tight_layout() # hide
nothing # hide
savefig("galewsky2.png") # hide
```

![Galewsky jet pyplot](galewsky2.png)

The jet broke up into many small eddies, but the turbulence is still confined to the northern hemisphere, cool!
How this may change when we add mountains (we had `NoOrography` above!), say Earth's orography, you may ask?
Let's try it out! We create an `EarthOrography` struct like so

```@example howtorun
orography = EarthOrography(spectral_grid)
```

It will read the orography from file as shown, and there are some smoothing options too, but let's not change them.
Same as before, create a model, initialize into a simulation, run. This time directly for 12 days so that we can
compare with the last plot
```@example howtorun
model = ShallowWaterModel(;spectral_grid, orography, initial_conditions)
simulation = initialize!(model)
run!(simulation,n_days=12,output=true)
```
This time the run got a new run id, `0002` in our case, but otherwise we do as before.

```@example howtorun
dsvor = NCDataset("run_0002/output.nc")["vor"]
```

```@example howtorun
vor = dsvor[:,:,1,49]
pcolormesh(lon,lat,vor')
xlabel("longitude")
ylabel("latitude")
tight_layout() # hide
nothing # hide
savefig("galewsky3.png") # hide
```
Which looks like

![Galewsky jet pyplot](galewsky3.png)

Interesting! The initial conditions have zero velocity in the southern hemisphere, but still, one can see
some imprint of the orography on vorticity. You can spot the coastline of Antarctica; the Andes and
Greenland are somewhat visible too. Mountains also completely changed the flow after 12 days,
probably not surprising!

## SpectralGrid

The life of every SpeedyWeather.jl simulation starts with a `SpectralGrid` object.
We have seen some examples above, now let's look into the details

```julia
help?> SpectralGrid
search: SpectralGrid

  Defines the horizontal spectral resolution and corresponding grid and the
  vertical coordinate for SpeedyWeather.jl. Options are

    •  NF::Type{<:AbstractFloat}: number format used throughout the model

    •  trunc::Int64: horizontal resolution as the maximum degree of
       spherical harmonics

    •  Grid::Type{<:SpeedyWeather.RingGrids.AbstractGrid}: horizontal
       grid used for calculations in grid-point space

    •  dealiasing::Float64: how to match spectral with grid resolution:
       dealiasing factor, 1=linear, 2=quadratic, 3=cubic grid

    •  radius::Float64: radius of the sphere [m]

    •  nlat_half::Int64: number of latitude rings on one hemisphere
       (Equator incl)

    •  npoints::Int64: total number of grid points in the horizontal

    •  nlev::Int64: number of vertical levels

    •  vertical_coordinates::SpeedyWeather.VerticalCoordinates:
       coordinates used to discretize the vertical

  nlat_half and npoints should not be chosen but are derived from trunc, Grid
  and dealiasing.
```

Say we wanted double precision (`Float64`), a spectral resolution of T42 on
a regular longitude-latitude grid (`FullClenshawGrid`) with cubic truncation
(`dealiasing=3`) and 4 vertical levels, we would do this by

```@example howtorun
spectral_grid = SpectralGrid(NF=Float64, trunc=42, Grid=FullClenshawGrid, dealiasing=3, nlev=4)
```

We don't specify the resolution of the grid (its `nlat_half` parameter) directly,
instead we chose a spectral truncation `trunc`
and through the `dealiasing` factor a grid resolution will be automatically chosen. 
Here T42 will be a matched with a 192x95 regular longitude-latitude grid
that has 18240 grid points in total. For details see [Matching spectral and grid resolution](@ref).

We could have also defined `SpectralGrid` on a smaller sphere than Earth,
or with a different vertical spacing
```@example howtorun
vertical_coordinates = SigmaCoordinates(0:0.2:1)
```
These are regularly spaced [Sigma coordinates](@ref), defined through their half levels.
```@example howtorun
spectral_grid = SpectralGrid(;vertical_coordinates,radius=1e6)
```

In the end, a `SpectralGrid` defines the physical domain of the simulation and its discretization.
This domain has to be a sphere because of the spherical harmonics, but it can have a different radius.
The discretization is for spectral, grid-point space and the vertical as this determines the size of many
arrays for preallocation, for which als the number format is essential to know.
That's why `SpectralGrid` is the beginning of every SpeedyWeather.jl simulation and that is why
it has to be passed on to (most) model components.

## References

[^1] Galewsky, Scott, Polvani, 2004. *An initial-value problem for testing numerical models of the global shallow-water equations*, Tellus A.
DOI: [10.3402/tellusa.v56i5.14436](https://doi.org/10.3402/tellusa.v56i5.14436)