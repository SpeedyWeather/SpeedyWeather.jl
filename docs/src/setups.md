# Model setups

The following is a collection of model setups, starting with an easy setup
of the [Barotropic vorticity equation](@ref) and continuing with more
complicated setups.

## 2D turbulence on a non-rotating sphere

!!! info "Setup script"
    ```julia
    using SpeedyWeather
    spectral_grid = SpectralGrid(trunc=63,nlev=1)
    still_earth = Earth(rotation=0)
    initial_conditions = StartWithRandomVorticity()
    model = BarotropicModel(;spectral_grid, initial_conditions, planet=still_earth)
    simulation = initialize!(model)
    run!(simulation,n_days=20)
    ```

We want to use the barotropic model to simulate some free-decaying 2D turbulence
on the sphere without rotation. We start by defining the `SpectralGrid` object.
To have a resolution of about 200km, we choose a spectral resolution of
T63 (see [Available horizontal resolutions](@ref)) and `nlev=1` vertical levels.
The `SpectralGrid` object will provide us with some more information
```@example barotropic_setup
using SpeedyWeather
spectral_grid = SpectralGrid(trunc=63,nlev=1)
```
Next step we create a planet that's like Earth but not rotating
```@example barotropic_setup
still_earth = Earth(rotation=0)
```
There are other options to create a planet but they are irrelevant for the
barotropic vorticity equations. We also want to specify the initial conditions,
randomly distributed vorticity is already defined
```@example barotropic_setup
using Random # hide
Random.seed!(1234) # hide
initial_conditions = StartWithRandomVorticity()
```
By default, the power of vorticity is spectrally distributed with ``k^{-3}``, ``k`` being the
horizontal wavenumber, and the amplitude is ``10^{-5}\text{s}^{-1}``.

Now we want to construct a `BarotropicModel`
with these
```@example barotropic_setup
model = BarotropicModel(;spectral_grid, initial_conditions, planet=still_earth)
nothing # hide
```
The `model` contains all the parameters, but isn't initialized yet, which we can do
with and then run it. The `run!` command will always return the prognostic variables, which, by default, are 
plotted for surface relative vorticity with a unicode plot. The resolution of the plot
is not necessarily representative but it lets us have a quick look at the result
```@example barotropic_setup
simulation = initialize!(model)
run!(simulation,n_days=20)
```

Woohoo! Something is moving! You could pick up where this simulation stopped by simply
doing `run!(simulation,n_days=50)` again. We didn't store any output, which
you can do by `run!(simulation,output=true)`, which will switch on NetCDF output
with default settings. More options on output in [NetCDF output](@ref).

## Shallow water with mountains

!!! info "Setup script"
    ```julia
    using SpeedyWeather
    spectral_grid = SpectralGrid(trunc=63,nlev=1)
    orography = NoOrography(spectral_grid)
    initial_conditions = ZonalJet()
    model = ShallowWaterModel(;spectral_grid, orography, initial_conditions)
    simulation = initialize!(model)
    run!(simulation,n_days=6)
    ```

As a second example, let's investigate the Galewsky et al.[^G04] test case for the shallow
water equations with and without mountains. As the shallow water system has also only
one level, we can reuse the `SpectralGrid` from Example 1.
```@example galewsky_setup
using SpeedyWeather
spectral_grid = SpectralGrid(trunc=63,nlev=1)
```
Now as a first simulation, we want to disable any orography, so we create a `NoOrography`
```@example galewsky_setup
orography = NoOrography(spectral_grid)
```
Although the orography is zero, you have to pass on `spectral_grid` so that it can
still initialize zero-arrays of the correct size and element type. Awesome.
This time the initial conditions should be set the the Galewsky et al.[^G04] zonal
jet, which is already defined as
```@example galewsky_setup
initial_conditions = ZonalJet()
```
The jet sits at 45˚N with a maximum velocity of 80m/s and a perturbation as described in their paper.
Now we construct a model, but this time a `ShallowWaterModel`
```@example galewsky_setup
model = ShallowWaterModel(;spectral_grid, orography, initial_conditions)
simulation = initialize!(model)
run!(simulation,n_days=6)
```
Oh yeah. That looks like the wobbly jet in their paper. Let's run it again for another 6 days
but this time also store [NetCDF output](@ref).
```@example galewsky_setup
run!(simulation,n_days=6,output=true)
```
The progress bar tells us that the simulation run got the identification "0001"
(which just counts up, so yours might be higher), meaning that
data is stored in the folder `/run_0001`, so let's plot that data properly (and not just using UnicodePlots).
```@example galewsky_setup
using PythonPlot, NCDatasets
ioff() # hide
ds = NCDataset("run_0001/output.nc")
ds["vor"]
```
Vorticity `vor` is stored as a lon x lat x vert x time array, we may want to look at the first time step,
which is the end of the previous simulation (time=6days) which we didn't store output for.
```@example galewsky_setup
t = 1
vor = Matrix{Float32}(ds["vor"][:,:,1,t]) # convert from Matrix{Union{Missing,Float32}} to Matrix{Float32}
lat = ds["lat"][:]
lon = ds["lon"][:]

fig,ax = subplots(1,1,figsize=(10,6))
ax.pcolormesh(lon,lat,vor')
ax.set_xlabel("longitude")
ax.set_ylabel("latitude")
ax.set_title("Relative vorticity")
tight_layout() # hide
savefig("galewsky1.png") # hide
nothing # hide
```
![Galewsky jet pyplot1](galewsky1.png)

You see that in comparison the unicode plot heavily coarse-grains the simulation, well it's unicode after all!
And now the last time step, that means time = 12days is

```@example galewsky_setup
t = ds.dim["time"]
vor = Matrix{Float32}(ds["vor"][:,:,1,t])
ax.pcolormesh(lon,lat,vor')
savefig("galewsky2.png") # hide
nothing # hide
```
![Galewsky jet pyplot2](galewsky2.png)

The jet broke up into many small eddies, but the turbulence is still confined to the northern hemisphere, cool!
How this may change when we add mountains (we had `NoOrography` above!), say Earth's orography, you may ask?
Let's try it out! We create an `EarthOrography` struct like so

```@example galewsky_setup
orography = EarthOrography(spectral_grid)
```

It will read the orography from file as shown, and there are some smoothing options too, but let's not change them.
Same as before, create a model, initialize into a simulation, run. This time directly for 12 days so that we can
compare with the last plot

```@example galewsky_setup
model = ShallowWaterModel(;spectral_grid, orography, initial_conditions)
simulation = initialize!(model)
run!(simulation,n_days=12,output=true)
```

This time the run got a new run id, which you see in the progress bar, but can also always check
after the `run!` call (the automatic run id is only determined just before the main time loop starts)
with `model.output.id`, but otherwise we do as before.
```@example galewsky_setup
id = model.output.id
```

```@example galewsky_setup
ds = NCDataset("run_$id/output.nc")
time = 49
vor = Matrix{Float32}(ds["vor"][:,:,1,time])

fig,ax = subplots(1,1,figsize=(10,6))
ax.pcolormesh(lon,lat,vor')
ax.set_xlabel("longitude")
ax.set_ylabel("latitude")
ax.set_title("Relative vorticity")
tight_layout() # hide
savefig("galewsky3.png") # hide
nothing # hide
```
![Galewsky jet pyplot3](galewsky3.png)

Interesting! The initial conditions have zero velocity in the southern hemisphere, but still, one can see
some imprint of the orography on vorticity. You can spot the coastline of Antarctica; the Andes and
Greenland are somewhat visible too. Mountains also completely changed the flow after 12 days,
probably not surprising!

## Polar jet streams in shallow water

Setup script:
```@example jet_stream_setup
using SpeedyWeather
spectral_grid = SpectralGrid(trunc=63,nlev=1)
forcing = JetStreamForcing(spectral_grid,latitude=60)
drag = QuadraticDrag(spectral_grid)
output = OutputWriter(spectral_grid,ShallowWater,output_dt=Hour(6),output_vars=[:u,:v,:pres,:orography])
model = ShallowWaterModel(;spectral_grid,output,drag,forcing)
simulation = initialize!(model)
model.feedback.verbose = false # hide
run!(simulation,n_days=20)   # discard first 20 days   
run!(simulation,n_days=20,output=true)
nothing # hide
```

We want to simulate polar jet streams in the shallow water model. We add a `JetStreamForcing`
that adds momentum at 60˚N to inject kinetic energy into the model. This energy needs to be removed
(the [diffusion](@ref diffusion) is likely not sufficient) through a drag, we have implemented
a `QuadraticDrag` and use the default drag coefficient. Outputting ``u,v,\eta`` (called `:pres`,
as it is the pressure equivalent in the shallow water system) we run 20 days without output
to give the system some time to adapt to the forcing. And visualize zonal wind after another
20 days with

```@example jet_stream_setup
using PythonPlot, NCDatasets
ioff() # hide

id = model.output.id
ds = NCDataset("run_$id/output.nc")
timestep = ds.dim["time"]
u = Matrix{Float32}(ds["u"][:,:,1,timestep])
lat = ds["lat"][:]
lon = ds["lon"][:]

fig,ax = subplots(1,1,figsize=(10,6))
q = ax.pcolormesh(lon,lat,u')
ax.set_xlabel("longitude")
ax.set_ylabel("latitude")
ax.set_title("Zonal wind [m/s]")
colorbar(q,ax=ax)
tight_layout() # hide
savefig("polar_jets.png") # hide
nothing # hide
```
![Polar jets pyplot](polar_jets.png)


## Gravity waves on the sphere

Setup script:
```@example gravity_wave_setup
using Random # hide
Random.seed!(1234) # hide
using SpeedyWeather
spectral_grid = SpectralGrid(trunc=127,nlev=1)
time_stepping = SpeedyWeather.Leapfrog(spectral_grid,Δt_at_T31=30)
implicit = SpeedyWeather.ImplicitShallowWater(spectral_grid,α=0.5)
orography = EarthOrography(spectral_grid,smoothing=false)
initial_conditions = SpeedyWeather.RandomWaves()
output = OutputWriter(spectral_grid,ShallowWater,output_dt=Hour(12),output_vars=[:u,:pres,:div,:orography])
model = ShallowWaterModel(;spectral_grid,orography,output,initial_conditions,implicit,time_stepping)
simulation = initialize!(model)
model.feedback.verbose = false # hide
run!(simulation,n_days=2,output=true)
nothing # hide
```

How are gravity waves propagating around the globe? We want to use the shallow water model
to start with some random perturbations of the interface displacement (the "sea surface height")
but zero velocity and let them propagate around the globe. We set the ``\alpha`` parameter
of the [semi-implicit time integration](@ref implicit_swm) to ``0.5`` to have a centred
implicit scheme which dampens the gravity waves less than a backward implicit scheme would do.
But we also want to keep orography, and particularly no smoothing on it, to have the orography
as rough as possible. The initial conditions are set to `RandomWaves` which set the spherical
harmonic coefficients of ``\eta`` to between given wavenumbers to some random values
```@example gravity_wave_setup
SpeedyWeather.RandomWaves()
```
so that the amplitude `A` is as desired, here 2000m. Our layer thickness is by default
```@example gravity_wave_setup
model.atmosphere.layer_thickness
```
8.5km so those waves are with an amplitude of 2000m quite strong.
But the semi-implicit time integration can handle that even with fairly large time steps of
```@example gravity_wave_setup
model.time_stepping.Δt_sec
```
seconds. Note that the gravity wave speed here is ``\sqrt{gH}`` so almost 300m/s.
Let us also output divergence, as gravity waves are quite pronounced in that variable.
But given the speed of gravity waves we don't have to integrate for long.
Visualise with


```@example gravity_wave_setup
using PythonPlot, NCDatasets
ioff() # hide

id = model.output.id
ds = NCDataset("run_$id/output.nc")
timestep = ds.dim["time"]
div = Matrix{Float32}(ds["div"][:,:,1,timestep])
lat = ds["lat"][:]
lon = ds["lon"][:]

fig,ax = subplots(1,1,figsize=(10,6))
ax.pcolormesh(lon,lat,div')
ax.set_xlabel("longitude")
ax.set_ylabel("latitude")
ax.set_title("Divergence")
tight_layout() # hide
savefig("gravity_waves.png") # hide
nothing # hide
```
![Gravity waves pyplot](gravity_waves.png)

Can you spot the Himalayas or the Andes?

## Jablonowski-Williamson baroclinic wave

```@example jablonowski
using SpeedyWeather
spectral_grid = SpectralGrid(trunc=31,nlev=8,Grid=FullGaussianGrid,dealiasing=3)
orography = ZonalRidge(spectral_grid)
initial_conditions = ZonalWind()
model = PrimitiveDryModel(;spectral_grid,orography,initial_conditions,physics=false)
simulation = initialize!(model)
model.feedback.verbose = false # hide
run!(simulation,n_days=9,output=true)
nothing # hide
```

The Jablonowski-Williamson baroclinic wave test case[^JW06] using the [Primitive equation model](@ref)
particularly the dry model, as we switch off all physics with `physics=false`.
We want to use 8 vertical levels, and a lower resolution of T31 on a [full Gaussian grid](@ref FullGaussianGrid).
The Jablonowski-Williamson initial conditions are in `ZonalWind`, the orography
is just a `ZonalRidge`. There is no forcing and the initial conditions are
baroclinically unstable which kicks off a wave propagating eastward.
This wave becomes obvious when visualised with

```@example jablonowski
using PythonPlot, NCDatasets
ioff() # hide

id = model.output.id
ds = NCDataset("run_$id/output.nc")
timestep = ds.dim["time"]
surface = ds.dim["lev"]
vor = Matrix{Float32}(ds["vor"][:,:,surface,timestep])
lat = ds["lat"][:]
lon = ds["lon"][:]

fig,ax = subplots(1,1,figsize=(10,6))
ax.pcolormesh(lon,lat,vor')
ax.set_xlabel("longitude")
ax.set_ylabel("latitude")
ax.set_title("Surface relative vorticity")
tight_layout() # hide
savefig("jablonowski.png") # hide
nothing # hide
```
![Jablonowski pyplot](jablonowski.png)


## References

[^G04]: Galewsky, Scott, Polvani, 2004. *An initial-value problem for testing numerical models of the global shallow-water equations*, Tellus A. DOI: [10.3402/tellusa.v56i5.14436](https://doi.org/10.3402/tellusa.v56i5.14436)
[^JW06]: Jablonowski, C. and Williamson, D.L. (2006), A baroclinic instability test case for atmospheric model dynamical cores. Q.J.R. Meteorol. Soc., 132: 2943-2975. [10.1256/qj.06.12](https://doi.org/10.1256/qj.06.12)