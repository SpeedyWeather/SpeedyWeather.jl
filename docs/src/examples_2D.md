# [Examples 2D](@id Examples)

The following is a collection of example model setups, starting with an easy setup
of the [Barotropic vorticity equation](@ref barotropic_vorticity_model) and
continuing with the [shallow water model](@ref shallow_water_model).

See also [Examples 3D](@ref) for examples with the primitive equation models.

## 2D turbulence on a non-rotating sphere

!!! info "Setup script to copy and paste"
    ```julia
    using SpeedyWeather
    spectral_grid = SpectralGrid(trunc=63, nlayers=1)
    still_earth = Earth(spectral_grid, rotation=0)
    initial_conditions = RandomVelocity()
    forcing = NoForcing()
    drag = NoDrag()
    model = BarotropicModel(spectral_grid; initial_conditions, planet=still_earth, forcing, drag)
    simulation = initialize!(model)
    run!(simulation, period=Day(20))
    ```

We want to use the barotropic model to simulate some free-decaying 2D turbulence
on the sphere without rotation. We start by defining the `SpectralGrid` object.
To have a resolution of about 200km, we choose a spectral resolution of
T63 (see [Available horizontal resolutions](@ref)) and `nlayers=1` vertical levels.
The `SpectralGrid` object will provide us with some more information
```@example barotropic_setup
using SpeedyWeather
spectral_grid = SpectralGrid(trunc=63, nlayers=1)
```
Next step we create a planet that's like Earth but not rotating. As a convention,
we always pass on the spectral grid object as the first argument to every other
model component we create.
```@example barotropic_setup
still_earth = Earth(spectral_grid, rotation=0)
```
There are other options to create a planet but they are irrelevant for the
barotropic vorticity equations. We also want to specify the initial conditions,
randomly distributed velocity is already defined
```@example barotropic_setup
initial_conditions = RandomVelocity()
```
By default, the velocity has an approximate amplitude as given and also higher
wavenumbers are truncated.

For free-decaying turbulence we switch off any forcing or drag with

```@example barotropic_setup
forcing = NoForcing()
drag = NoDrag()
```

Now we want to construct a `BarotropicModel` with these simply by passing
them as keyword arguments either with `key=argument` or just `; key` which matches
the argument name with the keyword
```@example barotropic_setup
model = BarotropicModel(spectral_grid; initial_conditions, planet=still_earth, forcing, drag)
nothing # hide
```
The `model` contains all the parameters, but isn't initialized yet, which we can do
with and then run it.

```@example barotropic_setup
simulation = initialize!(model)
run!(simulation, period=Day(20))
```

Woohoo! You could pick up where this simulation stopped by simply
doing `run!(simulation, period=Day(50))` again. We didn't store any output, which
you can do by `run!(simulation, output=true)`, which will switch on NetCDF output
with default settings. More options on output in [NetCDF output](@ref).

## Shallow water with mountains

!!! info "Setup script to copy and past"
    ```julia
    using SpeedyWeather
    spectral_grid = SpectralGrid(trunc=63, nlayers=1)
    orography = NoOrography(spectral_grid)
    initial_conditions = ZonalJet()
    model = ShallowWaterModel(spectral_grid; orography, initial_conditions)
    simulation = initialize!(model)
    run!(simulation, period=Day(6))
    ```

As a second example, let's investigate the Galewsky et al.[^G04] test case for the shallow
water equations with and without mountains. As the shallow water system has also only
one level, we can reuse the `SpectralGrid` from Example 1.
```@example galewsky_setup
using SpeedyWeather
spectral_grid = SpectralGrid(trunc=63, nlayers=1)
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
model = ShallowWaterModel(spectral_grid; orography, initial_conditions)
simulation = initialize!(model)
run!(simulation, period=Day(6))
```
Oh yeah. That looks like the wobbly jet in their paper. Let's run it again for another 6 days
but this time also store [NetCDF output](@ref).
```@example galewsky_setup
run!(simulation, period=Day(6), output=true)
```
The progress bar tells us that the simulation run got the identification "run_0001"
(which just counts up, so yours might be higher), meaning that
data is stored in the folder `run_0001`. In general we can check this also via
```@example galewsky_setup
model.output.run_folder
```

### Visualisation (manually)

So let's plot that data. `$id` in the following just means
that the string is interpolated to `run_0001` if this is the first unnamed run in your folder.
```@example galewsky_setup
using NCDatasets
ds = NCDataset("run_$id/output.nc")
ds["vor"]
```
Vorticity `vor` is stored as a lon x lat x vert x time array, we may want to look at the first time step,
which is the end of the previous simulation (time = 6days) which we didn't store output for.
```@example galewsky_setup
t = 1
vor = Matrix{Float32}(ds["vor"][:, :, 1, t]) # convert from Matrix{Union{Missing, Float32}} to Matrix{Float32}
lat = ds["lat"][:]
lon = ds["lon"][:]

using CairoMakie
heatmap(lon, lat, vor)
save("galewsky0.png", ans) # hide
nothing # hide
```
![Galewsky jet](galewsky0.png)

Here, we have unpacked the netCDF file using [NCDatasets.jl](https://github.com/Alexander-Barth/NCDatasets.jl)
and then plotted via `heatmap(lon, lat, vor)`. While you can do that to give you more control
on the plotting, SpeedyWeather.jl also defines an extension for Makie.jl, see [Extensions](@ref).
Because if our matrix `vor` here was an `AbstractGrid` (see [RingGrids](@ref)) then all
its geographic information (which grid point is where) would be implicitly known from the type.
From the netCDF file, however, you would need to use the longitude and latitude dimensions.

So we can also just do (`input_as=Matrix` here as all our grids use and expect a horizontal dimension
flattened into a vector by default)

```@example galewsky_setup
vor_grid = FullGaussianGrid(vor, input_as=Matrix)

using CairoMakie    # this will load the extension so that Makie can plot grids directly
heatmap(vor_grid, title="Relative vorticity [1/s]")
save("galewsky1.png", ans) # hide
nothing # hide
```
![Galewsky jet plot1](galewsky1.png)

Note that here you need to know which grid the data comes on (an error is thrown if `FullGaussianGrid(vor)`
is not size compatible). By default the output will be on the FullGaussianGrid, but if you
play around with other grids, you'd need to change this here,
see [NetCDF output](@ref) and [Output grid](@ref).

### Visualisation via Makie

We did want to showcase the usage of [NetCDF output](@ref) here, but from now on
we will use `heatmap` to plot data on our grids directly, without storing output first.
So for our current simulation, that means at time = 12 days, vorticity on the grid
is stored in the diagnostic variables and can be visualised with
(`[:, 1]` is horizontal x vertical dimension, so all grid points on the first and
only vertical layer)

```@example galewsky_setup
vor = simulation.diagnostic_variables.grid.vor_grid[:, 1]
heatmap(vor, title="Relative vorticity [1/s]")
save("galewsky2.png", ans) # hide
nothing # hide
```
![Galewsky jet](galewsky2.png)

The jet broke up into many small eddies, but the turbulence is still confined to the northern hemisphere, cool!

### Visualisation via UnicodePlots

Similar to the Makie extension that is loaded automatically with `using CairoMakie`
(or another backend like `using GLMakie`) we have defined an extension for
[UnicodePlots.jl](https://github.com/JuliaPlots/UnicodePlots.jl)

```@example galewsky_setup
using UnicodePlots
UnicodePlots.heatmap(vor)       # UnicodePlots. only needed as CairoMakie was loaded first!
```

Note that if you first load Makie and then UnicodePlots then calling `heatmap`
without further specifying the scope like `UnicodePlots.heatmap` will continue
to use the `heatmap` from Makie. Add `UnicodePlots.` or `CairoMakie.` if
needed.

UnicodePlots look good in the Julia REPL, especially with the more colours
a terminal like [iTerm](https://iterm2.com/) provides. But these plots are
literally made up of unicode characters so the most you can expect is somthing like

<img src="https://github.com/SpeedyWeather/SpeedyWeather.jl/assets/25530332/a04fbb10-1cc1-4f77-93f2-7bdf047f277d" width="450"><br>

However, here in the documentation they are usually vertically spaced as the line
spacing is by default higher than in the REPL.
A similar issue arises in Jupyter Notebooks by default. Well, they are unicode
after all!

### Adding mountains

How this may change when we add mountains (we had `NoOrography` above!), say Earth's orography, you may ask?
Let's try it out! We create an `EarthOrography` struct like so

```@example galewsky_setup2
using SpeedyWeather # hide
spectral_grid = SpectralGrid(trunc=63, nlayers=1) # hide
initial_conditions = ZonalJet() # hide                    
orography = EarthOrography(spectral_grid)
```

It will read the orography from file as shown (only at `initialize!(model)`), and there are some smoothing
options too, but let's not change them. Same as before, create a model, initialize into a simulation, run.
This time directly for 12 days so that we can compare with the last plot

```@example galewsky_setup2
model = ShallowWaterModel(spectral_grid; orography, initial_conditions)
simulation = initialize!(model)
run!(simulation, period=Day(12), output=true)
```

This time the run got a new `run_number`, which you see in the progress bar, but can again always check
after the `run!` call (the automatic `run_number` is only determined just before the main time loop starts)
with `model.output.run_folder`, but otherwise we do as before.
```@example galewsky_setup2
run_folder = model.output.run_folder
```

```@example galewsky_setup2
using NCDatasets
ds = NCDataset("$run_folder/output.nc")
```

While you could plot the [NetCDF output](@ref) manually as before, 
we'll be plotting directly from the current state of the `simulation` using
the Makie extension

```@example galewsky_setup2
using CairoMakie
vor = simulation.diagnostic_variables.grid.vor_grid[:, 1]   # 1 to index surface
heatmap(vor, title="Relative vorticity [1/s]")
save("galewsky3.png", ans) # hide
nothing # hide
```
![Galewsky jet](galewsky3.png)

Interesting! One can clearly see some imprint of the orography on vorticity and there is especially
more vorticity in the southern hemisphere. You can spot the coastline of Antarctica; the Andes and
Greenland are somewhat visible too. Mountains also completely changed the flow after 12 days,
probably not surprising!

## Polar jet streams in shallow water

Setup script to copy and paste:
```@example jet_stream_setup
using SpeedyWeather
spectral_grid = SpectralGrid(trunc=63, nlayers=1)

forcing = JetStreamForcing(spectral_grid, latitude=60)
drag = LinearVorticityDrag(spectral_grid)

model = ShallowWaterModel(spectral_grid; drag, forcing)
simulation = initialize!(model)
run!(simulation, period=Day(40))
nothing # hide
```

We want to simulate polar jet streams in the shallow water model. We add a `JetStreamForcing`
that adds momentum at 60˚N and 60˚S an to inject kinetic energy into the model. This energy needs to be removed
(the [diffusion](@ref diffusion) is likely not sufficient) through a drag, we have implemented
a `LinearVorticityDrag` and use the default drag coefficient. Then visualize zonal wind after
40 days with

```@example jet_stream_setup
using CairoMakie

u = simulation.diagnostic_variables.grid.u_grid[:, 1]
heatmap(u, title="Zonal wind [m/s]")
save("polar_jets.png", ans) # hide
nothing # hide
```
![Polar jets plot](polar_jets.png)


## Gravity waves on the sphere

Setup script to copy and paste:
```@example gravity_wave_setup
using Random # hide
Random.seed!(1234) # hide
using SpeedyWeather
spectral_grid = SpectralGrid(trunc=127, nlayers=1)

# model components
implicit = ImplicitShallowWater(spectral_grid, α=0.5)
orography = EarthOrography(spectral_grid, smoothing=false)
initial_conditions = RandomWaves(lmin=10, lmax=30)      # between wavenumber 10 and 30

# construct, initialize, run
model = ShallowWaterModel(spectral_grid; orography, initial_conditions, implicit)
simulation = initialize!(model)
run!(simulation, period=Day(2))
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
RandomWaves()
```
so that the amplitude `A` is as desired, here 2000m. Our layer thickness in meters is by default
```@example gravity_wave_setup
model.atmosphere.layer_thickness
```
so those waves are with an amplitude of 2000m quite strong.
But the semi-implicit time integration can handle that even with fairly large time steps of
```@example gravity_wave_setup
model.time_stepping.Δt_sec
```
seconds. Note that the gravity wave speed here is ``\sqrt{gH}`` so almost 300m/s,
given the speed of gravity waves we don't have to integrate for long.
Visualise the dynamic layer thickness ``h = \eta + H + H_b`` 
(interface displacement ``\eta``, layer thickness at rest ``H`` and orography ``H_b``)
with

```@example gravity_wave_setup
using CairoMakie

H = model.atmosphere.layer_thickness
Hb = model.orography.orography
η = simulation.diagnostic_variables.grid.pres_grid
h = @. η + H - Hb   # @. to broadcast grid + scalar - grid

heatmap(h, title="Dynamic layer thickness h", colormap=:oslo)
save("gravity_waves.png", ans) # hide
nothing # hide
```
![Gravity waves plot](gravity_waves.png)

Mountains like the Himalayas or the Andes are quite obvious because the atmospheric layer
is much thinner there. The pressure gradient is relative to ``z=0`` so in a fluid
at rest the mountains would just "reach into" the fluid, thinning the layer the higher
the mountain. As the atmosphere here is not at rest the layer thickness is not perfectly
(anti-)correlated with orography but almost so.

## References

[^G04]: Galewsky, Scott, Polvani, 2004. *An initial-value problem for testing numerical models of the global shallow-water equations*, Tellus A. DOI: [10.3402/tellusa.v56i5.14436](https://doi.org/10.3402/tellusa.v56i5.14436)
