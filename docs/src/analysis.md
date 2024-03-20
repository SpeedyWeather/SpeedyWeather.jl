# Analysing a simulation

While you can analyze a SpeedyWeather simulation through its [NetCDF output](@ref),
as most users will be used to with other models, you can also reuse a lot of
functionality from SpeedyWeather interactively for analysis. This makes
SpeedyWeather beyond being a model for the atmospheric general circulation
also a library with many functions for the analysis of simulations.
Often this also avoids the two language problem that you will face if you
run a simulation with a model in one language but then do the data
analysis in another, treating the model as a blackbox although it likely
has many of the functions you will need for analysis already defined.
With SpeedyWeather we try to avoid this and are working towards a 
more unified approach in atmospheric modelling where simulation
and analysis are done interactively with the same library: SpeedyWeather.jl.

## Advantages of online analysis

Now you could run a SpeedyWeather simulation, and analyse the [NetCDF output](@ref)
but that comes with several issues related to accuracy

- If you use a reduced grid for the simulation, then the output will (by default) be
interpolated on a full grid. This interpolation comes introduces an error.

- Computing integrals over gridded data on the sphere by weighting every grid point
according to its area is not the most accurate numerical integration.

- Computing gradients over gridded data comes with similar issues. While our
[RingGrids](@ref) are always equidistant in longitude, they are not necessarily
in latitude.

The first point you can avoid by running a simulation on one of the full grids that
are implemented, see [SpectralGrid](@ref). But that also impacts the simulation
and for various reasons we don't run on full grids by default.
The second point you can address by defining a more advanced numerical integration
scheme, but that likely requires you to depend on external libraries and then,
well, you could also just depend on SpeedyWeather.jl directly, because we
have to do these computations internally anyway. Similar for the third point,
many gradients have to be computed on every time step and we do that with
spectral transforms to reduce the discretization error.

The following contains a (hopefully growing) list of examples
of how a simulation can be analysed interactively. We call this
_online_ analysis because you are directly using the functionality from
the model as if it was a library.

## Mass conservation

In the absence of sources and sinks for the interface displacement ``\eta``
in the [shallow water equations](@ref shallow_water_model), total mass
(or equivalently volume as density is constant) is conserved.
The total volume is defined as the integral of the dynamic layer thickness
``h = \eta + H - H_b`` (``H`` is the layer thickness at rest, ``H_b`` is orography)
over the surface ``A`` of the sphere

```math
\iint h dA = \iint \eta dA + \iint H dA - \iint H_b dA
```

to check for conservation we want to assess that

```math
\frac{\partial}{\partial t} \iint h dA = 0
```

And because ``V = \iint H dA - \iint H_b dA``, the total volume at rest,
is a constant (``H`` is a global constant, the orography ``H_b`` does not change with time)
we can just check whether ``\iint \eta dA`` changes over time.
Instead of computing this integral in grid-point space, we use the spectral
``\eta`` whereby the coefficient of the first spherical harmonic (the ``l = m = 0`` mode,
or wavenumber 0, see [Spherical Harmonic Transform](@ref)) encodes the global average.

```@example analysis
using SpeedyWeather
spectral_grid = SpectralGrid(trunc=31, nlev=1)
model = ShallowWaterModel(;spectral_grid)
simulation = initialize!(model)
```

Now we check ``\eta_{0,0}`` the ``l = m = 0`` coefficent of the inital conditions
of that simulation with

```@example analysis
simulation.prognostic_variables.surface.timesteps[1].pres[1]
```

`[1]` pulls the first element of the underlying [LowerTriangularMatrix](@ref lowertriangularmatrices)
which is the coefficient of the ``l = m = 0`` mode.
Its imaginary part is always zero (which is true for any zonal harmonic ``m=0`` as its
imaginary part would just unnecessarily rotate something zonally constant in zonal direction),
so you can `real` it. Also for spherical harmonic transforms there is a norm of the sphere
by which you have to divide to get your mean value in the original units

```@example analysis
a = model.spectral_transform.norm_sphere    # = 2√π = 3.5449078
η_mean = real(simulation.prognostic_variables.surface.timesteps[1].pres[1]) / a
```

So the initial conditions in this simulation are such that the global mean interface displacement
is that value in meters. You would need to multiply by the area of the sphere
``4\pi r^2`` (radius ``r``) to get the actual integral from above, but because that doesn't
change with time either, we just want to check that `η_mean` doesn't change with time.
Which is equivalent to ``\partial_t \iint \eta dA = 0`` and so volume conservation and because density is constant
also mass conservation. Let's check what happens after the simulation ran for some days

```@example analysis
model.feedback.verbose = false # hide
run!(simulation, period=Day(10))

# now we check η_mean again
η_mean_later = real(simulation.prognostic_variables.surface.timesteps[1].pres[1]) / a
```

which is _exactly_ the same. So mass is conserved, woohoo. 

Insight from a numerical perspective: The tendency of ``\eta`` is
``\partial_t \eta = -\nabla \cdot (\mathbf{u} h)`` which is a divergence of a flux.
Calculating the divergence in spherical harmonics always sets the ``l=m=0`` mode to zero exactly
(the gradient of a periodic variable has zero mean) so the time integration here is always with
an exactly zero tendency.

## Energy conservation

The total energy in the shallow water equation is the sum of the kinetic energy
``\frac{1}{2}(u^2 + v^2)`` and the potential energy ``gh`` integrated over the
total volume (times ``h`` for the vertical then integrated over the sphere``\iint dA``).

```math
\iint \frac{h}{2}(u^2 + v^2) + gh^2 dA
```

In contrast to the [Mass conservation](@ref) which, with respect to the
spectral transform, is a linear calculation, here we need to multiply
variables, which is done in grid-point space. Then we can transform to spectral
space for the global integral. Let us define a `total_energy` function as

```@example analysis
using SpeedyWeather
function total_energy(u, v, η, model)
    
    h = zero(u)
    E = zero(u)                             # allocate grid variable
    
    H = model.atmosphere.layer_thickness
    Hb = model.orography.orography
    g = model.planet.gravity
    
    @. h = η + H - Hb
    @. E = h/2*(u^2 + v^2) + g*h^2

    # transform to spectral, take l=m=0 mode at [1] and normalize for mean
    E_mean = real(spectral(E)[1]) / model.spectral_transform.norm_sphere
end
```

So at the current state of our simulation we have a total energy
(per square meter as we haven't multiplied by the surface area)

```@example analysis
# flat copies for convenience
u = simulation.diagnostic_variables.layers[1].grid_variables.u_grid
v = simulation.diagnostic_variables.layers[1].grid_variables.v_grid
η = simulation.diagnostic_variables.surface.pres_grid

TE = total_energy(u, v, η, model)
```

with units of ``m^3 s^{-2}`` (multiplying by surface area of the sphere
and density of the fluid would turn it into joule = ``kg m^2 s^{-2}``).
Now let us continue the simulation
```@example analysis
run!(simulation, period=Day(10))

# we don't need to reassign u, v, η as they were flat copies
# pointing directly to the diagnostic variables inside the simulation
# which got updated during run!
TE_later = total_energy(u, v, η, model)
```
So the total energy has somewhat changed, it decreased to
```@example analysis
TE_later/TE
```
of its previous value over 10 days. While technically energy should
be conserved in an unforced system, numerically this is rarely
exactly the case. We need some [Horizontal diffusion](@ref diffusion)
for numerical stability and also the time integration is dissipative
due to temporal filtering, see [Time integration](@ref leapfrog).

## Potential vorticity

Potential vorticity in the shallow water equations is defined as

```math
q = \frac{f + \zeta}{h}
```

with ``f`` the Coriolis parameter, relative vorticity ``\zeta``
and ``h`` the layer thickness as before. We can calculate this
conveniently directly on the model grid (whichever you chose)
as

```@example analysis
# vorticity
ζ = simulation.diagnostic_variables.layers[1].grid_variables.vor_grid
f = coriolis(ζ)     # create f on that grid

# layer thickness
η = simulation.diagnostic_variables.surface.pres_grid
h = zero(η)
H = model.atmosphere.layer_thickness
Hb = model.orography.orography
@. h = η + H - Hb

# potential vorticity
q = zero(ζ)
@. q = (f + ζ) / h
nothing # hide
```

and we can compare the relative vorticity field to
```@example analysis
plot(ζ)
```
the potential vorticity
```@example analysis
plot(q)
```

## Absolute angular momentum

More to follow ...

## Circulation

More to follow ...

## Enstrophy

More to follow ...