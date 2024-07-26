# Orography

Orography (in height above the surface) forms the surface boundary of the
lowermost layer in SpeedyWeather. 

In the [shallow-water equations](@ref shallow_water_model) the orography
``H_b`` enters the equations when computing the layer thickness ``h = \eta + H_0 - H_b``
for the volume fluxes ``\mathbf{u}h`` in the continuity equation.
Here, the orography is used in meters above the surface which shortens
``h`` over mountains. The orography here is needed in grid-point space.

In the [primitive equations](@ref primitive_equation_model) the orography enters
the equations when computing the [Geopotential](@ref). So actually required here
is the surface geopotential ``\Phi_s = gz_s`` where ``z_s``
is the orography height in meters as used in the shallow-water equations too
``z_s = H_b``.
However, the primitive equations require the orography in spectral
space as the geopotential calculation is a linear operation in the
horizontal and can therefore be applied in either grid-point or spectral space.
The latter is more convenient as SpeedyWeather solves the equations
to avoid additional transforms.

In the current formulation of the
[barotropic vorticity equations](@ref barotropic_vorticity_model) there is no
orography. In fact, the field `model.orography` is not defined for
`model::BarotropicModel`.

## Orographies implemented

Currently implemented are
```@example orography
using InteractiveUtils # hide
using SpeedyWeather
subtypes(SpeedyWeather.AbstractOrography)
```
which are 

- ``\Phi_s = z_s = H_b = 0`` for `NoOrography`
- For `ZonalRidge` the zonal ridge from the Jablonowski and Williamson initial conditions, see [Jablonowski-Williamson baroclinic wave](@ref)
- For `EarthOrography` a high-resolution orography is loaded and interpolated to the resolution as defined by `spectral_grid`.

all orographies need to be created with `spectral_grid::SpectralGrid` as the first argument,
so that the respective fields for `geopot_surf`, i.e. ``\Phi_s`` and `orography`, i.e. ``H_b``
can be allocated in the right size and number format.

## Earth's orography

Earth's orography can be created with (here we use a resolution of T85, about 165km globally)

```@example orography
using SpeedyWeather
spectral_grid = SpectralGrid(trunc=85)
orography = EarthOrography(spectral_grid)
```

but note that only allocates the orography, it does not actually load
and interpolate the orography which happens at the `initialize!`
step. Visualised with

```@example orography
model = PrimitiveDryModel(;spectral_grid, orography)
initialize!(orography, model)   # happens also in simulation = initialize!(model)

using CairoMakie
heatmap(orography.orography, title="Earth's orography at T85 resolution, no smoothing")
save("earth_orography.png", ans) # hide
nothing # hide
```
![EarthOrography](earth_orography.png)

typing `?EarthOrography` shows the various options that are provided.
An orogaphy at T85 resolution that is as smooth as it would be at T42
(controlled by the `smoothing_fraction`, the fraction of highest wavenumbers
which are the top half here, about T43 to T85) for example can be created with

```@example orography
orography = EarthOrography(spectral_grid, smoothing=true, smoothing_fraction=0.5)
initialize!(orography, model)

heatmap(orography.orography, title="Earth's orography at T85 resolution, smoothed to T42")
save("earth_orography_smooth.png", ans) # hide
nothing # hide
```
![EarthOrography_smooth](earth_orography_smooth.png)

## Load orography from file

The easiest to load another orography from a netCDF file is to reuse the
`EarthOrography`, e.g.

```julia
mars_orography = EarthOrography(spectal_grid, 
                                path="path/to/my/orography",
                                file="mars_orography.nc",
                                file_Grid=FullClenshawGrid)
```
the orography itself need to come on one of the full grids
SpeedyWeather defines, i.e. `FullGaussianGrid` or `FullClenshawGrid`
(a regular lat-lon grid, see [FullClenshawGrid](@ref FullClenshawGrid)),
which you can specify. Best to inspect the correct orientation with
`plot(mars_orography.orography)` (where the first is whatever name you chose here).
You can use smoothing as above.

## Changing orography manually

You can also change orography manually, that means by mutating the elements
in either `orography.orography` (to set it for the shallow-water model)
or `orography.geopot_surf` (for the primitive equations, but this is in
spectral space, advanced!). After the orography has been initialised
(for testing to `initialize!(orography, model)` but in most cases when 
you do `simulation = initialize!(model)`), for example you can do

```@example orography
sort!(orography.orography)
nothing # hide
```

to move all mountains to the south pole, or

```@example orography
orography.orography[1] = 100
```

to set the grid point `1` (at 0˚E, on the first ring around the north pole)
to a height of 100m. Whichever way you tweak the orography manually
you want to reflect this in the surface geopotential `geopot_surf` which is
used in the primitive equations by

```@example orography
transform!(orography.geopot_surf, orography.orography, model.spectral_transform)
orography.geopot_surf .*= model.planet.gravity
spectral_truncation!(orography.geopot_surf)
```

In the first line, the surface geopotential is still missing the gravity,
which is multiplied in the second line. The `spectral_truncation!`
only removes the ``l_{max}+1`` degree of the spherical harmonics that
scalar fields do not use, see [One more degree for spectral fields](@ref).

## Defining a new orography type

You can also define a new orography like we defined `ZonalRidge` or `EarthOrography`.
The following explains what's necessary for this. The new `MyOrography` has to be defined as
(`mutable` or not, but always with `@kwdef`)

```@example orography
@kwdef struct MyOrography{NF, Grid<:RingGrids.AbstractGrid{NF}} <: SpeedyWeather.AbstractOrography
    # optional, any parameters as fields here, e.g.
    constant_height::Float64 = 100
    # add some other parameters with default values

    # mandatory, every <:AbstractOrography needs those (same name, same type)
    orography::Grid                                 # in grid-point space [m]
    geopot_surf::LowerTriangularMatrix{Complex{NF}} # in spectral space *gravity [m^2/s^2]
end
```

for convenience with a generator function is automatically defined for all `AbstractOrography`

```@example orography
my_orography = MyOrography(spectral_grid, constant_height=200)
```

Now we have to extend the `initialize!` function. The first argument has to be
`::MyOrography` i.e. the new type we just defined, the second argument has to be
`::ModelSetup` although you could constrain it to `::ShallowWater` for example
but then it cannot be used for primitive equations.

```@example orography
function SpeedyWeather.initialize!(
    orog::MyOrography,      # first argument as to be ::MyOrography, i.e. your new type
    model::ModelSetup,      # second argument, use anything from model read-only
)
    (; orography, geopot_surf) = orog   # unpack

    # maybe use lat, lon coordinates (in degree or radians)
    (; latds, londs, lats, lons) = model.geometry

    # change here the orography grid [m], e.g.
    orography .= orography.constant_height

    # then also calculate the surface geopotential for primitive equations
    # given orography we just set
    transform!(geopot_surf, orography, model.spectral_transform)
    geopot_surf .*= model.planet.gravity
    spectral_truncation!(geopot_surf)
    return nothing
end
```

