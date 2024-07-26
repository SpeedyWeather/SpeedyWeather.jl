## Gradient operators

[SpeedyTransforms](@ref) also includes many gradient operators to take derivatives in
spherical harmonics. These are in particular ``\nabla, \nabla \cdot, \nabla \times,
\nabla^2, \nabla^{-2}``. We call them `divergence`, `curl`, `∇`, `∇²`, `∇⁻²`
(as well as their in-place versions with `!`) within the limits of unicode characters
and Julia syntax. These functions are defined for inputs being spectral coefficients
(i.e. `LowerTriangularMatrix`) or gridded fields (i.e. `<:AbstractGrid`) and
also allow as an additional argument a spectral transform object
(see [SpectralTransform](@ref SpectralTransform)) which avoids recalculating it
under the hood.

!!! info "SpeedyTransforms assumes a unit sphere"
    The gradient operators in SpeedyTransforms generally assume a sphere of radius ``R=1``.
    For the transforms themselves that does not make a difference, but the gradient operators
    `divergence`, `curl`, `∇`, `∇²`, `∇⁻²` omit the radius scaling unless you provide the optional
    keyword `radius` (or you can do `./= radius` manually). 
    Also note that meridional derivates in spectral space expect a ``\cos^{-1}(\theta)`` scaling.
    Details are always outlined in the respective docstrings, `?∇` for example.


The actually implemented operators are,
in contrast to the mathematical [Derivatives in spherical coordinates](@ref)
due to reasons of scaling as follows. Let the implemented operators be
``\hat{\nabla}`` etc.

```math
\hat{\nabla} A = \left(\frac{\partial A}{\partial \lambda}, \cos(\theta)\frac{\partial A}{\partial \theta} \right) =
R\cos(\theta)\nabla A
```
So the zonal derivative omits the radius and the ``\cos^{-1}(\theta)`` scaling.
The meridional derivative adds a ``\cos(\theta)`` due to a recursion relation
being defined that way, which, however, is actually convenient because the whole
operator is therefore scaled by ``R\cos(\theta)``. The curl and divergence operators
expect the input velocity fields to be scaled by ``\cos^{-1}(\theta)``, i.e.

```math
\begin{aligned}
\hat{\nabla} \cdot (\cos^{-1}(\theta)\mathbf{u}) &= \frac{\partial u}{\partial \lambda} +
\cos\theta\frac{\partial v}{\partial \theta} = R\nabla \cdot \mathbf{u}, \\
\hat{\nabla} \times (\cos^{-1}(\theta)\mathbf{u}) &= \frac{\partial v}{\partial \lambda} -
\cos\theta\frac{\partial u}{\partial \theta} = R\nabla \times \mathbf{u}.
\end{aligned}
```

And the Laplace operators omit a ``R^2`` (radius ``R``) scaling, i.e.

```math
\hat{\nabla}^{-2}A = \frac{1}{R^2}\nabla^{-2}A , \quad \hat{\nabla}^{2}A = R^2\nabla^{2}A
```
## Gradient `∇`

We illustrate the usage of the gradient function `∇`. Let us create some fake
data `G` on the grid first 

```@example gradient
using SpeedyWeather, CairoMakie

# create some data with wave numbers 0,1,2,3,4
trunc = 64                  # 1-based maximum degree of spherical harmonics
L = randn(LowerTriangularMatrix{ComplexF32}, trunc, trunc)
spectral_truncation!(L, 5)              # remove higher wave numbers
G = transform(L)
heatmap(G, title="Some fake data G")    # requires `using CairoMakie`
save("gradient_data.png", ans) # hide
nothing # hide
```
![Gradient data](gradient_data.png)

Now we can take the gradient as follows
```@example gradient
dGdx, dGdy = ∇(G)
nothing # hide
```
this transforms internally back to spectral space takes the gradients in
zonal and meridional direction, transforms to grid-point space again
und unscales the coslat-scaling on the fly but assumes a radius of 1
as the keyword argument `radius` was not provided.
Use `∇(G, radius=6.371e6)` for a gradient on Earth in units of "data unit"
divided by meters.

```@example gradient
heatmap(dGdx, title="dG/dx on the unit sphere")
save("dGdx.png", ans) # hide
heatmap(dGdy, title="dG/dy on the unit sphere")
save("dGdy.png", ans) # hide
nothing # hide
```
![dGdx](dGdx.png)
![dGdy](dGdy.png)

## Geostrophy

Now, we want to use the following example to illustrate a more complex
use of the gradient operators: We have ``u, v`` and want to calculate
``\eta`` in the shallow water system from it following geostrophy.
Analytically we have
```math
-fv = -g\partial_\lambda \eta, \quad fu = -g\partial_\theta \eta
```
which becomes, if you take the divergence of these two equations
```math
\zeta = \frac{g}{f}\nabla^2 \eta
```
Meaning that if we start with ``u, v`` we can obtain the relative vorticity
``\zeta`` and, using Coriolis parameter ``f`` and gravity ``g``, invert
the Laplace operator to obtain displacement ``\eta``. How to do this with
SpeedyTransforms? 

Let us start by generating some data
```@example gradient
spectral_grid = SpectralGrid(trunc=31, nlayers=1)
forcing = SpeedyWeather.JetStreamForcing(spectral_grid)
drag = QuadraticDrag(spectral_grid)
model = ShallowWaterModel(; spectral_grid, forcing, drag)
model.feedback.verbose = false # hide
simulation = initialize!(model);
run!(simulation, period=Day(30))
nothing # hide
```

Now pretend you only have `u, v` to get vorticity (which is actually the prognostic variable in the model,
so calculated anyway...).
```@example gradient
u = simulation.diagnostic_variables.grid.u_grid[:,1]
v = simulation.diagnostic_variables.grid.v_grid[:,1]
vor = curl(u, v, radius = spectral_grid.radius)
nothing # hide
```
Here, `u, v` are the grid-point velocity fields, and the function `curl` takes in either
`LowerTriangularMatrix`s (no transform needed as all gradient operators act in spectral space),
or, as shown here, arrays of the same grid and size. In this case, the function actually
runs through the following steps
```@example gradient
RingGrids.scale_coslat⁻¹!(u)
RingGrids.scale_coslat⁻¹!(v)

S = SpectralTransform(u, one_more_degree=true)
us = transform(u, S)
vs = transform(v, S)

vor = curl(us, vs, radius = spectral_grid.radius)
```
(Copies of) the velocity fields are unscaled by the cosine of latitude (see above),
then transformed into spectral space, and the `curl` has the keyword argument `radius`
to divide internally by the radius (if not provided it assumes a unit sphere).
We always unscale vector fields by the cosine of latitude if they are provided to
`curl` or `divergence` in spectral as you can only do this scaling
effectively in grid-point space. The methods accepting arguments as grids generally
do this for you. If in doubt, check the docstrings, `?∇` for example.

## One more degree for spectral fields

The `SpectralTransform` in general takes a `one_more_degree` keyword argument,
if otherwise the returned `LowerTriangularMatrix` would be of size 32x32, setting this
to true would return 33x32. The reason is that while most people would expect
square lower triangular matrices for a triangular spectral truncation, all vector quantities
always need one more degree (= one more row) because of a recursion relation in the
meridional gradient. So as we want to take the curl of `us, vs` here, they need this
additional degree, but in the returned lower triangular matrix this row is set to zero.

!!! info "One more degree for vector quantities"
    All gradient operators expect the input lower triangular matrices of shape ``(N+1) \times N``.
    This one more degree of the spherical harmonics is required for the meridional derivative.
    Scalar quantities contain this degree too for size compatibility but they should not
    make use of it. Use `spectral_truncation` to add or remove this degree manually.

## Example: Geostrophy (continued)

Now we transfer `vor` into grid-point space, but specify that we want it on the grid
that we also used in `spectral_grid`. The Coriolis parameter for a grid like `vor_grid`
is obtained, and we do the following for ``f\zeta/g``.

```@example gradient
vor_grid = transform(vor, Grid=spectral_grid.Grid)
f = coriolis(vor_grid)      # create Coriolis parameter f on same grid with default rotation
g = model.planet.gravity
fζ_g = @. vor_grid * f / g  # in-place and element-wise
nothing # hide
```
Now we need to apply the inverse Laplace operator to ``f\zeta/g`` which we do as follows

```@example gradient
fζ_g_spectral = transform(fζ_g, one_more_degree=true)

R = spectral_grid.radius
η = SpeedyTransforms.∇⁻²(fζ_g_spectral) * R^2
η_grid = transform(η, Grid=spectral_grid.Grid)
nothing # hide
```
Note the manual scaling with the radius ``R^2`` here. We now compare the results
```@example gradient
using CairoMakie
heatmap(η_grid, title="Geostrophic interface displacement η [m]")
save("eta_geostrophic.png", ans) # hide
nothing # hide
```
![Geostrophic eta](eta_geostrophic.png)

Which is the interface displacement assuming geostrophy.
The actual interface displacement contains also ageostrophy
```@example gradient
η_grid2 = simulation.diagnostic_variables.grid.pres_grid
heatmap(η_grid2, title="Interface displacement η [m] with ageostrophy")
save("eta_ageostrophic.png", ans) # hide
nothing # hide
```
![Ageostrophic eta](eta_ageostrophic.png)

Strikingly similar! The remaining differences are the ageostrophic motions but
also note that the mean is off. This is because geostrophy only use/defines the gradient
of ``\eta`` not the absolute values itself. Our geostrophic ``\eta_g`` has by construction
a mean of zero (that is how we define the inverse Laplace operator) but the actual ``\eta``
is some 1400m higher.
