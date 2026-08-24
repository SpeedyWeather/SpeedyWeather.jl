# Visualization

SpeedyWeather.jl provides flexible visualization tools through its [RingGrids](https://github.com/SpeedyWeather/RingGrids.jl) module, which includes Makie and GeoMakie extensions for both 2D heatmap plots and interactive 3D globe visualizations.

## [Heatmap Plots](@id heatmap_plots)

The `heatmap` function provides a convenient way to visualize 2D fields on the sphere.

For full fields (which can be reshaped into matrices), creating a heatmap is straightforward:

```@example visualization
using SpeedyWeather
using CairoMakie, GeoMakie
using Random: Xoshiro

rng = Xoshiro(0)
grid = FullGaussianGrid(24)
field = randn(rng, grid)
heatmap(field)
```

The heatmap automatically:
- Sets the correct aspect ratio (2:1 for global maps)
- Labels longitude ticks in degrees East (0˚E, 60˚E, 120˚E, ...)
- Labels latitude ticks in degrees North (-60˚N, -30˚N, 0˚N, ...)
- Adds a colorbar
- Generates a descriptive title

For reduced grids (like `OctahedralGaussianGrid` or `HEALPixGrid`), the field is automatically interpolated to a full grid before plotting:

```@example visualization
grid = OctahedralGaussianGrid(24)
field = randn(rng, grid)
heatmap(field)
```

You can also pass additional keyword arguments to customize the appearance:

```@example visualization
grid = FullGaussianGrid(24)
field = randn(rng, grid)
heatmap(field; colormap = :thermal, colorrange = (-2, 2), title = "Temperature Anomaly")
```

Mutating variants of `heatmap` are also provided which can be used to plot into an externally defined `Figure` or `Axis`:

```@example visualization
grid = FullGaussianGrid(24)
fig = Figure(size = (800, 400))
field = randn(rng, grid)
ax, hm = heatmap!(fig[1, 1], field; title = "Time step 1")
field = randn(rng, grid)
heatmap!(ax, field; title = "Time step 2")
fig
```

Note that the mutating `heatmap!` returns different types depending on the type of the first argument:
- `heatmap!(pos::GridPosition, field; ...)` returns `(ax, hm)` where `ax` is the `Axis` and `hm` is the `HeatMap`
- `heatmap!(ax::Axis, field; ...)` returns the `HeatMap` object

## [Globe Plots](@id globe_plots)

For interactive 3D visualizations, SpeedyWeather.jl provides `globe` plots through the GeoMakie extension. These work with both `GLMakie` (interactive) and `CairoMakie` (static).

We can visualize the grid itself by passing the grid type and number of rings to `globe`:

```@example visualization
globe(FullGaussianGrid, 24; interactive = false)
```

Or plot data on the globe:

```@example visualization
grid = FullGaussianGrid(24)
field = randn(rng, grid)
globe(field; interactive = false)
```

The `globe` function accepts several keyword arguments:

```@example visualization
globe(
    OctahedralGaussianGrid, 24;
    interactive = false,
    title = "Octahedral Gaussian Grid",
    color = :blue,
    faces = true,
    centers = true,
    coastlines = true,
    background = true
)
```

Available options:
- `interactive::Bool` - Use 3D rotation/zoom (requires `GLMakie`)
- `title::String` - Plot title
- `color` - Color for grid points and lines
- `faces::Bool` - Show grid cell boundaries
- `centers::Bool` - Show grid points
- `coastlines::Bool` - Add coastlines
- `background::Bool` - Add Earth background image

Like heatmaps, globe plots have mutating variants that allow :

```@example visualization
fig = Figure(size = (800, 800))
ax, plotdata = globe!(fig[1, 1], FullGaussianGrid, 24; interactive = false)
```

`plotdata` is a NamedTuple with keys `:background`, `:centers`, `:faces`, and `:coastlines`.

The mutating `globe!` also returns different values depending on the arguments:
- `globe!(pos::GridPosition, Grid, nlat_half; ...)` returns `(ax, plotdata)` where `plotdata` contains the plotted elements
- `globe!(pos::GridPosition, field; ...)` returns `(ax, poly)` where `poly` is a `Poly` plot
- `globe!(ax::Axis, Grid, nlat_half; ...)` returns `plotdata` NamedTuple
- `globe!(ax::Axis, field; ...)` returns the `Poly` plot object

## Multi-panel `Figure`s

The mutating functions like `heatmap!` are especially useful for composing multiple plots into a single `Figure`, since it returns both the `Axis` and `HeatMap` needed to share components across panels.

```@example visualization
grid = OctahedralGaussianGrid(48)
u_field = randn(rng, grid)
v_field = randn(rng, grid)

fig = Figure(size = (1000, 300))
_, hm1 = heatmap!(fig[1, 1], u_field; title = "Zonal wind",     colormap = :balance, colorrange = (-3, 3))
_, hm2 = heatmap!(fig[1, 2], v_field; title = "Meridional wind", colormap = :balance, colorrange = (-3, 3))
Colorbar(fig[1, 3], hm, ticklabelsize = 10)
resize_to_layout!(fig)
```

## Tips

Use `heatmap!` and `globe!` when updating plots repeatedly—they skip figure and axis
allocation. When interpolation from a reduced grid happens many times (e.g. in an animation
loop), pre-interpolate once with `RingGrids.interpolate` and plot the resulting full field directly.

For the Makie backend, `GLMakie` supports interactive zoom and rotation in 3D globe plots,
`CairoMakie` produces static publication-quality figures, and `UnicodePlots` gives quick
terminal previews.

## API Reference

```@autodocs
Modules = [RingGrids]
Pages = ["RingGridsMakieExt.jl", "RingGridsGeoMakieExt.jl"]
```

## References

- [Makie.jl](https://makie.juliaplots.org/) - Flexible plotting library
- [GeoMakie.jl](https://makie.juliaplots.org/stable/backends/makie_extensions/geoplotting/) - Geographic plotting extension
- [RingGrids.jl](https://github.com/SpeedyWeather/RingGrids.jl) - Ring-based grid definitions
