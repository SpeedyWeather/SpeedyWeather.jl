# Visualization

SpeedyWeather.jl provides flexible visualization tools through its [RingGrids](https://github.com/SpeedyWeather/RingGrids.jl) module, which includes Makie and GeoMakie extensions for both 2D heatmap plots and interactive 3D globe visualizations.

## [Heatmap Plots](@id heatmap_plots)

The `heatmap` function provides a convenient way to visualize 2D fields on the sphere. It automatically handles the transformation from the ring-based grid representation to a format suitable for plotting.

### Basic Usage

For full fields (which can be reshaped into matrices), creating a heatmap is straightforward:

```@example heatmap_basic
using SpeedyWeather
using CairoMakie   # or GLMakie for interactive plots

# Create a grid and field
grid = FullGaussianGrid(24)
field = randn(grid)

# Plot with automatic figure creation
fig = heatmap(field)
```

The heatmap automatically:
- Sets the correct aspect ratio (2:1 for global maps)
- Labels longitude ticks in degrees East (0˚E, 60˚E, 120˚E, ...)
- Labels latitude ticks in degrees North (-60˚N, -30˚N, 0˚N, ...)
- Adds a colorbar
- Generates a descriptive title

### Reduced Grids

For reduced grids (like `OctahedralGaussianGrid` or `HEALPixGrid`), the field is automatically interpolated to a full grid before plotting:

```@example heatmap_reduced
using SpeedyWeather
using CairoMakie

# Reduced grid
grid = OctahedralGaussianGrid(24)
field = randn(grid)

# Automatic interpolation to full grid for plotting
fig = heatmap(field)
```

### Customizing Heatmaps

You can pass additional keyword arguments to customize the appearance:

```@example heatmap_custom
using SpeedyWeather
using CairoMakie

grid = FullGaussianGrid(24)
field = randn(grid)

fig = heatmap(
    field;
    colormap = :thermal,
    colorrange = (-2, 2),
    title = "Temperature Anomaly"
)
```

### Mutating Variants for Efficient Updates

For animations or interactive applications where you need to update field data frequently, use the mutating variants `heatmap!`. These reuse the existing axis and only update the color data, which is much more efficient than creating new figures.

#### Updating with a GridPosition

```@example heatmap_mutating
using SpeedyWeather
using CairoMakie

# Create initial figure and axis
fig = Figure(size = (800, 400))

# First plot - creates the axis
grid = FullGaussianGrid(24)
field1 = randn(grid)
ax, hm = heatmap!(fig[1, 1], field1; title = "Time step 1")

# Update with new data - reuses the same axis
field2 = randn(grid)
heatmap!(fig[1, 1], field2; title = "Time step 2")

# Or update directly with an existing axis
field3 = randn(grid)
hm = heatmap!(ax, field3; title = "Time step 3")
```

#### Return Values

The mutating variants return useful objects:
- `heatmap!(pos::GridPosition, field; ...)` returns `(ax, hm)` where `ax` is the `Axis` and `hm` is the `HeatMap`
- `heatmap!(ax::Axis, field; ...)` returns the `HeatMap` object

This allows you to further customize the plot or extract the data later.

#### Animation Example

Here's an example of how to create an animation by updating field data:

```@example heatmap_animation
using SpeedyWeather
using CairoMakie

# Setup
grid = FullGaussianGrid(24)
fig = Figure(size = (800, 400))
ax, hm = heatmap!(fig[1, 1], randn(grid); title = "Frame 1")

# Update in a loop (would save frames in a real animation)
for t in 1:10
    field = randn(grid) * exp(-t * 0.1)
    heatmap!(ax, field; title = "Frame $t")
    # In practice: save(fig, "frame_$(lpad(t, 3, '0')).png")
end
```

## [Globe Plots](@id globe_plots)

For interactive 3D visualizations, SpeedyWeather.jl provides `globe` plots through the GeoMakie extension. These work with both `GLMakie` (interactive) and `CairoMakie` (static).

### Basic Globe Plots

Visualize the grid structure itself:

```@example globe_basic
using SpeedyWeather
using GLMakie, GeoMakie

# Interactive globe showing grid points and faces
globe(FullGaussianGrid, 24)
```

Or plot data on the globe:

```@example globe_field
using SpeedyWeather
using GLMakie, GeoMakie

grid = FullGaussianGrid(24)
field = randn(grid)

# Plot field as colored polygons
globe(field)
```

### Customization Options

The `globe` function accepts several keyword arguments:

```@example globe_custom
using SpeedyWeather
using CairoMakie, GeoMakie

grid = OctahedralGaussianGrid(24)

# Static globe with custom options
globe(
    grid, 24;
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

### Mutating Globe Variants

Like heatmaps, globe plots have mutating variants for efficient updates:

```@example globe_mutating
using SpeedyWeather
using GLMakie, GeoMakie

# Create initial globe
fig = Figure(size = (800, 800))
ax, objects = globe!(fig[1, 1], FullGaussianGrid, 24)

# Update with field data
grid = FullGaussianGrid(24)
field = randn(grid)
p = globe!(ax, field)

# The return values give you access to the plotted objects
# `objects` is a NamedTuple with keys: background, centers, faces, coastlines
# `p` is the PolyPlot object for the field
```

#### Return Values

- `globe!(pos::GridPosition, Grid, nlat_half; ...)` returns `(ax, objects)` where `objects` contains the plotted elements
- `globe!(pos::GridPosition, field; ...)` returns `(ax, poly_plot)`
- `globe!(ax::Axis, Grid, nlat_half; ...)` returns `objects` NamedTuple
- `globe!(ax::Axis, field; ...)` returns the `PolyPlot` object

## [Visualizing Simulation Output](@id sim_visualization)

When running SpeedyWeather simulations, you'll typically want to visualize the output variables. Here's a complete example:

```@example simulation_plot
using SpeedyWeather
using CairoMakie, GeoMakie

# Setup a simple simulation
arch = CPU()
spectral_grid = SpectralGrid(truncation = 32, dealiasing = 3)
model = PrimitiveDryModel(spectral_grid)
simulation = initialize!(model)

# Run for a few timesteps
run!(simulation, period = Hour(6))

# Extract and plot a variable
vars = simulation.model.vars
u_grid = vars.grid.u[:, 1]  # Surface zonal wind

# Create a figure with both heatmap and globe
fig = Figure(size = (1200, 400))

# Heatmap view
ax1, hm1 = heatmap!(fig[1, 1], u_grid; 
    title = "Zonal Wind (Heatmap)",
    colormap = :balance,
    colorrange = (-50, 50)
)

# Globe view  
ax2, p2 = globe!(fig[1, 2], RingGrids.Field(u_grid, spectral_grid.grid.grid);
    title = "Zonal Wind (Globe)",
    colormap = :balance
)

resize_to_layout!(fig)
```

## Working with Higher-Dimensional Fields

Fields can have multiple vertical layers and/or time steps. The heatmap functions automatically select the first slice for additional dimensions:

```@example higher_dims
using SpeedyWeather
using CairoMakie

grid = FullGaussianGrid(24)

# 3D field (horizontal + 10 vertical layers)
field_3d = randn(grid, 10)

# Automatically selects the first layer
fig = heatmap(field_3d)

# To plot a specific layer, index it first
fig = heatmap(field_3d[:, 5])  # 5th layer
```

For animations of simulation outputs over time, use the mutating variants with time-indexed slices:

```@example time_animation
using SpeedyWeather
using CairoMakie

grid = FullGaussianGrid(24)
field_4d = randn(grid, 5)  # 5 time steps

fig = Figure()
ax, hm = heatmap!(fig[1, 1], field_4d[:, 1]; title = "t = 1")

for t in 2:5
    heatmap!(ax, field_4d[:, t]; title = "t = $t")
end
```

## Performance Tips

1. **Use mutating variants for updates**: When creating animations or updating plots interactively, always use `heatmap!` or `globe!` instead of the non-mutating versions.

2. **Reuse figures**: Create a figure once and update its content rather than creating new figures in loops.

3. **Choose the right backend**: 
   - `GLMakie` for interactive exploration and 3D globes
   - `CairoMakie` for publication-quality static images
   - `UnicodePlots` for quick terminal visualization (limited functionality)

4. **Reduced grids are interpolated**: Plotting reduced grids automatically interpolates to full grids. For very high resolutions, consider pre-interpolating if you're making many plots.

## API Reference

```@autodocs
Modules = [RingGrids]
Pages = ["RingGridsMakieExt.jl", "RingGridsGeoMakieExt.jl"]
```

## References

- [Makie.jl](https://makie.juliaplots.org/) - Flexible plotting library
- [GeoMakie.jl](https://makie.juliaplots.org/stable/backends/makie_extensions/geoplotting/) - Geographic plotting extension
- [RingGrids.jl](https://github.com/SpeedyWeather/RingGrids.jl) - Ring-based grid definitions
