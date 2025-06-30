module SpeedyWeatherGeoMakieExt

using SpeedyWeather
using GeoMakie
using SpeedyWeather.DocStringExtensions
using SpeedyWeather.NCDatasets

SpeedyWeather.globe(SG::SpectralGrid; kwargs...) = globe(SG.grid; kwargs...)
SpeedyWeather.globe(geometry::Geometry; kwargs...) = globe(geometry.grid; kwargs...)
SpeedyWeather.globe(grid::AbstractGrid; kwargs...) = globe(typeof(grid), grid.nlat_half; kwargs...)

"""($TYPEDSIGNATURES)
Create a 3D interactive globe plot of the grid `Grid` at resolution `nlat_half` displaying
cell centers and faces. Optionally, add coastlines and a background image of the Earth."""
function SpeedyWeather.globe(
    Grid::Type{<:AbstractGrid},
    nlat_half::Integer;
    interactive::Bool = true,
    title::String = "$(RingGrids.get_nlat(Grid, nlat_half))-ring $Grid",
    color = :black,
    faces::Bool = true,
    centers::Bool = true,
    coastlines::Bool = true,
    background::Bool = true,
)

    fig = Figure(size=(800, 800));

    if interactive
        transf = GeoMakie.Geodesy.ECEFfromLLA(GeoMakie.Geodesy.WGS84())
        ax = LScene(fig[1,1], show_axis=false);
    else
        ax = GeoAxis(fig[1, 1],
            title= title,
            dest = "+proj=ortho +lon_0=30 +lat_0=45")
    end

    # background image
    if background
        bg = meshimage!(ax, -180..180, -90..90, rotr90(GeoMakie.earth()); npoints = 100, z_level = -10_000);
        interactive && (bg.transformation.transform_func[] = transf)
    end

    # cell centers, i.e. the grid points
    if centers
        londs, latds = RingGrids.get_londlatds(Grid, nlat_half)
        c = scatter!(ax, londs, latds, markersize=5; color)
        interactive && (c.transformation.transform_func[] = transf)
    end

    # cell faces, a vector of NTuple{2, T}, concatenated all vertices for each grid point
    if faces
        # add nan after every face to avoid lines linking grid cells
        faces = RingGrids.get_gridcell_polygons(Grid, nlat_half, add_nan=true)
        f = lines!(ax, vec(faces); color)
        interactive && (f.transformation.transform_func[] = transf)
    end

    # coastlines
    if coastlines
        cl = lines!(GeoMakie.coastlines(50); color, linewidth=1)
        interactive && (cl.transformation.transform_func[] = transf)
    end

    # Makie stuff
    if interactive
        cc = cameracontrols(ax.scene)
        cc.settings.mouse_translationspeed[] = 0.0
        cc.settings.zoom_shift_lookat[] = false
        Makie.update_cam!(ax.scene, cc)
    else
        hidedecorations!(ax)
    end

    return fig
end

import GeoMakie.Makie.GeometryBasics: Polygon

"""($TYPEDSIGNATURES)
Create a 3D interactive globe plot of the data in `grid` displayed as polygons bounded by
the cell faces. Optionally, add coastlines (default true)."""
function SpeedyWeather.globe(
    field::AbstractField2D;
    interactive::Bool = true,
    title::String = "$(RingGrids.get_nlat(field))-ring $(typeof(field))",
    colormap = :viridis,
    coastlines::Bool = true,
)

    fig = Figure(size=(800, 800));
    
    if interactive
        transf = GeoMakie.Geodesy.ECEFfromLLA(GeoMakie.Geodesy.WGS84())
        ax = LScene(fig[1,1], show_axis=false);
    else
        ax = GeoAxis(fig[1, 1],
            title= title,
            dest = "+proj=ortho +lon_0=30 +lat_0=45")
    end

    faces = RingGrids.get_gridcell_polygons(field.grid)
    polygons = [Polygon(faces[:, ij]) for ij in axes(faces, 2)]
    p = poly!(ax, polygons, color=field.data; colormap)
    interactive && (p.transformation.transform_func[] = transf)

    if coastlines
        c = lines!(GeoMakie.coastlines(50); color=:white, linewidth=1, alpha=0.7)
        interactive && (c.transformation.transform_func[] = transf)
    end

    # Makie stuff
    if interactive
        cc = cameracontrols(ax.scene)
        cc.settings.mouse_translationspeed[] = 0.0
        cc.settings.zoom_shift_lookat[] = false
        Makie.update_cam!(ax.scene, cc)
    else
        hidedecorations!(ax)
    end

    return fig
end


"""
$(TYPEDSIGNATURES)
Create an animation from a SpeedyWeather simulation NetCDF output file. Needs a backend like 
CairoMakie or GtkMakie to be loaded at the time of calling this function.

# Arguments
- `nc_file::NCDataset`: NetCDF dataset containing SpeedyWeather simulation data
- `variable::String = "temp"`: Variable to animate (e.g., "temp", "pres", "vor", "u", "v", "humid")
- `level::Int = 1`: Vertical level to plot (for 3D variables)
- `transient_timesteps::Int = 0`: Number of timesteps to skip at the beginning of the animation
- `output_file::String = "animation.mp4"`: Path to save the animation
- `colormap = :viridis`: Colormap to use for the animation
- `framerate::Int = 15`: Frame rate of the animation
- `title::String = ""`: Title for the animation (if empty, will use variable name)
- `colorrange = nothing`: Range for the colorbar (nothing for automatic)
- `projection = "+proj=robin"`: Map projection to use (e.g., "+proj=robin", "+proj=ortho")
- `figure_size = (800, 600)`: Size of the figure
- `geoaxis_kwargs = (:xgridvisible => false, :ygridvisible => false)`: Keyword arguments to pass to the GeoAxis

# Example
```julia
using SpeedyWeather, GeoMakie, CairoMakie

# Animate temperature at level 1
animate("output.nc", variable="temp", level=1)

# Animate surface pressure with a specific colormap
animate(nc_file, variable="pres", colormap=:thermal)

# Create an animation with Orthographic projection
animate(nc_file, variable="temp", projection="+proj=ortho +lon_0=0 +lat_0=30")
```
"""
function SpeedyWeather.animate(
    ds::NCDataset;
    variable::String = "temp",
    level::Int = 1,
    transient_timesteps::Int = 0,
    output_file::String = "animation.mp4",
    colormap = :viridis,
    framerate::Int = 15,
    title::String = "",
    colorrange = nothing,
    projection = "+proj=robin",
    figure_size = (800, 600),
    coastlines::Bool = true,
    geoaxis_kwargs = (:xgridvisible => false, :ygridvisible => false, :titlesize => 20)
)
  
    # Get dimensions
    lon = ds["lon"][:]
    lat = ds["lat"][:]
    time = ds["time"][:]
    
    # Get time units for proper labeling
    time_units = ds["time"].attrib["units"]
    
    # Get variable metadata
    var_long_name = get(ds[variable].attrib, "long_name", variable)
    var_units = get(ds[variable].attrib, "units", "")
    
    # Set title if not provided
    if isempty(title)
        title = var_long_name
    end
    
    # Check if the variable is 3D (has a vertical dimension)
    is_3d = "layer" in dimnames(ds[variable])
    
    # Create the figure
    fig = Figure(size = figure_size)
    
    # Create the axis
    ax = GeoAxis(
            fig[1, 1],
            title = title,
            dest = projection;
            geoaxis_kwargs...
        )
    
    tsteps = Observable(transient_timesteps + 1)

    # Create data based on tsteps
    if is_3d
        # For 3D variables, extract the specified level
        data = @lift ds[variable][:, :, level, $tsteps]
    else
        # For 2D variables
        data = @lift ds[variable][:, :, $tsteps]
    end
    
    # Determine color range if not specified
    if isnothing(colorrange)
        # Sample data to determine a good color range
        sample_indices = min(10, length(time))

        sample_data = is_3d ? ds[variable][:, :, level, 1:sample_indices] : ds[variable][:, :, 1:sample_indices]
        
        # Calculate min and max across samples
        sample_data = filter(!isnan, sample_data)
        
        if !isempty(sample_data)
            data_min, data_max = extrema(sample_data)
            # Add a small buffer to the range
            range_buffer = (data_max - data_min) * 0.05
            colorrange = (data_min - range_buffer, data_max + range_buffer)
        else
            colorrange = (-1, 1)  # Fallback if all data is NaN
        end
    end
    
    # Create the surface plot
    hm = surface!(ax, lon, lat, data; colormap=colormap, colorrange=colorrange)
    
    # Add colorbar
    cb = Colorbar(fig[1, 2], hm, label = var_units)
    
    # Add time label
    time_label = Label(fig[2, 1:2], "Time: $(time[1]) $(time_units)")

    # Add coastlines 
    if coastlines
        lines!(GeoMakie.coastlines(); color=:white)
    end
    
    # Create the animation
    record(fig, output_file, (transient_timesteps + 1):length(time); framerate=framerate) do frame
        # Update the Observable
        tsteps[] = frame

        # Update the time label
        time_label.text = "Time: $(time[frame]) $(time_units)"
    end
    
    return output_file
end

"""
$(TYPEDSIGNATURES)
Create an animation from a NetCDF file. Needs a backend like 
CairoMakie or GtkMakie to be loaded at the time of calling this function.
Takes the same keyword arguments as [`SpeedyWeather.animate`](@ref).
"""
function SpeedyWeather.animate(
    file_path::String;
    kwargs...
)

    file = NCDataset(file_path)
    output_file = animate(file; kwargs...)
    close(file)

    return output_file
end

"""
$(TYPEDSIGNATURES)
Create an animation from a SpeedyWeather Simulation object. 
Needs the output of a simulation with NetCDF output enabled and a backend like 
CairoMakie or GtkMakie to be loaded at the time of calling this function.
Takes the same keyword arguments as [`SpeedyWeather.animate`](@ref).
"""
function SpeedyWeather.animate(
    simulation::Simulation;
    kwargs...
)
    # Extract the NetCDF file path from the simulation object
    if simulation.model.output.active
        nc_file = joinpath(simulation.model.output.run_path, simulation.model.output.filename)
    else
        error("NetCDF output is not active. Make sure the simulation has been run with NetCDF output enabled.")
    end
    
    # Check if the NetCDF file exists
    if !isfile(nc_file)
        error("NetCDF file $(nc_file) not found. Make sure the simulation has been run with NetCDF output enabled.")
    end
    
    # Call animate with the extracted file path
    return animate(
        nc_file;
        kwargs...
    )
end

end # module

