module RingGridsGeoMakieExt

using RingGrids
using GeoMakie
using RingGrids.DocStringExtensions

import GeoMakie.Makie.GeometryBasics: Polygon, Point

RingGrids.globe(grid::AbstractGrid; kwargs...) = globe(typeof(grid), grid.nlat_half; kwargs...)

function default_title(Grid::Type{<:RingGrids.AbstractGrid}, nlat_half::Integer)
    Grid_ = RingGrids.nonparametric_type(Grid)
    return "$(RingGrids.get_nlat(Grid, nlat_half))-ring $Grid"
end

function default_title(field::RingGrids.AbstractField)
    Grid = RingGrids.nonparametric_type(field.grid)
    NF = eltype(field)
    return "$(RingGrids.get_nlat(field))-ring Field{$NF} on $Grid"
end

# ============================================================================
# Non-mutating globe variants
# ============================================================================

"""
$(TYPEDSIGNATURES)
Create a 3D interactive globe plot of the grid `Grid` at resolution `nlat_half` displaying
cell centers and faces. Optionally, add coastlines and a background image of the Earth."""
function RingGrids.globe(
        Grid::Type{<:AbstractGrid},
        nlat_half::Integer;
        interactive::Bool = true,
        title::String = default_title(Grid, nlat_half),
        color = :black,
        faces::Bool = true,
        centers::Bool = true,
        coastlines::Bool = true,
        background::Bool = true,
    )
    fig = Figure(size = (800, 800))
    globe!(fig[1, 1], Grid, nlat_half; interactive, title, color, faces, centers, coastlines, background)
    return fig
end

"""
$(TYPEDSIGNATURES)
Create a 3D interactive globe plot of the data in `field` displayed as polygons bounded by
the cell faces. Optionally, add coastlines (default true)."""
function RingGrids.globe(
        field::AbstractField2D;
        interactive::Bool = true,
        title::String = default_title(field),
        colormap = :viridis,
        coastlines::Bool = true,
    )
    fig = Figure(size = (800, 800))
    globe!(fig[1, 1], field; interactive, title, colormap, coastlines)
    return fig
end

# ============================================================================
# Mutating globe variants
# ============================================================================

"""
$(TYPEDSIGNATURES)
Mutating variant of `globe` for RingGrids grid types. Creates or updates a globe plot
in the given `GridPosition`. Returns the `Axis` and a NamedTuple of the plotted objects.
"""
function RingGrids.globe!(
        pos::Makie.GridPosition,
        Grid::Type{<:AbstractGrid},
        nlat_half::Integer;
        interactive::Bool = true,
        title::String = default_title(Grid, nlat_half),
        color = :black,
        faces::Bool = true,
        centers::Bool = true,
        coastlines::Bool = true,
        background::Bool = true,
        axis_kwargs = (;),
    )
    if interactive
        transf = GeoMakie.Geodesy.ECEFfromLLA(GeoMakie.Geodesy.WGS84())
        ax = LScene(pos; show_axis = false, axis_kwargs...)
    else
        ax = GeoAxis(
            pos;
            title,
            dest = "+proj=ortho +lon_0=30 +lat_0=45",
            axis_kwargs...
        )
    end

    objects = _plot_globe_grid!(
        ax, Grid, nlat_half; interactive, color, faces, centers, coastlines, background, transf
    )

    # Makie setup
    if interactive
        cc = cameracontrols(ax.scene)
        cc.settings.mouse_translationspeed[] = 0.0
        cc.settings.zoom_shift_lookat[] = false
        Makie.update_cam!(ax.scene, cc)
    else
        hidedecorations!(ax)
    end

    return ax, objects
end

"""
$(TYPEDSIGNATURES)
Mutating variant of `globe` for RingGrids `Field`s. Creates or updates a globe plot
in the given `GridPosition`. Returns the `Axis` and the `PolyPlot` object.
"""
function RingGrids.globe!(
        pos::Makie.GridPosition,
        field::AbstractField2D;
        interactive::Bool = true,
        title::String = default_title(field),
        colormap = :viridis,
        coastlines::Bool = true,
        axis_kwargs = (;),
    )
    if interactive
        transf = GeoMakie.Geodesy.ECEFfromLLA(GeoMakie.Geodesy.WGS84())
        ax = LScene(pos; show_axis = false, axis_kwargs...)
    else
        transf = nothing
        ax = GeoAxis(
            pos;
            title,
            dest = "+proj=ortho +lon_0=30 +lat_0=45",
            axis_kwargs...
        )
    end

    p = _plot_globe_field!(ax, field; interactive, colormap, coastlines, transf)

    # Makie setup
    if interactive
        cc = cameracontrols(ax.scene)
        cc.settings.mouse_translationspeed[] = 0.0
        cc.settings.zoom_shift_lookat[] = false
        Makie.update_cam!(ax.scene, cc)
    else
        hidedecorations!(ax)
    end

    return ax, p
end

"""
$(TYPEDSIGNATURES)
Mutating variant of `globe` that plots directly into an existing `Axis`.
Returns a NamedTuple of the plotted objects.
"""
function RingGrids.globe!(
        ax::Makie.AbstractAxis,
        Grid::Type{<:AbstractGrid},
        nlat_half::Integer;
        interactive::Bool = true,
        color = :black,
        faces::Bool = true,
        centers::Bool = true,
        coastlines::Bool = true,
        background::Bool = true,
    )
    if interactive
        transf = GeoMakie.Geodesy.ECEFfromLLA(GeoMakie.Geodesy.WGS84())
    else
        transf = nothing
    end

    objects = _plot_globe_grid!(
        ax, Grid, nlat_half; interactive, color, faces, centers, coastlines, background, transf
    )

    return objects
end

"""
$(TYPEDSIGNATURES)
Mutating variant of `globe` that plots directly into an existing `Axis`.
Returns the `PolyPlot` object.
"""
function RingGrids.globe!(
        ax::Makie.AbstractAxis,
        field::AbstractField2D;
        interactive::Bool = true,
        colormap = :viridis,
        coastlines::Bool = true,
    )
    if interactive
        transf = GeoMakie.Geodesy.ECEFfromLLA(GeoMakie.Geodesy.WGS84())
    else
        transf = nothing
    end

    p = _plot_globe_field!(ax, field; interactive, colormap, coastlines, transf)

    return p
end

# ============================================================================
# Helper functions for globe plotting
# ============================================================================

"""
    _plot_globe_grid!(ax, Grid, nlat_half; kwargs...)

Internal helper to plot grid cells and points on a globe. Returns a NamedTuple
of the plotted objects.
"""
function _plot_globe_grid!(
        ax,
        Grid::Type{<:AbstractGrid},
        nlat_half::Integer;
        interactive::Bool,
        color,
        faces::Bool,
        centers::Bool,
        coastlines::Bool,
        background::Bool,
        transf,
    )
    objects = NamedTuple()

    # background image
    if background
        bg = meshimage!(ax, -180 .. 180, -90 .. 90, rotr90(GeoMakie.earth()); npoints = 100, z_level = -10_000)
        interactive && (bg.transformation.transform_func[] = transf)
        objects = merge(objects, (; background = bg))
    end

    # cell centers, i.e. the grid points
    if centers
        londs, latds = RingGrids.get_londlatds(Grid, nlat_half)
        c = scatter!(ax, londs, latds, markersize = 5; color)
        interactive && (c.transformation.transform_func[] = transf)
        objects = merge(objects, (; centers = c))
    end

    # cell faces, a vector of NTuple{2, T}, concatenated all vertices for each grid point
    if faces
        # add nan after every face to avoid lines linking grid cells
        grid_faces = RingGrids.get_gridcell_polygons(Grid, nlat_half, add_nan = true)
        f = lines!(ax, vec(grid_faces); color)
        interactive && (f.transformation.transform_func[] = transf)
        objects = merge(objects, (; faces = f))
    end

    # coastlines
    if coastlines
        cl = lines!(GeoMakie.coastlines(50); color, linewidth = 1)
        interactive && (cl.transformation.transform_func[] = transf)
        objects = merge(objects, (; coastlines = cl))
    end

    return objects
end

"""
    _plot_globe_field!(ax, field; kwargs...)

Internal helper to plot a field on a globe as polygons. Returns the `PolyPlot` object.
"""
function _plot_globe_field!(
        ax,
        field::AbstractField2D;
        interactive::Bool,
        colormap,
        coastlines::Bool,
        transf,
    )
    faces = RingGrids.get_gridcell_polygons(field.grid)
    polygons = [Polygon(Point.(faces[:, ij])) for ij in axes(faces, 2)]
    p = poly!(ax, polygons, color = field.data; colormap)
    interactive && (p.transformation.transform_func[] = transf)

    if coastlines
        c = lines!(GeoMakie.coastlines(50); color = :white, linewidth = 1, alpha = 0.7)
        interactive && (c.transformation.transform_func[] = transf)
    end

    return p
end

end # module
