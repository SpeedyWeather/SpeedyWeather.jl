module RingGridsGeoMakieExt

using RingGrids
using GeoMakie
using RingGrids.DocStringExtensions

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

"""($TYPEDSIGNATURES)
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

import GeoMakie.Makie.GeometryBasics: Polygon, Point

"""($TYPEDSIGNATURES)
Create a 3D interactive globe plot of the data in `grid` displayed as polygons bounded by
the cell faces. Optionally, add coastlines (default true)."""
function RingGrids.globe(
    field::AbstractField2D;
    interactive::Bool = true,
    title::String = default_title(field),
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
    polygons = [Polygon(Point.(faces[:, ij])) for ij in axes(faces, 2)]
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

end # module

