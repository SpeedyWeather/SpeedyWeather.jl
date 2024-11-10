module SpeedyWeatherGeoMakieExt

using SpeedyWeather
using GeoMakie

using DocStringExtensions

include("faces.jl")

SpeedyWeather.globe(SG::SpectralGrid) = globe(SG.Grid, SG.nlat_half)
SpeedyWeather.globe(geometry::Geometry{NF, Grid}) where {NF, Grid} = globe(Grid, geometry.nlat_half)

"""($TYPEDSIGNATURES)
Create a 3D interactive globe plot of the grid `Grid` at resolution `nlat_half` displaying
cell centers and faces. Optionally, add coastlines and a background image of the Earth."""
function SpeedyWeather.globe(
    Grid::Type{<:AbstractGridArray},
    nlat_half::Integer;
    color = :black,
    faces::Bool = true,
    centers::Bool = true,
    coastlines::Bool = true,
    background::Bool = true,
)
    transf = GeoMakie.Geodesy.ECEFfromLLA(GeoMakie.Geodesy.WGS84())

    fig = Figure(size=(800, 800));
    ax = LScene(fig[1,1], show_axis=false);

    # background image
    if background
        bg = meshimage!(ax, -180..180, -90..90, rotr90(GeoMakie.earth()); npoints = 100, z_level = -10_000);
        bg.transformation.transform_func[] = transf
    end

    # cell centers, i.e. the grid points
    if centers
        latds, londs = RingGrids.get_latdlonds(Grid, nlat_half)
        c = scatter!(ax, londs, latds, markersize=5; color)
        c.transformation.transform_func[] = transf
    end

    # cell faces, a vector of Point2, concatenated all vertices for each grid point
    if faces
        # add nan after every face to avoid lines linking grid cells
        faces = get_faces(Grid, nlat_half, add_nan=true)
        f = lines!(ax, vec(faces); color)
        f.transformation.transform_func[] = transf
    end

    # coastlines
    if coastlines
        cl = lines!(GeoMakie.coastlines(50); color, linewidth=1)
        cl.transformation.transform_func[] = transf
    end

    # Makie stuff
    cc = cameracontrols(ax.scene)
    cc.settings.mouse_translationspeed[] = 0.0
    cc.settings.zoom_shift_lookat[] = false
    Makie.update_cam!(ax.scene, cc)

    return fig
end

import Makie.GeometryBasics: Polygon

"""($TYPEDSIGNATURES)
Create a 3D interactive globe plot of the data in `grid` displayed as polygons bounded by
the cell faces. Optionally, add coastlines (default true)."""
function SpeedyWeather.globe(
    grid::AbstractGrid;
    colormap = :viridis,
    coastlines::Bool = true,
)
    transf = GeoMakie.Geodesy.ECEFfromLLA(GeoMakie.Geodesy.WGS84())

    fig = Figure(size=(800, 800));
    ax = LScene(fig[1,1], show_axis=false);

    faces = get_faces(grid)
    polygons = [Polygon(faces[:, ij]) for ij in axes(faces, 2)]
    p = poly!(ax, polygons, color=grid.data; colormap)
    p.transformation.transform_func[] = transf

    if coastlines
        c = lines!(GeoMakie.coastlines(50); color=:white, linewidth=1, alpha=0.7)
        c.transformation.transform_func[] = transf
    end

    # Makie stuff
    cc = cameracontrols(ax.scene)
    cc.settings.mouse_translationspeed[] = 0.0
    cc.settings.zoom_shift_lookat[] = false
    Makie.update_cam!(ax.scene, cc)

    return fig
end

end # module