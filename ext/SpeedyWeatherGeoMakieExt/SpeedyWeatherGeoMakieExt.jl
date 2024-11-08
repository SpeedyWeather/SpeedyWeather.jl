module SpeedyWeatherGeoMakieExt

using SpeedyWeather
using GeoMakie

using DocStringExtensions

include("faces.jl")

SpeedyWeather.globe(SG::SpectralGrid) = globe(SG.Grid, SG.nlat_half)
SpeedyWeather.globe(geometry::Geometry{NF, Grid}) where {NF, Grid} = globe(Grid, geometry.nlat_half)

function SpeedyWeather.globe(Grid::Type{<:AbstractGridArray}, nlat_half::Integer)

    transf = GeoMakie.Geodesy.ECEFfromLLA(GeoMakie.Geodesy.WGS84())

    fig = Figure(size=(800, 800));
    ax = LScene(fig[1,1], show_axis=false);

    # background image
    bg = meshimage!(ax, -180..180, -90..90, rotr90(GeoMakie.earth()); npoints = 100, z_level = -10_000);
    bg.transformation.transform_func[] = transf

    # cell centers
    color = :black
    latds, londs = RingGrids.get_latdlonds(Grid, nlat_half)
    s = scatter!(ax, londs, latds, markersize=5; color)
    s.transformation.transform_func[] = transf

    # cell faces, a vector of Point2, concanated all vertices for each grid point
    faces = get_faces(Grid, nlat_half)
    lp = lines!(ax, vec(faces); color)

    # coastlines
    lpc = lines!(GeoMakie.coastlines(50); color, linewidth=1)
    lpc.transformation.transform_func[] = transf

    # Makie stuff
    cc = cameracontrols(ax.scene)
    cc.settings.mouse_translationspeed[] = 0.0
    cc.settings.zoom_shift_lookat[] = false
    Makie.update_cam!(ax.scene, cc)

    return fig
end

end # module