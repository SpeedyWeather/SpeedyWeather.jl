module SpeedyWeatherGeoMakieExt

using SpeedyWeather
using GeoMakie

include("faces.jl")

SpeedyWeather.globe(Grid::Type{<:AbstractGridArray}, nlat_half::Integer) = SpeedyWeather.globe(SpectralGrid(; Grid, nlat_half))
SpeedyWeather.globe(SG::SpectralGrid) = globe(Geometry(SG))

function SpeedyWeather.globe(geometry::Geometry)

    faces, facesr = _faces(geometry)

    transf = GeoMakie.Geodesy.ECEFfromLLA(GeoMakie.Geodesy.WGS84())

    fig = Figure(size=(800, 800));
    ax = LScene(fig[1,1], show_axis=false);

    # background image
    bg = meshimage!(ax, -180..180, -90..90, rotr90(GeoMakie.earth()); npoints = 100, z_level = -10_000);
    bg.transformation.transform_func[] = transf

    # cell centers
    color = :black
    s = scatter!(ax, geometry.londs, geometry.latds, markersize=5; color)
    s.transformation.transform_func[] = transf

    # cell faces
    for i in 1:geometry.nlon_max÷4
        for offset in [0, 90, 180, 270]
            lp1 = lines!(ax, facesr[i, :, 1] .+ offset, facesr[i, :, 2]; color)
            lp2 = lines!(ax, offset .- faces[i, :, 1], faces[i, :, 2]; color)

            lp1.transformation.transform_func[] = transf
            lp2.transformation.transform_func[] = transf
        end
    end

    # add equator
    lpe = lines!(ax, 0:360, zeros(361); color)
    lpe.transformation.transform_func[] = transf

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