module SpeedyWeatherGeoMakieExt

using SpeedyWeather
using GeoMakie, Geodesy

function _faces(geometry::Geometry{NF, <:OctaminimalGaussianArray}) where NF

    faces = zeros(geometry.nlon_max÷4+1, geometry.nlat+2, 2)
    fill!(faces, NaN)
    facesr = deepcopy(faces)

    for (j, latd) in enumerate(geometry.latd)
        nlon = RingGrids.get_nlon_per_ring(geometry.Grid, geometry.nlat_half, j)
        lond = RingGrids.get_lon_per_ring(geometry.Grid, geometry.nlat_half, j) * (360/2π)
        nlon ÷= 4
        dϕ = 90/nlon
            
        for i in 1:nlon
            faces[i,  j+1,  1] = lond[i] + dϕ/2
            faces[i,  j+1,  2] = latd
            facesr[i, j+1,  :] = faces[i, j+1, :]
        end

        facesr[nlon+1, 1, :] .= NaN      # remove to avoid overlapping grid lines at edges

        faces[end-1, j+1, 1] = 90
        faces[end-1, j+1, 2] = latd
    end

    # continue to the poles
    faces[:, 1, 1] = faces[:, 2, 1]
    faces[:, 1, 2] .= 90

    faces[:, end, 1] = faces[:, end-1, 1]
    faces[:, end, 2] .= -90;

    return faces, facesr
end

SpeedyWeather.globe(Grid::Type{<:AbstractGridArray}, nlat_half::Integer) = SpeedyWeather.globe(SpectralGrid(; Grid, nlat_half))
SpeedyWeather.globe(SG::SpectralGrid) = globe(Geometry(SG))

function SpeedyWeather.globe(geometry::Geometry)

    faces, facesr = _faces(geometry)

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