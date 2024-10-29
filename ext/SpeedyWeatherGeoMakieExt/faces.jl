function _faces(geometry::Geometry{NF, Grid}) where {NF, Grid<:OctaminimalGaussianArray}

    faces = zeros(geometry.nlon_max÷4+1, geometry.nlat+2, 2)
    fill!(faces, NaN)
    facesr = deepcopy(faces)

    for (j, latd) in enumerate(geometry.latd)
        nlon = RingGrids.get_nlon_per_ring(Grid, geometry.nlat_half, j)
        lond = RingGrids.get_lon_per_ring(Grid, geometry.nlat_half, j) * (360/2π)
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
    faces[:, end, 2] .= -90

    return faces, facesr
end

function _faces(geometry::Geometry{NF, Grid}) where {NF, Grid<:Union{OctahedralGaussianArray, OctahedralClenshawArray}}

    faces = zeros(geometry.nlon_max÷4+1, geometry.nlat+2, 2)
    fill!(faces, NaN)
    facesr = deepcopy(faces)

    for (j, latd) in enumerate(geometry.latd)
        nlon = RingGrids.get_nlon_per_ring(Grid, geometry.nlat_half, j)
        lond = RingGrids.get_lon_per_ring(Grid, geometry.nlat_half, j) * (360/2π)
        nlon ÷= 4
        dϕ = 90/nlon
            
        for i in 1:nlon+1
            faces[i,  j+1,  1] = lond[i] + dϕ/2
            faces[i,  j+1,  2] = latd
            facesr[i, j+1,  :] = faces[i, j+1, :]
        end

        facesr[nlon+1, 1, :] .= NaN      # remove to avoid overlapping grid lines at edges
        facesr[nlon+1, j+1, :] .= NaN    # remove to avoid overlapping grid lines at edges
    end

    # continue to the poles
    faces[:, 1, 1] = faces[:, 2, 1]
    faces[:, 1, 2] .= 90

    faces[:, end, 1] = faces[:, end-1, 1]
    faces[:, end, 2] .= -90
end