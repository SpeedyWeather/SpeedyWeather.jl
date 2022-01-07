"""
Compute spectral geopotential `geopot` from spectral temperature `temp`
and spectral surface geopotential `geopot_surf` (orography*gravity).
"""
function geopotential!( geopot::Array{Complex{NF},3},       # spectral geopotential
                        geopot_surf::Array{Complex{NF},2},  # spectral surface geopotential
                        temp::Array{Complex{NF},3},         # spectral absolute Temperature
                        G::GeoSpectral{NF}) where {NF<:AbstractFloat}

    mx,nx,nlev = size(geopot)

    @boundscheck size(geopot) == size(temp) || throw(BoundsError())
    @boundscheck (mx,nx) == size(geopot_surf)   || throw(BoundsError())

    @unpack xgeop1, xgeop2, lapserate_correction = G.geometry

    # BOTTOM LAYER
    # is last index k=end, integration over half a layer
    for j in 1:nx
        for i in 1:mx
            geopot[i,j,end] = geopot_surf[i,j] + xgeop1[end]*temp[i,j,end]
        end
    end

    # OTHER LAYERS
    # integrate two half-layers from bottom to top
    for k in nlev-1:-1:1
        for j in 1:nx
            for i in 1:mx
                geopot[i,j,k] = geopot[i,j,k+1] + xgeop2[k+1]*temp[i,j,k+1] + xgeop1[k]*temp[i,j,k]
            end
        end
    end

    # LAPSERATE CORRECTION IN THE FREE TROPOSPHERE (>nlev)
    # TODO only for spectral coefficients 1,: ?
    for k in 2:nlev-1
        for j in 1:nx
            geopot[1,j,k] = geopot[1,j,k] +
                    lapserate_correction[k-1]*(temp[1,j,k+1] - temp[1,j,k-1])
        end
    end
end
