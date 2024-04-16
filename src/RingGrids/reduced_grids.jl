abstract type AbstractReducedGridArray{T, N, ArrayType <: AbstractArray{T, N}} <: AbstractGridArray{T, N, ArrayType} end
const AbstractReducedGrid{T} = AbstractReducedGridArray{T, 1, Vector{T}}

# all reduced grids have their maximum number of longitude points around the equator, i.e. j = nlat_half
get_nlon_max(Grid::Type{<:AbstractReducedGridArray}, nlat_half::Integer) = get_nlon_per_ring(Grid, nlat_half, nlat_half)

function get_colatlons(Grid::Type{<:AbstractReducedGridArray}, nlat_half::Integer)
    
    colat = get_colat(Grid, nlat_half)
    nlat = get_nlat(Grid, nlat_half)
    
    npoints = get_npoints(Grid, nlat_half)
    colats = zeros(npoints)                 # preallocate arrays
    lons = zeros(npoints)

    ij = 1                                  # continuous index
    for j in 1:nlat                         # populate arrays ring by ring
        lon = get_lon_per_ring(Grid, nlat_half, j)
        nlon = length(lon)

        colats[ij:ij+nlon-1] .= colat[j]
        lons[ij:ij+nlon-1] .= lon

        ij += nlon
    end

    return colats, lons
end