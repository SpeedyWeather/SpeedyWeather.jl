full_grid_type(G::Type{<:AbstractFullGrid}) = G

## SIZE
get_nlon_max(Grid::Type{<:AbstractFullGrid}, nlat_half::Integer) = get_nlon(Grid, nlat_half)
get_nlon_per_ring(Grid::Type{<:AbstractFullGrid}, nlat_half::Integer, j::Integer) =
    get_nlon(Grid, nlat_half)
matrix_size(Grid::Type{<:AbstractFullGrid}, nlat_half::Integer) =
    (get_nlon(Grid, nlat_half), get_nlat(Grid, nlat_half))

## INDEXING
"""$(TYPEDSIGNATURES) `UnitRange` for every grid point of full grid `Grid` of resolution `nlat_half`
on ring `j` (`j=1` is closest ring around north pole, `j=nlat` around south pole)."""
function each_index_in_ring(
    Grid::Type{<:AbstractFullGrid},     # function for full grids
    j::Integer,                         # ring index north to south
    nlat_half::Integer,
)
    @boundscheck 0 < j <= get_nlat(Grid, nlat_half) || throw(BoundsError)    # valid ring index?
    nlon = get_nlon(Grid, nlat_half)    # number of longitudes per ring (const)
    index_1st = (j-1)*nlon + 1          # first in-ring index i
    index_end = j*nlon                  # last in-ring index i  
    return index_1st:index_end          # range of js in ring
end

# precompute ring indices for full grids
function each_index_in_ring!(   
    rings::AbstractVector,
    Grid::Type{<:AbstractFullGrid},
    nlat_half::Integer,
)
    nlat = length(rings)                # number of latitude rings
    @boundscheck nlat == get_nlat(Grid, nlat_half) || throw(BoundsError)

    nlon = get_nlon(Grid, nlat_half)    # number of longitudes
    index_end = 0                       
    @inbounds for j in 1:nlat
        index_1st = index_end + 1       # 1st index is +1 from prev ring's last index
        index_end += nlon               # only calculate last index per ring
        rings[j] = index_1st:index_end  # write UnitRange to rings vector
    end
end

## COORDINATES
function get_lon(Grid::Type{<:AbstractFullGrid}, nlat_half::Integer)
    lon = get_lond(Grid, nlat_half) # in degrees
    lon .*= 2π/360                  # convert to radians in-place
    return lon
end

function get_londlatds(Grid::Type{<:AbstractFullGrid}, nlat_half::Integer)

    lond = get_lond(Grid, nlat_half)        # vector of longitudes [0, 2π)
    latd = get_latd(Grid, nlat_half)        # vector of latitudes [90, -90]
    nlon = get_nlon(Grid, nlat_half)        # number of longitudes
    nlat = get_nlat(Grid, nlat_half)        # number of latitudes

    npoints = get_npoints(Grid, nlat_half)  # total number of grid points
    londs = zeros(npoints)
    latds = zeros(npoints)                  # preallocate

    for j in 1:nlat                         # populate preallocated colats, lons
        for i in 1:nlon
            ij = i + (j-1)*nlon             # continuous index ij
            londs[ij] = lond[i]
            latds[ij] = latd[j]
        end
    end

    return londs, latds
end