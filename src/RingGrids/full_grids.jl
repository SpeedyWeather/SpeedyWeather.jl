abstract type AbstractFullGridArray{T, N, ArrayType <: AbstractArray{T, N}} <: AbstractGridArray{T, N, ArrayType} end

"""
abstract type AbstractFullGrid{T} <: AbstractGrid{T} end

An `AbstractFullGrid` is a horizontal grid with a constant number of longitude
points across latitude rings. Different latitudes can be used, Gaussian latitudes,
equi-angle latitdes, or others."""
const AbstractFullGrid{T} = AbstractFullGridArray{T, 1, Vector{T}}

full_grid(G::Type{<:AbstractFullGridArray}) = G


get_nlon(Grid::Type{<:AbstractFullGrid}, nlat_half::Integer) = get_nlon_max(Grid, nlat_half)
get_nlon_max(::Type{<:AbstractFullGrid}, nlat_half::Integer) = 4nlat_half
get_nlon_per_ring(Grid::Type{<:AbstractFullGrid}, nlat_half::Integer, j::Integer) = 
    get_nlon(Grid, nlat_half)

function get_lon(Grid::Type{<:AbstractFullGrid}, nlat_half::Integer)
    nlat_half == 0 && return Float64[]
    nlon = get_nlon(Grid, nlat_half)
    return collect(range(0, 2π-π/nlon, step=2π/nlon))
end

function get_lond(Grid::Type{<:AbstractFullGrid}, nlat_half::Integer)
    lon = get_lon(Grid, nlat_half)
    lon .*= 360/2π      # convert to lond in-place
    return lon          # = lond
end

# convert an AbstractMatrix to the full grids, and vice versa
(Grid::Type{<:AbstractFullGrid})(M::AbstractMatrix{T}) where T = Grid{T}(vec(M))
Base.Matrix(grid::AbstractFullGrid{T}) where T = Matrix{T}(reshape(grid.data, :, get_nlat(grid)))
matrix_size(grid::AbstractFullGrid) = (get_nlon_max(grid), get_nlat(grid))
matrix_size(Grid::Type{<:AbstractFullGrid}, n::Integer) = (get_nlon_max(Grid, n), get_nlat(Grid, n))

function get_colatlons(Grid::Type{<:AbstractFullGrid}, nlat_half::Integer)

    colat = get_colat(Grid, nlat_half)       # vector of colats [0, π]
    lon = get_lon(Grid, nlat_half)           # vector of longitudes [0, 2π)
    nlon = get_nlon(Grid, nlat_half)         # number of longitudes
    nlat = get_nlat(Grid, nlat_half)         # number of latitudes

    npoints = get_npoints(Grid, nlat_half)   # total number of grid points
    colats = zeros(npoints)                 # preallocate
    lons = zeros(npoints)

    for j in 1:nlat                         # populate preallocated colats, lons
        for i in 1:nlon
            ij = i + (j-1)*nlon             # continuous index ij
            colats[ij] = colat[j]
            lons[ij] = lon[i]
        end
    end

    return colats, lons
end

function each_index_in_ring(
    Grid::Type{<:AbstractFullGridArray},    # function for full grids
    j::Integer,                             # ring index north to south
    nlat_half::Integer,
)
    @boundscheck 0 < j <= get_nlat(Grid, nlat_half) || throw(BoundsError)    # valid ring index?
    nlon = 4nlat_half               # number of longitudes per ring (const)
    index_1st = (j-1)*nlon + 1      # first in-ring index i
    index_end = j*nlon              # last in-ring index i  
    return index_1st:index_end      # range of js in ring
end

function each_index_in_ring!(   
    rings::Vector{<:UnitRange{<:Integer}},
    Grid::Type{<:AbstractFullGridArray},
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