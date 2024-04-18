"""Subtype of `AbstractGridArray` for arrays of rings grids that have a reduced number
of longitude points towards the poles, i.e. they are not "full", see `AbstractFullGridArray`.
Data on these grids cannot be represented as matrix and has to be unravelled into a vector,
ordered 0 to 360˚E then north to south, ring by ring. Examples for reduced grids are 
the octahedral Gaussian or Clenshaw grids, or the HEALPix grid."""
abstract type AbstractReducedGridArray{T, N, ArrayType <: AbstractArray{T, N}} <: AbstractGridArray{T, N, ArrayType} end

"""Horizontal abstract type for all `AbstractReducedGridArray` with `N=1` (i.e. horizontal only)
and `ArrayType` of `Vector{T}` with element type `T`."""
const AbstractReducedGrid{T} = AbstractReducedGridArray{T, 1, Vector{T}}

# only defined for full grids
get_lon(::Type{<:AbstractReducedGridArray}, nlat_half::Integer) = Float64[]
get_lond(::Type{<:AbstractReducedGridArray}, nlat_half::Integer) = Float64[]

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