"""A `HEALPixArray` is an array of HEALPix grids, subtyping `AbstractReducedGridArray`.
First dimension of the underlying `N`-dimensional `data` represents the horizontal dimension,
in ring order (0 to 360˚E, then north to south), other dimensions are used for the vertical and/or
time or other dimensions. The resolution parameter of the horizontal grid is `nlat_half`
(number of latitude rings on one hemisphere, Equator included) which has to be even
(non-fatal error thrown otherwise) which is less strict than the original HEALPix
formulation (only powers of two for nside = nlat_half/2). Ring indices are precomputed in `rings`.

A HEALPix grid has 12 faces, each `nside`x`nside` grid points, each covering the same area
of the sphere. They start with 4 longitude points on the northern-most ring,
increase by 4 points per ring in the "polar cap" (the top half of the 4 northern-most faces)
but have a constant number of longitude points in the equatorial belt. The southern hemisphere
is symmetric to the northern, mirrored around the Equator. HEALPix grids have a ring on the
Equator. For more details see Górski et al. 2005, DOI:10.1086/427976. 

`rings` are the precomputed ring indices, for nlat_half = 4 it is
`rings = [1:4, 5:12, 13:20, 21:28, 29:36, 37:44, 45:48]`. So the first ring has indices
1:4 in the unravelled first dimension, etc. For efficient looping see `eachring` and `eachgrid`.
Fields are
$(TYPEDFIELDS)"""
struct HEALPixArray{T, N, ArrayType <: AbstractArray{T, N}} <: AbstractReducedGridArray{T, N, ArrayType}
    data::ArrayType                 # data array, ring by ring, north to south
    nlat_half::Int                  # number of latitudes on one hemisphere
    rings::Vector{UnitRange{Int}}   # TODO make same array type as data?

    HEALPixArray(data::A, nlat_half, rings) where {A <: AbstractArray{T, N}} where {T, N} =
        check_inputs(data, nlat_half, rings, HEALPixArray) ?
        new{T, N, A}(data, nlat_half, rings) :
        error_message(data, nlat_half, rings, HEALPixArray, T, N, A)
end

## TYPES
const HEALPixGrid{T} = HEALPixArray{T, 1, Vector{T}}
nonparametric_type(::Type{<:HEALPixArray}) = HEALPixArray
horizontal_grid_type(::Type{<:HEALPixArray}) = HEALPixGrid
full_array_type(::Type{<:HEALPixArray}) = FullHEALPixArray

"""An `HEALPixArray` but constrained to `N=1` dimensions (horizontal only) and data is a `Vector{T}`."""
HEALPixGrid

## SIZE
nlat_odd(::Type{<:HEALPixArray}) = true
get_npoints2D(::Type{<:HEALPixArray}, nlat_half::Integer) = 3*nlat_half^2
get_nlat_half(::Type{<:HEALPixArray}, npoints2D::Integer) = round(Int, sqrt(npoints2D/3))

"""$(TYPEDSIGNATURES) Number of longitude points for ring `j` on `Grid` of resolution
`nlat_half`."""
function get_nlon_per_ring(Grid::Type{<:HEALPixArray}, nlat_half::Integer, j::Integer)
    nlat = get_nlat(Grid, nlat_half)
    @assert 0 < j <= nlat "Ring $j is outside H$nlat_half grid."
    return min(4j, 2nlat_half, 8nlat_half-4j)
end

"""$(TYPEDSIGNATURES) The original `Nside` resolution parameter of the HEALPix grids.
The number of grid points on one side of each (square) face.
While we use `nlat_half` across all ring grids, this function translates this to
Nside. Even `nlat_half` only."""
function nside_healpix(nlat_half::Integer)
    iseven(nlat_half) || @error "Odd nlat_half=$nlat_half not supported for HEALPix."
    return nlat_half÷2
end

# for future reordering the HEALPix ring order into a matrix consisting of the
# 12 square faces of a HEALPix grid.
function matrix_size(::Type{<:HEALPixArray}, nlat_half::Integer)
    nside = nside_healpix(nlat_half)
    return (5nside, 5nside)
end

## COORDINATES
function get_colat(::Type{<:HEALPixArray}, nlat_half::Integer)
    nlat_half == 0 && return Float64[]
    
    nlat = get_nlat(HEALPixArray, nlat_half)
    nside = nside_healpix(nlat_half)
    colat = zeros(nlat)
    
    # Górski et al 2005 eq 4 and 8
    for j in 1:nside        colat[j] = acos(1-j^2/3nside^2)                 end     # north polar cap
    for j in nside+1:3nside colat[j] = acos(4/3 - 2j/3nside)                end     # equatorial belt
    for j in 3nside+1:nlat  colat[j] = acos((2nlat_half-j)^2/3nside^2-1)    end     # south polar cap

    return colat
end

function get_lon_per_ring(Grid::Type{<:HEALPixArray}, nlat_half::Integer, j::Integer)
    nside = nside_healpix(nlat_half)
    nlon = get_nlon_per_ring(Grid, nlat_half, j)

    # Górski et al 2005 eq 5 and 9
    # s = 1 for polar caps, s=2, 1, 2, 1, ... in the equatorial zone
    s = (j < nside) || (j >= 3nside) ? 1 : ((j - nside) % 2 + 1)
    lon = [π/(nlon÷2)*(i - s/2) for i in 1:nlon]
    return lon
end

## INDEXING
function each_index_in_ring(::Type{<:HEALPixArray},
                            j::Integer,                     # ring index north to south
                            nlat_half::Integer)             # resolution param
    nside = nside_healpix(nlat_half)

    @boundscheck 0 < j < 4nside || throw(BoundsError)       # ring index valid?
    if j < nside                                            # northern polar cap
        index_1st = 2j*(j-1) + 1                            # first in-ring index i
        index_end = 2j*(j+1)                                # last in-ring index i
    elseif j <= 3nside                                      # equatorial zone with const nlon
        n = 2nside^2-2nside                                 # points in polar cap
        nlon = 4nside                                       # points on latitude rings
        j = j - nside + 1                                   # offset ring index into eq zone
        index_1st = n + (j-1)*nlon + 1                      # add constant nlon per ring
        index_end = n + j*nlon
    else                                                    # southern polar cap
        n = 12nside^2                                       # total number of points
        j = 4nside - j                                      # count ring index from south pole
        index_1st = n - 2j*(j+1) + 1                        # count backwards
        index_end = n - 2j*(j-1)
    end
    return index_1st:index_end                              # range of i's in ring
end

function each_index_in_ring!(   rings::Vector{<:UnitRange{<:Integer}},
                                Grid::Type{<:HEALPixArray},
                                nlat_half::Integer) # resolution param

    nlat = length(rings)
    @boundscheck nlat == get_nlat(Grid, nlat_half) || throw(BoundsError)
    
    # HEALPix not defined for odd nlat_half, the last rings would not be written
    @boundscheck iseven(nlat_half) || throw(BoundsError)

    index_end = 0
    nside = nside_healpix(nlat_half)                # side length of a basepixel

    # North polar cap
    @inbounds for j in 1:nside-1
        index_1st = index_end + 1                   # 1st index is +1 from prev ring's last index
        index_end += 4j                             # add number of grid points per ring
        rings[j] = index_1st:index_end              # turn into UnitRange
    end

    # Equatorial belt
    nlon_max = get_nlon_max(Grid, nlat_half)        # number of grid points on belt
    @inbounds for j in nside:3nside
        index_1st = index_end + 1                   # 1st index is +1 from prev ring's last index
        index_end += nlon_max                       # nlon constant in belt
        rings[j] = index_1st:index_end              # turn into UnitRange
    end

    # South polar cap
    @inbounds for (j, j_mirrored) in zip(   3nside+1:nlat,  # South only
                                            nside-1:-1:1)   # mirror index

        index_1st = index_end + 1                   # 1st index is +1 from prev ring's last index
        index_end += 4j_mirrored                    # add number of grid points per ring
        rings[j] = index_1st:index_end              # turn into UnitRange
    end
end