"""An `OctaHEALPixArray` is an array of OctaHEALPix grids, subtyping `AbstractReducedGridArray`.
First dimension of the underlying `N`-dimensional `data` represents the horizontal dimension,
in ring order (0 to 360˚E, then north to south), other dimensions are used for the vertical and/or
time or other dimensions. The resolution parameter of the horizontal grid is `nlat_half`
(number of latitude rings on one hemisphere, Equator included) and the ring indices are
precomputed in `rings`.

An OctaHEALPix grid has 4 faces, each `nlat_half x nlat_half` in size,
covering 90˚ in longitude, pole to pole. As part of the HEALPix family of grids,
the grid points are equal area. They start with 4 longitude points on the northern-most ring,
increase by 4 points per ring  towards the Equator with one ring on the Equator before reducing
the number of points again towards the south pole by 4 per ring. There is no equatorial belt for
OctaHEALPix grids. The southern hemisphere is symmetric to the northern, mirrored around the Equator.
OctaHEALPix grids have a ring on the Equator. For more details see
Górski et al. 2005, DOI:10.1086/427976, the OctaHEALPix grid belongs to the family of
HEALPix grids with Nθ = 1, Nφ = 4 but is not explicitly mentioned therein.

`rings` are the precomputed ring indices, for nlat_half = 3 (in contrast to HEALPix this can be odd)
it is `rings = [1:4, 5:12, 13:24, 25:32, 33:36]`. For efficient looping see `eachring` and `eachgrid`.
Fields are
$(TYPEDFIELDS)"""
struct OctaHEALPixArray{T, N, ArrayType <: AbstractArray{T, N}} <: AbstractReducedGridArray{T, N, ArrayType}
    data::ArrayType                 # data array, ring by ring, north to south
    nlat_half::Int                  # number of latitudes on one hemisphere
    rings::Vector{UnitRange{Int}}   # TODO make same array type as data?

    OctaHEALPixArray(data::A, nlat_half, rings) where {A <: AbstractArray{T, N}} where {T, N} =
        check_inputs(data, nlat_half, rings, OctaHEALPixArray) ?
        new{T, N, A}(data, nlat_half, rings) :
        error_message(data, nlat_half, rings, OctaHEALPixArray, T, N, A)
end

## TYPES
const OctaHEALPixGrid{T} = OctaHEALPixArray{T, 1, Vector{T}}
nonparametric_type(::Type{<:OctaHEALPixArray}) = OctaHEALPixArray
horizontal_grid_type(::Type{<:OctaHEALPixArray}) = OctaHEALPixGrid
full_array_type(::Type{<:OctaHEALPixArray}) = FullOctaHEALPixArray

"""An `OctaHEALPixArray` but constrained to `N=1` dimensions (horizontal only) and data is a `Vector{T}`."""
OctaHEALPixGrid

## SIZE
nlat_odd(::Type{<:OctaHEALPixArray}) = true
get_npoints2D(::Type{<:OctaHEALPixArray}, nlat_half::Integer) = 4*nlat_half^2
get_nlat_half(::Type{<:OctaHEALPixArray}, npoints2D::Integer) = round(Int, sqrt(npoints2D/4))

# number of longitude 
function get_nlon_per_ring(Grid::Type{<:OctaHEALPixArray}, nlat_half::Integer, j::Integer)
    nlat = get_nlat(Grid, nlat_half)
    @assert 0 < j <= nlat "Ring $j is outside P$nlat_half grid."
    # j = j > nlat_half ? nlat - j + 1 : j      # flip north south due to symmetry
    return min(4j, 8nlat_half-4j)
end

matrix_size(::Type{OctaHEALPixGrid}, nlat_half::Integer) = (2nlat_half, 2nlat_half)

## COORDINATES
function get_colat(::Type{<:OctaHEALPixArray}, nlat_half::Integer)
    nlat_half == 0 && return Float64[]
    colat = zeros(get_nlat(OctaHEALPixArray, nlat_half))
    for j in 1:nlat_half
        # Górski et al. 2005 eq 4 but without the 1/3 and Nside=nlat_half
        colat[j] = acos(1-(j/nlat_half)^2)  # northern hemisphere
        colat[2nlat_half-j] = π - colat[j]  # southern hemisphere
    end
    return colat
end

function get_lon_per_ring(Grid::Type{<:OctaHEALPixArray}, nlat_half::Integer, j::Integer)
    nlon = get_nlon_per_ring(Grid, nlat_half, j)
    # equidistant longitudes with equal offsets from 0˚ and 360˚,
    # e.g. 45, 135, 225, 315 for nlon=4
    return collect(π/nlon:2π/nlon:2π)
end

## INDEXING
function each_index_in_ring(::Type{<:OctaHEALPixArray},     # function for OctaHEALPix grids
                            j::Integer,                     # ring index north to south
                            nlat_half::Integer)             # resolution param

    @boundscheck 0 < j < 2nlat_half || throw(BoundsError)   # ring index valid?
    if j <= nlat_half                                       # northern hemisphere incl Equator
        index_1st = 2j*(j-1) + 1                            # first in-ring index i
        index_end = 2j*(j+1)                                # last in-ring index i
    else                                                    # southern hemisphere 
        n = 4nlat_half^2                                    # total number of points
        j = 2nlat_half - j                                  # count ring index from south pole
        index_1st = n - 2j*(j+1) + 1                        # count backwards
        index_end = n - 2j*(j-1)
    end
    return index_1st:index_end                              # range of i's in ring
end

function each_index_in_ring!(   rings::Vector{<:UnitRange{<:Integer}},
                                Grid::Type{<:OctaHEALPixArray},
                                nlat_half::Integer) # resolution param

    nlat = length(rings)
    @boundscheck nlat == get_nlat(Grid, nlat_half) || throw(BoundsError)

    index_end = 0
    @inbounds for j in 1:nlat_half                  # North incl Eq only
        index_1st = index_end + 1                   # 1st index is +1 from prev ring's last index
        index_end += 4j                             # add number of grid points per ring
        rings[j] = index_1st:index_end              # turn into UnitRange
    end
    @inbounds for (j, j_rev) in zip(nlat_half+1:nlat,       # South only
                                    nlat-nlat_half:-1:1)    # reverse index

        index_1st = index_end + 1                   # 1st index is +1 from prev ring's last index
        index_end += 4j_rev                         # add number of grid points per ring
        rings[j] = index_1st:index_end              # turn into UnitRange
    end
end

## CONVERSION
Base.Matrix(G::OctaHEALPixGrid{T}; kwargs...) where T = Matrix!(zeros(T, matrix_size(G)...), G; kwargs...)

"""
    Matrix!(M::AbstractMatrix,
            G::OctaHEALPixGrid;
            quadrant_rotation=(0, 1, 2, 3),
            matrix_quadrant=((2, 2), (1, 2), (1, 1), (2, 1)),
            )

Sorts the gridpoints in `G` into the matrix `M` without interpolation.
Every quadrant of the grid `G` is rotated as specified in `quadrant_rotation`,
0 is no rotation, 1 is 90˚ clockwise, 2 is 180˚ etc. Grid quadrants are counted
eastward starting from 0˚E. The grid quadrants are moved into the matrix quadrant
(i, j) as specified. Defaults are equivalent to centered at 0˚E and a rotation
such that the North Pole is at M's midpoint."""
Matrix!(M::AbstractMatrix, G::OctaHEALPixGrid; kwargs...) = Matrix!((M, G); kwargs...)

"""
    Matrix!(MGs::Tuple{AbstractMatrix{T}, OctaHEALPixGrid}...; kwargs...)

Like `Matrix!(::AbstractMatrix, ::OctaHEALPixGrid)` but for simultaneous
processing of tuples `((M1, G1), (M2, G2), ...)` with matrices `Mi` and grids `Gi`.
All matrices and grids have to be of the same size respectively."""
function Matrix!(   MGs::Tuple{AbstractMatrix{T}, OctaHEALPixGrid}...;
                    quadrant_rotation::NTuple{4, Integer}=(0, 1, 2, 3),     # = 0˚, 90˚, 180˚, 270˚ anti-clockwise
                    matrix_quadrant::NTuple{4, Tuple{Integer, Integer}}=((2, 2), (1, 2), (1, 1), (2, 1)),
                    ) where T
                    
    ntuples = length(MGs)

    # check that the first (matrix, grid) tuple has corresponding sizes
    M, G = MGs[1]
    m, n = size(M)
    @boundscheck m == n || throw(BoundsError)
    @boundscheck m == 2*G.nlat_half || throw(BoundsError)

    for MG in MGs   # check that all matrices and all grids are of same size
        Mi, Gi = MG
        @boundscheck size(Mi) == size(M) || throw(BoundsError)
        @boundscheck size(Gi) == size(G) || throw(BoundsError)
    end

    for q in matrix_quadrant    # check always in 2x2
        sr, sc = q
        @boundscheck ((sr in (1, 2)) && (sc in (1, 2))) || throw(BoundsError)
    end

    rings = eachring(G)         # index ranges for all rings
    nlat_half = G.nlat_half     # number of latitude rings on one hemisphere incl Equator
    nside = nlat_half           # side length of a basepixel matrix

    # sort grid indices from G into matrix M
    # 1) loop over each grid point per ring
    # 2) determine quadrant (0, 1, 2, 3) via modulo
    # 3) get longitude index iq within quadrant
    # 4) determine corresponding indices r, c in matrix M

    @inbounds for (j, ring) in enumerate(rings)
        nlon = length(ring)                         # number of grid points in ring
        for ij in ring                              # continuous index in grid
            i = ij-ring[1]                          # 0-based index in ring
            grid_quadrant = floor(Int, mod(4*i/nlon, 4))  # either 0, 1, 2, 3
            iq = i - grid_quadrant*(nlon÷4)             # 0-based index i relative to quadrant
            r = min(j, nlat_half) - iq               # row in matrix m (1-based)
            c = (iq+1) + max(0, j-nlat_half)         # column in matrix m (1-based)

            # rotate indices in quadrant
            r, c = rotate_matrix_indices(r, c, nside, quadrant_rotation[grid_quadrant+1])

            # shift grid quadrant to matrix quadrant
            sr, sc = matrix_quadrant[grid_quadrant+1]
            r += (sr-1)*nside                       # shift row into matrix quadrant
            c += (sc-1)*nside                       # shift column into matrix quadrant

            for (Mi, Gi) in MGs                      # for every (matrix, grid) tuple
                Mi[r, c] = convert(T, Gi[ij])         # convert data and copy over
            end
        end
    end

    ntuples == 1 && return M
    return Tuple(Mi for (Mi, Gi) in MGs)
end