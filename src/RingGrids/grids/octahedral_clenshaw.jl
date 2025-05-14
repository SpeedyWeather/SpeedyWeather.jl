"""An `OctahedralClenshawArray` is an array of octahedral grids, subtyping `AbstractReducedGridArray`,
that use equidistant latitudes for each ring, the same as for `FullClenshawArray`.
First dimension of the underlying `N`-dimensional `data` represents the horizontal dimension,
in ring order (0 to 360˚E, then north to south), other dimensions are used for the vertical and/or
time or other dimensions. The resolution parameter of the horizontal grid is `nlat_half`
(number of latitude rings on one hemisphere, Equator included) and the ring indices are
precomputed in `rings`.

These grids are called octahedral (same as for the `OctahedralGaussianArray` which only uses different
latitudes) because after starting with 20 points on the first ring around the north pole (default) they
increase the number of longitude points for each ring by 4, such that they can be conceptually thought
of as lying on the 4 faces of an octahedron on each hemisphere. Hence, these grids have 20, 24, 28, ...
longitude points for ring 1, 2, 3, ... Clenshaw grids have a ring on the Equator which has 16 + 4nlat_half
longitude points before reducing the number of longitude points per ring by 4 towards the southern-most
ring `j = nlat`. `rings` are the precomputed ring indices, the the example above
`rings = [1:20, 21:44, 45:72, ...]`. For efficient looping see `eachring` and `eachgrid`.
Fields are
$(TYPEDFIELDS)"""
struct OctahedralClenshawGrid{A, V} <: AbstractReducedGrid
    nlat_half::Int                  # number of latitudes on one hemisphere
    architecture::A                 # information about device, CPU/GPU
    rings::V                        # precomputed ring indices
end

nonparametric_type(::Type{<:OctahedralClenshawGrid}) = OctahedralClenshawGrid
full_grid_type(::Type{<:OctahedralClenshawGrid}) = FullClenshawGrid

# SIZE
nlat_odd(::Type{<:OctahedralClenshawGrid}) = true
npoints_pole(::Type{<:OctahedralClenshawGrid}) = 16
npoints_added_per_ring(::Type{<:OctahedralClenshawGrid}) = 4

function get_npoints(::Type{<:OctahedralClenshawGrid}, nlat_half::Integer)
    m, o = npoints_added_per_ring(OctahedralClenshawGrid), npoints_pole(OctahedralClenshawGrid)
    return max(0, m*nlat_half^2 + 2o*nlat_half - o)     # to avoid negative for nlat_half = 0
end

function get_nlat_half(::Type{<:OctahedralClenshawGrid}, npoints::Integer)
    m, o = npoints_added_per_ring(OctahedralClenshawGrid), npoints_pole(OctahedralClenshawGrid)
    return round(Int, -o/m + sqrt(((o/m)^2 + (o+npoints)/m)))
end

function get_nlon_per_ring(Grid::Type{<:OctahedralClenshawGrid}, nlat_half::Integer, j::Integer)
    nlat = get_nlat(Grid, nlat_half)
    @assert 0 < j <= nlat "Ring $j is outside O$nlat_half grid."
    m, o = npoints_added_per_ring(OctahedralClenshawGrid), npoints_pole(OctahedralClenshawGrid)
    j = j > nlat_half ? nlat - j + 1 : j      # flip north south due to symmetry
    return o + m*j
end

matrix_size(grid::Grid) where {Grid<:OctahedralClenshawGrid} = matrix_size(Grid, grid.nlat_half)
function matrix_size(::Type{<:OctahedralClenshawGrid}, nlat_half::Integer)
    m, o = npoints_added_per_ring(OctahedralClenshawGrid), npoints_pole(OctahedralClenshawGrid)
    m != 4 && @warn "This algorithm has not been generalised for m!=4."
    N = (o + 4nlat_half)÷2
    return (N, N)
end

## COORDINATES
get_latd(::Type{<:OctahedralClenshawGrid}, nlat_half::Integer) = get_latd(FullClenshawGrid, nlat_half)
function get_lond_per_ring(Grid::Type{<:OctahedralClenshawGrid}, nlat_half::Integer, j::Integer)
    nlon = get_nlon_per_ring(Grid, nlat_half, j)
    return collect(0:360/nlon:360-180/nlon)
end

## QUADRATURE
get_quadrature_weights(::Type{<:OctahedralClenshawGrid}, nlat_half::Integer) =
    clenshaw_curtis_weights(nlat_half)

## INDEXING
function each_index_in_ring!(   rings::Vector{<:UnitRange{<:Integer}},
                                Grid::Type{<:OctahedralClenshawGrid},
                                nlat_half::Integer) # resolution param

    nlat = length(rings)
    @boundscheck nlat == get_nlat(Grid, nlat_half) || throw(BoundsError)
    m, o = npoints_added_per_ring(OctahedralClenshawGrid), npoints_pole(OctahedralClenshawGrid)

    index_end = 0
    @inbounds for j in 1:nlat_half                  # North incl Eq only
        index_1st = index_end + 1                   # 1st index is +1 from prev ring's last index
        index_end += o + m*j                        # add number of grid points per ring
        rings[j] = index_1st:index_end              # turn into UnitRange
    end
    @inbounds for (j, j_mirrored) in zip(   nlat_half+1:nlat,       # South only
                                            nlat-nlat_half:-1:1)    # reverse index

        index_1st = index_end + 1                   # 1st index is +1 from prev ring's last index
        index_end += o + m*j_mirrored               # add number of grid points per ring
        rings[j] = index_1st:index_end              # turn into UnitRange
    end
end

## CONVERSION
# Base.Matrix(G::OctahedralClenshawGrid{T}; kwargs...) where T =
#     Matrix!(zeros(T, matrix_size(G)...), G; kwargs...)

# """
# $(TYPEDSIGNATURES)
# Sorts the gridpoints in `G` into the matrix `M` without interpolation.
# Every quadrant of the grid `G` is rotated as specified in `quadrant_rotation`,
# 0 is no rotation, 1 is 90˚ clockwise, 2 is 180˚ etc. Grid quadrants are counted
# eastward starting from 0˚E. The grid quadrants are moved into the matrix quadrant
# (i, j) as specified. Defaults are equivalent to centered at 0˚E and a rotation
# such that the North Pole is at M's midpoint."""
# Matrix!(M::AbstractMatrix, G::OctahedralClenshawGrid; kwargs...) = Matrix!((M, G); kwargs...)

# """
# $(TYPEDSIGNATURES)
# Like `Matrix!(::AbstractMatrix, ::OctahedralClenshawGrid)` but for simultaneous
# processing of tuples `((M1, G1), (M2, G2), ...)` with matrices `Mi` and grids `Gi`.
# All matrices and grids have to be of the same size respectively."""
# function Matrix!(   MGs::Tuple{AbstractMatrix{T}, OctahedralClenshawGrid}...;
#                     quadrant_rotation::NTuple{4, Integer}=(0, 1, 2, 3),     # = 0˚, 90˚, 180˚, 270˚ anti-clockwise
#                     matrix_quadrant::NTuple{4, Tuple{Integer, Integer}}=((2, 2), (1, 2), (1, 1), (2, 1)),
#                     ) where T
    
#     # TODO make m, o dependent
#     m, o = npoints_added_per_ring(OctahedralGaussianArray), npoints_pole(OctahedralGaussianArray)
#     m != 4 || o != 16 && @warn "This algorithm has not been generalised for m!=4, o!=16."

#     ntuples = length(MGs)

#     # check that the first (matrix, grid) tuple has corresponding sizes
#     M, G = MGs[1]
#     m, n = size(M)
#     @boundscheck m == n || throw(BoundsError)
#     @boundscheck m == 2*(4+G.nlat_half) || throw(BoundsError)

#     for MG in MGs   # check that all matrices and all grids are of same size
#         Mi, Gi = MG
#         @boundscheck size(Mi) == size(M) || throw(BoundsError)
#         @boundscheck size(Gi) == size(G) || throw(BoundsError)
#     end

#     for q in matrix_quadrant    # check always in 2x2
#         sr, sc = q
#         @boundscheck ((sr in (1, 2)) && (sc in (1, 2))) || throw(BoundsError)
#     end

#     rings = eachring(G)         # index ranges for all rings
#     nlat_half = G.nlat_half     # number of latitude rings on one hemisphere incl Equator
#     nside = 4+G.nlat_half       # side length of a basepixel matrix

#     # sort grid indices from G into matrix M
#     # 1) loop over each grid point per ring
#     # 2) determine quadrant (0, 1, 2, 3) via modulo
#     # 3) get longitude index iq within quadrant
#     # 4) determine corresponding indices r, c in matrix M

#     @inbounds for (j, ring) in enumerate(rings)
#         nlon = length(ring)                         # number of grid points in ring
#         for ij in ring                              # continuous index in grid
#             i = ij-ring[1]                          # 0-based index in ring
#             grid_quadrant = floor(Int, mod(4*i/nlon, 4))  # either 0, 1, 2, 3
#             iq = i - grid_quadrant*(nlon÷4)         # 0-based index i relative to quadrant
#             r = 4+min(j, nlat_half) - iq            # row in matrix m (1-based)
#             c = (iq+1) + max(0, j-nlat_half)        # column in matrix m (1-based)

#             # rotate indices in quadrant
#             r, c = rotate_matrix_indices(r, c, nside, quadrant_rotation[grid_quadrant+1])

#             # shift grid quadrant to matrix quadrant
#             sr, sc = matrix_quadrant[grid_quadrant+1]
#             r += (sr-1)*nside                       # shift row into matrix quadrant
#             c += (sc-1)*nside                       # shift column into matrix quadrant

#             for (Mi, Gi) in MGs                     # for every (matrix, grid) tuple
#                 Mi[r, c] = convert(T, Gi[ij])       # convert data and copy over
#             end
#         end
#     end

#     ntuples == 1 && return M
#     return Tuple(Mi for (Mi, Gi) in MGs)
# end