"""
    G = OctahedralClenshawGrid{T}

An Octahedral Clenshaw grid that uses `nlat` equi-spaced latitudes. Like FullClenshawGrid, the central
latitude ring is on the Equator. Like OctahedralGaussianGrid, the number of longitude points per
latitude ring decreases towards the poles. Starting with 20 equi-spaced longitude points (starting at 0˚E)
on the rings around the poles, each latitude ring towards the equator has consecuitively 4 more points,
one for each face of the octahedron. E.g. 20, 24, 28, 32, ...nlon-4, nlon, nlon, nlon-4, ..., 32, 28, 24, 20.
The maximum number of longitue points is `nlon`. The values of all grid points are stored in a vector
field `v` that unravels the data 0 to 360˚, then ring by ring, which are sorted north to south."""
struct OctahedralClenshawGrid{T} <: AbstractOctahedralGrid{T}
    data::Vector{T}    # data vector, ring by ring, north to south
    nlat_half::Int  # number of latitudes on one hemisphere (incl Equator)

    # check that `nlat_half` match the vector `v` length
    OctahedralClenshawGrid{T}(data::AbstractVector, nlat_half::Integer) where T = length(data) == npoints_octahedral(nlat_half, true) ?
    new(data, nlat_half) : error("$(length(data))-element Vector{$(eltype(data))}"*
    "cannot be used to create a O$(nlat_half) OctahedralClenshawGrid{$T}.")
end

nonparametric_type(::Type{<:OctahedralClenshawGrid}) = OctahedralClenshawGrid

# infer nlat_half from data vector length, infer parametric type from eltype of data
OctahedralClenshawGrid{T}(data::AbstractVector) where T = OctahedralClenshawGrid{T}(data,
                                                        nlat_half_octahedral(length(data), true))
OctahedralClenshawGrid(data::AbstractVector, n::Integer...) = OctahedralClenshawGrid{eltype(data)}(data, n...)

nlat_odd(::Type{<:OctahedralClenshawGrid}) = true
get_npoints(::Type{<:OctahedralClenshawGrid}, nlat_half::Integer) = npoints_octahedral(nlat_half, true)
get_colat(::Type{<:OctahedralClenshawGrid}, nlat_half::Integer) = get_colat(FullClenshawGrid, nlat_half)
get_quadrature_weights(::Type{<:OctahedralClenshawGrid}, nlat_half::Integer) = clenshaw_curtis_weights(nlat_half)
full_grid(::Type{<:OctahedralClenshawGrid}) = FullClenshawGrid    # the full grid with same latitudes

matrix_size(G::OctahedralClenshawGrid) = (2*(4+G.nlat_half), 2*(4+G.nlat_half))
matrix_size(::Type{OctahedralClenshawGrid}, nlat_half::Integer) = (2*(4+nlat_half), 2*(4+nlat_half))
Base.Matrix(G::OctahedralClenshawGrid{T}; kwargs...) where T = Matrix!(zeros(T, matrix_size(G)...), G; kwargs...)

"""
    Matrix!(M::AbstractMatrix,
            G::OctahedralClenshawGrid;
            quadrant_rotation=(0, 1, 2, 3),
            matrix_quadrant=((2, 2), (1, 2), (1, 1), (2, 1)),
            )

Sorts the gridpoints in `G` into the matrix `M` without interpolation.
Every quadrant of the grid `G` is rotated as specified in `quadrant_rotation`,
0 is no rotation, 1 is 90˚ clockwise, 2 is 180˚ etc. Grid quadrants are counted
eastward starting from 0˚E. The grid quadrants are moved into the matrix quadrant
(i, j) as specified. Defaults are equivalent to centered at 0˚E and a rotation
such that the North Pole is at M's midpoint."""
Matrix!(M::AbstractMatrix, G::OctahedralClenshawGrid; kwargs...) = Matrix!((M, G); kwargs...)

"""
    Matrix!(MGs::Tuple{AbstractMatrix{T}, OctahedralClenshawGrid}...; kwargs...)

Like `Matrix!(::AbstractMatrix, ::OctahedralClenshawGrid)` but for simultaneous
processing of tuples `((M1, G1), (M2, G2), ...)` with matrices `Mi` and grids `Gi`.
All matrices and grids have to be of the same size respectively."""
function Matrix!(   MGs::Tuple{AbstractMatrix{T}, OctahedralClenshawGrid}...;
                    quadrant_rotation::NTuple{4, Integer}=(0, 1, 2, 3),     # = 0˚, 90˚, 180˚, 270˚ anti-clockwise
                    matrix_quadrant::NTuple{4, Tuple{Integer, Integer}}=((2, 2), (1, 2), (1, 1), (2, 1)),
                    ) where T
                    
    ntuples = length(MGs)

    # check that the first (matrix, grid) tuple has corresponding sizes
    M, G = MGs[1]
    m, n = size(M)
    @boundscheck m == n || throw(BoundsError)
    @boundscheck m == 2*(4+G.nlat_half) || throw(BoundsError)

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
    nside = 4+G.nlat_half       # side length of a basepixel matrix

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
            r = 4+min(j, nlat_half) - iq             # row in matrix m (1-based)
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
    return Tuple(Mi for (Mi, Gi) in MGs)             # 
end