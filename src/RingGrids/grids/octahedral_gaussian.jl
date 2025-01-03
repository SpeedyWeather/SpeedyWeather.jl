"""An `OctahedralGaussianArray` is an array of octahedral grids, subtyping `AbstractReducedGridArray`,
that use Gaussian latitudes for each ring. First dimension of the underlying `N`-dimensional `data`
represents the horizontal dimension, in ring order (0 to 360˚E, then north to south),
other dimensions are used for the vertical and/or time or other dimensions.
The resolution parameter of the horizontal grid is `nlat_half` (number of latitude rings
on one hemisphere, Equator included) and the ring indices are precomputed in `rings`.

These grids are called octahedral because after starting with 20 points on the first ring around the
north pole (default) they increase the number of longitude points for each ring by 4, such that
they can be conceptually thought of as lying on the 4 faces of an octahedron on each hemisphere.
Hence, these grids have 20, 24, 28, ... longitude points for ring 1, 2, 3, ... There is no
ring on the Equator and the two rings around it have 16 + 4nlat_half longitude points before
reducing the number of longitude points per ring by 4 towards the southern-most ring
j = nlat. `rings` are the precomputed ring indices, the the example above
`rings = [1:20, 21:44, 45:72, ...]`. For efficient looping see `eachring` and `eachgrid`.
Fields are
$(TYPEDFIELDS)"""
struct OctahedralGaussianArray{T, N, ArrayType <: AbstractArray{T, N}} <: AbstractReducedGridArray{T, N, ArrayType}
    data::ArrayType                 # data array, ring by ring, north to south
    nlat_half::Int                  # number of latitudes on one hemisphere
    rings::Vector{UnitRange{Int}}   # TODO make same array type as data?

    OctahedralGaussianArray(data::A, nlat_half, rings) where {A <: AbstractArray{T, N}} where {T, N} =
        check_inputs(data, nlat_half, rings, OctahedralGaussianArray) ?
        new{T, N, A}(data, nlat_half, rings) :
        error_message(data, nlat_half, rings, OctahedralGaussianArray, T, N, A)
end

# TYPES
const OctahedralGaussianGrid{T} = OctahedralGaussianArray{T, 1, Vector{T}}
nonparametric_type(::Type{<:OctahedralGaussianArray}) = OctahedralGaussianArray
horizontal_grid_type(::Type{<:OctahedralGaussianArray}) = OctahedralGaussianGrid
full_array_type(::Type{<:OctahedralGaussianArray}) = FullGaussianArray

"""An `OctahedralGaussianArray` but constrained to `N=1` dimensions (horizontal only) and data is a `Vector{T}`."""
OctahedralGaussianGrid

# SIZE
nlat_odd(::Type{<:OctahedralGaussianArray}) = false

"""$(TYPEDSIGNATURES) [EXPERIMENTAL] additional number of longitude points on the first and last ring.
Change to 0 to start with 4 points on the first ring."""
npoints_pole(::Type{<:OctahedralGaussianArray}) = 16

"""$(TYPEDSIGNATURES) [EVEN MORE EXPERIMENTAL] number of longitude points added (removed) for every ring
towards the Equator (on the southern hemisphere towards the south pole)."""
npoints_added_per_ring(::Type{<:OctahedralGaussianArray}) = 4

function get_npoints2D(::Type{<:OctahedralGaussianArray}, nlat_half::Integer)
    m, o = npoints_added_per_ring(OctahedralGaussianArray), npoints_pole(OctahedralGaussianArray)
    return m*nlat_half^2 + (2o+m)*nlat_half
end

function get_nlat_half(::Type{<:OctahedralGaussianArray}, npoints2D::Integer)
    m, o = npoints_added_per_ring(OctahedralGaussianArray), npoints_pole(OctahedralGaussianArray)
    return round(Int, -(2o + m)/2m + sqrt(((2o+m)/2m)^2 + npoints2D/m))
end

function get_nlon_per_ring(Grid::Type{<:OctahedralGaussianArray}, nlat_half::Integer, j::Integer)
    nlat = get_nlat(Grid, nlat_half)
    @assert 0 < j <= nlat "Ring $j is outside O$nlat_half grid."
    m, o = npoints_added_per_ring(OctahedralGaussianArray), npoints_pole(OctahedralGaussianArray)
    j = j > nlat_half ? nlat - j + 1 : j      # flip north south due to symmetry
    return o + m*j
end

# maybe define at some point for Matrix(::OctahedralGaussianGrid)
# matrix_size(G::OctahedralGaussianGrid) = (2*(4+G.nlat_half), 2*(4+G.nlat_half+1))

## COORDINATES
get_latd(::Type{<:OctahedralGaussianArray}, nlat_half::Integer) = get_latd(FullGaussianArray, nlat_half)
function get_lond_per_ring(Grid::Type{<:OctahedralGaussianArray}, nlat_half::Integer, j::Integer)
    nlon = get_nlon_per_ring(Grid, nlat_half, j)
    return collect(0:360/nlon:360-180/nlon)
end

## QUADRATURE
get_quadrature_weights(::Type{<:OctahedralGaussianArray}, nlat_half::Integer) = gaussian_weights(nlat_half)

## INDEXING
"""$(TYPEDSIGNATURES) precompute a `Vector{UnitRange{Int}} to index grid points on
every ring `j` (elements of the vector) of `Grid` at resolution `nlat_half`.
See `eachring` and `eachgrid` for efficient looping over grid points."""
function each_index_in_ring!(   rings::Vector{<:UnitRange{<:Integer}},
                                Grid::Type{<:OctahedralGaussianArray},
                                nlat_half::Integer) # resolution param

    nlat = length(rings)
    @boundscheck nlat == get_nlat(Grid, nlat_half) || throw(BoundsError)
    m, o = npoints_added_per_ring(OctahedralGaussianArray), npoints_pole(OctahedralGaussianArray)

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