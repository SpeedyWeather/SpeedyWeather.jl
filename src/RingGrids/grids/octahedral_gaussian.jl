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
full_grid(::Type{<:OctahedralGaussianArray}) = FullGaussianArray

# SIZE
nlat_odd(::Type{<:OctahedralGaussianArray}) = false
npoints_pole(::Type{<:OctahedralGaussianArray}) = 16
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
get_colat(::Type{<:OctahedralGaussianArray}, nlat_half::Integer) = get_colat(FullGaussianArray, nlat_half)
function get_lon_per_ring(Grid::Type{<:OctahedralGaussianArray}, nlat_half::Integer, j::Integer)
    nlon = get_nlon_per_ring(Grid, nlat_half, j)
    return collect(0:2π/nlon:2π-π/nlon)
end

## QUADRATURE
get_quadrature_weights(::Type{<:OctahedralGaussianArray}, nlat_half::Integer) = gaussian_weights(nlat_half)

## INDEXING
function each_index_in_ring(Grid::Type{<:OctahedralGaussianArray},
                            j::Integer,                     # ring index north to south
                            nlat_half::Integer)             # resolution param

    nlat = get_nlat(Grid, nlat_half)
    
    # TODO make m, o dependent
    m, o = npoints_added_per_ring(OctahedralGaussianArray), npoints_pole(OctahedralGaussianArray)
    m != 4 || o != 16 && @warn "This algorithm has not been generalised for m!=4, o!=16."

    @boundscheck 0 < j <= nlat || throw(BoundsError)        # ring index valid?
    if j <= nlat_half                                       # northern hemisphere incl Equator
        index_1st = 2j*(j+7) - 15                           # first in-ring index i
        index_end = 2j*(j+9)                                # last in-ring index i
    else                                                    # southern hemisphere excl Equator
        j = nlat - j + 1                                    # mirror ring index around Equator
        n = get_npoints2D(Grid, nlat_half) + 1              # number of grid points + 1
        index_1st = n - 2j*(j+9)                            # count backwards
        index_end = n - (2j*(j+7) - 15)
    end
    return index_1st:index_end                              # range of i's in ring
end

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

# """
#     G = OctahedralGaussianGrid{T}

# An Octahedral Gaussian grid that uses `nlat` Gaussian latitudes, but a decreasing number of longitude
# points per latitude ring towards the poles. Starting with 20 equi-spaced longitude points (starting at 0˚E)
# on the rings around the poles, each latitude ring towards the equator has consecuitively 4 more points,
# one for each face of the octahedron. E.g. 20, 24, 28, 32, ...nlon-4, nlon, nlon, nlon-4, ..., 32, 28, 24, 20.
# The maximum number of longitue points is `nlon`. The values of all grid points are stored in a vector
# field `v` that unravels the data 0 to 360˚, then ring by ring, which are sorted north to south.
# """
# # OctahedralGaussianGrid