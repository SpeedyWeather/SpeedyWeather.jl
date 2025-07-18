"""An `OctahedralGaussianGrid` is a discretization of the sphere that uses Gaussian latitudes for each latitude ring.
As a reduced grid it has a different number of longitude points on every latitude ring. 
These grids are called octahedral because after starting with 20 points on the first ring around the north pole
they increase the number of longitude points for each ring by 4, such that they can be conceptually thought of
as lying on the 4 faces of an octahedron on each hemisphere. Hence, these grids have 20, 24, 28, ... longitude points
for ring 1, 2, 3, ... There is no ring on the Equator and the two rings around it have 16 + 4nlat_half longitude points
before reducing the number of longitude points per ring by 4 towards the southern-most ring j = nlat.
The first point on every ring is at longitude 0˚E, i.e. no offset is applied to the longitude values (in contrast to HEALPix grids).

The first dimension of data on this grid (a `Field`) represents the horizontal dimension,
in ring order (0 to 360˚E, then north to south), other dimensions can be used for the vertical and/or
time or other dimensions.  Note that a `Grid` does not contain any data, it only describes
the discretization of the space, see `Field` for a data on a `Grid`.  But a "grid" only defines the
two horizontal dimensions, two fields, one 2D and one 3D, possibly different ArrayTypes or element types,
can share the same grid which just defines the discretization and the architecture (CPU/GPU) the grid is on.

The resolution parameter of the horizontal grid is `nlat_half` (number of latitude rings on one hemisphere,
Equator included) `rings` are the precomputed ring indices, the the example above `rings = [1:20, 21:44, 45:72, ...]`. 
`whichring` is a precomputed vector of ring indices for each grid point ij, i.e. `whichring[ij]` gives the ring index j of grid point ij.
For efficient looping see `eachring` and `eachgrid`.
Fields are
$(TYPEDFIELDS)"""
struct OctahedralGaussianGrid{A, V, W} <: AbstractReducedGrid{A}
    nlat_half::Int                  # number of latitudes on one hemisphere
    architecture::A                 # information about device, CPU/GPU
    rings::V                        # precomputed ring indices (ring j -> grid point ij)
    whichring::W                    # precomputed ring index for each grid point ij
end

# TYPES
nonparametric_type(::Type{<:OctahedralGaussianGrid}) = OctahedralGaussianGrid
full_grid_type(::Type{<:OctahedralGaussianGrid}) = FullGaussianGrid

# FIELD
const OctahedralGaussianField{T, N} = Field{T, N, ArrayType, Grid} where {ArrayType, Grid<:OctahedralGaussianGrid}

# define grid_type (i) without T, N, (ii) with T, (iii) with T, N but not with <:?Field
# to not have precendence over grid_type(::Type{Field{...})
grid_type(::Type{OctahedralGaussianField}) = OctahedralGaussianGrid
grid_type(::Type{OctahedralGaussianField{T}}) where T = OctahedralGaussianGrid
grid_type(::Type{OctahedralGaussianField{T, N}}) where {T, N} = OctahedralGaussianGrid

function Base.showarg(io::IO, F::Field{T, N, ArrayType, Grid}, toplevel) where {T, N, ArrayType, Grid<:OctahedralGaussianGrid{A}} where A <: AbstractArchitecture
    print(io, "OctahedralGaussianField{$T, $N}")
    toplevel && print(io, " on ", nonparametric_type(ArrayType))
    toplevel && print(io, " on ", A)
end 

# SIZE
nlat_odd(::Type{<:OctahedralGaussianGrid}) = false
npoints_pole(::Type{<:OctahedralGaussianGrid}) = 16
npoints_added_per_ring(::Type{<:OctahedralGaussianGrid}) = 4

function get_npoints(::Type{<:OctahedralGaussianGrid}, nlat_half::Integer)
    m, o = npoints_added_per_ring(OctahedralGaussianGrid), npoints_pole(OctahedralGaussianGrid)
    return m*nlat_half^2 + (2o+m)*nlat_half
end

function get_nlat_half(::Type{<:OctahedralGaussianGrid}, npoints2D::Integer)
    m, o = npoints_added_per_ring(OctahedralGaussianGrid), npoints_pole(OctahedralGaussianGrid)
    return round(Int, -(2o + m)/2m + sqrt(((2o+m)/2m)^2 + npoints2D/m))
end

function get_nlon_per_ring(Grid::Type{<:OctahedralGaussianGrid}, nlat_half::Integer, j::Integer)
    nlat = get_nlat(Grid, nlat_half)
    @assert 0 < j <= nlat "Ring $j is outside O$nlat_half grid."
    m, o = npoints_added_per_ring(OctahedralGaussianGrid), npoints_pole(OctahedralGaussianGrid)
    j = j > nlat_half ? nlat - j + 1 : j      # flip north south due to symmetry
    return o + m*j
end

# maybe define at some point for Matrix(::OctahedralGaussianGrid)
# matrix_size(G::OctahedralGaussianGrid) = (2*(4+G.nlat_half), 2*(4+G.nlat_half+1))

## COORDINATES
get_latd(::Type{<:OctahedralGaussianGrid}, nlat_half::Integer) = get_latd(FullGaussianGrid, nlat_half)
function get_lond_per_ring(Grid::Type{<:OctahedralGaussianGrid}, nlat_half::Integer, j::Integer)
    nlon = get_nlon_per_ring(Grid, nlat_half, j)
    return collect(0:360/nlon:360-180/nlon)
end

## QUADRATURE
get_quadrature_weights(::Type{<:OctahedralGaussianGrid}, nlat_half::Integer) = gaussian_weights(nlat_half)

## INDEXING
"""$(TYPEDSIGNATURES) precompute a `Vector{UnitRange{Int}} to index grid points on
every ring `j` (elements of the vector) of `Grid` at resolution `nlat_half`.
See `eachring` and `eachgrid` for efficient looping over grid points."""
function each_index_in_ring!(   rings,
                                Grid::Type{<:OctahedralGaussianGrid},
                                nlat_half::Integer) # resolution param

    nlat = length(rings)
    @boundscheck nlat == get_nlat(Grid, nlat_half) || throw(BoundsError)
    m, o = npoints_added_per_ring(OctahedralGaussianGrid), npoints_pole(OctahedralGaussianGrid)

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

Adapt.@adapt_structure OctahedralGaussianGrid