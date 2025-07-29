"""An `OctaminimalGaussianGrid` is a discretization of the sphere that uses Gaussian latitudes for each latitude ring.
It is very similar to the `OctahedralGaussianGrid`, with a few differences:

1. It starts with 4 longitude points around the poles, so 4, 8, 12, 16, ... on ring 1, 2, 3, 4, ... reducing the number
of total number of grid points over the `OctahedralGaussianGrid`, hence the name "octaminimal" (octahedral minimal).
2. The first longitude points on every ring have an offset, such that the cell face between the first and last longitude point
   on every ring is at 0˚E (like the HEALPix grid)

Fields are
$(TYPEDFIELDS)"""
struct OctaminimalGaussianGrid{A, V, W} <: AbstractReducedGrid{A}
    nlat_half::Int                  # number of latitudes on one hemisphere
    architecture::A                 # information about device, CPU/GPU
    rings::V                        # precomputed ring indices
    whichring::W                    # precomputed ring index for each grid point ij
end

# TYPES
Architectures.nonparametric_type(::Type{<:OctaminimalGaussianGrid}) = OctaminimalGaussianGrid
full_grid_type(::Type{<:OctaminimalGaussianGrid}) = FullGaussianGrid

# FIELD
const OctaminimalGaussianField{T, N} = Field{T, N, ArrayType, Grid} where {ArrayType, Grid<:OctaminimalGaussianGrid}

# define grid_type (i) without T, N, (ii) with T, (iii) with T, N but not with <:?Field
# to not have precendence over grid_type(::Type{Field{...})
grid_type(::Type{OctaminimalGaussianField}) = OctaminimalGaussianGrid
grid_type(::Type{OctaminimalGaussianField{T}}) where T = OctaminimalGaussianGrid
grid_type(::Type{OctaminimalGaussianField{T, N}}) where {T, N} = OctaminimalGaussianGrid

function Base.showarg(io::IO, F::Field{T, N, ArrayType, Grid}, toplevel) where {T, N, ArrayType, Grid<:OctaminimalGaussianGrid{A}} where A <: AbstractArchitecture
    print(io, "OctaminimalGaussianField{$T, $N}")
    toplevel && print(io, " on ", nonparametric_type(ArrayType))
    toplevel && print(io, " on ", A)
end

# SIZE
nlat_odd(::Type{<:OctaminimalGaussianGrid}) = false
npoints_pole(::Type{<:OctaminimalGaussianGrid}) = 0
npoints_added_per_ring(::Type{<:OctaminimalGaussianGrid}) = 4

function get_npoints(::Type{<:OctaminimalGaussianGrid}, nlat_half::Integer)
    m, o = npoints_added_per_ring(OctaminimalGaussianGrid), npoints_pole(OctaminimalGaussianGrid)
    return m*nlat_half^2 + (2o+m)*nlat_half
end

function get_nlat_half(::Type{<:OctaminimalGaussianGrid}, npoints::Integer)
    m, o = npoints_added_per_ring(OctaminimalGaussianGrid), npoints_pole(OctaminimalGaussianGrid)
    return round(Int, -(2o + m)/2m + sqrt(((2o+m)/2m)^2 + npoints/m))
end

function get_nlon_per_ring(Grid::Type{<:OctaminimalGaussianGrid}, nlat_half::Integer, j::Integer)
    nlat = get_nlat(Grid, nlat_half)
    @assert 0 < j <= nlat "Ring $j is outside O$nlat_half grid."
    m, o = npoints_added_per_ring(OctaminimalGaussianGrid), npoints_pole(OctaminimalGaussianGrid)
    j = j > nlat_half ? nlat - j + 1 : j      # flip north south due to symmetry
    return o + m*j
end

# maybe define at some point for Matrix(::OctaminimalGaussianGrid)
# matrix_size(G::OctaminimalGaussianGrid) = (2*(4+G.nlat_half), 2*(4+G.nlat_half+1))

## COORDINATES
get_latd(::Type{<:OctaminimalGaussianGrid}, nlat_half::Integer) = get_latd(FullGaussianGrid, nlat_half)
function get_lond_per_ring(Grid::Type{<:OctaminimalGaussianGrid}, nlat_half::Integer, j::Integer)
    nlon = get_nlon_per_ring(Grid, nlat_half, j)
    return collect(180/nlon:360/nlon:360)       # use HEALPix definition for longitudes
end

## QUADRATURE
get_quadrature_weights(::Type{<:OctaminimalGaussianGrid}, nlat_half::Integer) = gaussian_weights(nlat_half)

## INDEXING
"""$(TYPEDSIGNATURES) precompute a `Vector{UnitRange{Int}} to index grid points on
every ring `j` (elements of the vector) of `Grid` at resolution `nlat_half`.
See `eachring` and `eachgrid` for efficient looping over grid points."""
function each_index_in_ring!(   rings,
                                Grid::Type{<:OctaminimalGaussianGrid},
                                nlat_half::Integer) # resolution param

    nlat = length(rings)
    @boundscheck nlat == get_nlat(Grid, nlat_half) || throw(BoundsError)
    m, o = npoints_added_per_ring(OctaminimalGaussianGrid), npoints_pole(OctaminimalGaussianGrid)

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

Adapt.@adapt_structure OctaminimalGaussianGrid