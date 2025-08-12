"""A `FullClenshawGrid` is a discretization of the sphere, a full grid, subtyping `AbstractFullGrid`,
using equidistant latitudes for each ring (a regular lon-lat grid). These require the
Clenshaw-Curtis quadrature in the spectral transform, hence the name. One ring is on the equator,
total number of rings is odd, no rings on the north or south pole.

The first dimension of data on this grid (a `Field`) represents the horizontal dimension,
in ring order (0 to 360˚E, then north to south), other dimensions can be used for the vertical and/or
time or other dimensions.  Note that a `Grid` does not contain any data, it only describes
the discretization of the space, see `Field` for a data on a `Grid`.  But a "grid" only defines the
two horizontal dimensions, two fields, one 2D and one 3D, possibly different ArrayTypes or element types,
can share the same grid which just defines the discretization and the architecture (CPU/GPU) the grid is on.

The resolution parameter of the horizontal grid is `nlat_half` (number of latitude rings on one hemisphere,
Equator included) and the ring indices are precomputed in `rings`.
$(TYPEDFIELDS)"""
struct FullClenshawGrid{A, V, W} <: AbstractFullGrid{A}
    nlat_half::Int      # number of latitudes on one hemisphere
    architecture::A     # information about device, CPU/GPU
    rings::V            # precomputed ring indices
    whichring::W        # precomputed ring index for each grid point ij
end

# TYPES
Architectures.nonparametric_type(::Type{<:FullClenshawGrid}) = FullClenshawGrid

# FIELD
const FullClenshawField{T, N} = Field{T, N, ArrayType, Grid} where {ArrayType, Grid<:FullClenshawGrid}

# define grid_type (i) without T, N, (ii) with T, (iii) with T, N but not with <:?Field
# to not have precendence over grid_type(::Type{Field{...})
grid_type(::Type{FullClenshawField}) = FullClenshawGrid
grid_type(::Type{FullClenshawField{T}}) where T = FullClenshawGrid
grid_type(::Type{FullClenshawField{T, N}}) where {T, N} = FullClenshawGrid

function Base.showarg(io::IO, F::Field{T, N, ArrayType, Grid}, toplevel) where {T, N, ArrayType, Grid<:FullClenshawGrid{A}} where A <: AbstractArchitecture
    print(io, "FullClenshawField{$T, $N}")
    toplevel && print(io, " on ", nonparametric_type(ArrayType))
    toplevel && print(io, " on ", A)
end

# SIZE
nlat_odd(::Type{<:FullClenshawGrid}) = true
get_npoints(::Type{<:FullClenshawGrid}, nlat_half::Integer) = 8 * nlat_half^2 - 4nlat_half
get_nlat_half(::Type{<:FullClenshawGrid}, npoints::Integer) = round(Int, 1/4 + sqrt(1/16 + npoints/8))
get_nlon(::Type{<:FullClenshawGrid}, nlat_half::Integer) = 4nlat_half

## COORDINATES
get_latd(::Type{<:FullClenshawGrid}, nlat_half::Integer) = [90 - 90j/nlat_half for j in 1:2nlat_half-1]
get_lond(::Type{<:FullClenshawGrid}, nlat_half::Integer) = get_lond(FullGaussianGrid, nlat_half)

# QUADRATURE
get_quadrature_weights(::Type{<:FullClenshawGrid}, nlat_half::Integer) = clenshaw_curtis_weights(nlat_half)

Adapt.@adapt_structure FullClenshawGrid