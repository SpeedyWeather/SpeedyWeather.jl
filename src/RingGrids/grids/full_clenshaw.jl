"""A `FullClenshawGrid` is a discretization of the sphere, a full grid, subtyping `AbstractFullGrid`,
using equidistant latitudes for each ring (a regular lon-lat grid). These require the
Clenshaw-Curtis quadrature in the spectral transform, hence the name. One ring is on the equator,
total number of rings is odd, no rings on the north or south pole. First dimension of the
underlying `N`-dimensional `data` represents the horizontal dimension, in ring order
(0 to 360˚E, then north to south), other dimensions are used for the vertical and/or
time or other dimensions. The resolution parameter of the horizontal grid is `nlat_half`
(number of latitude rings on one hemisphere, Equator included) and the ring indices are
precomputed in `rings`. Note that a `Grid` does not contain any data, it only describes
the discretization of the space, see `Field` for a data on a `Grid`.
$(TYPEDFIELDS)"""
struct FullClenshawGrid{A, V} <: AbstractFullGrid
    nlat_half::Int      # number of latitudes on one hemisphere
    architecture::A     # information about device, CPU/GPU
    rings::V            # precomputed ring indices
end

# TYPES
nonparametric_type(::Type{<:FullClenshawGrid}) = FullClenshawGrid

const FullClenshawField{T, N} = Field{T, N, A, G} where {A, G<:FullClenshawGrid}
Base.show(io::IO, F::Type{<:FullClenshawField{T, N}}) where {T, N} = print(io, "FullClenshawField{$T, $N}")

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