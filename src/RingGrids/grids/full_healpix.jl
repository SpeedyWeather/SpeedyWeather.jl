"""A `FullHEALPixGrid` is like a `HEALPixGrid` but with every latitude ring having the same number of longitude
points (a full grid). This grid is mostly defined for output to minimize the interpolation needed from a
HEALPixGrid to a full grid. A `FullHEALPixGrid` has none of the equal-area properties of the `HEALPixGrid`.
It only shares the latitudes with the `HEALPixGrid` but uses the longitudes from the `FullGaussianGrid`
without offset, i.e. the first longitude point on every ring is at 0˚E.
$(TYPEDFIELDS)"""
struct FullHEALPixGrid{A, V, W} <: AbstractFullGrid{A}
    nlat_half::Int      # number of latitudes on one hemisphere
    architecture::A     # information about device, CPU/GPU
    rings::V            # precomputed ring indices
    whichring::W        # precomputed ring index for each grid point ij
end

Architectures.nonparametric_type(::Type{<:FullHEALPixGrid}) = FullHEALPixGrid

# FIELD
const FullHEALPixField{T, N} = Field{T, N, ArrayType, Grid} where {ArrayType, Grid<:FullHEALPixGrid}

# define grid_type (i) without T, N, (ii) with T, (iii) with T, N but not with <:?Field
# to not have precendence over grid_type(::Type{Field{...})
grid_type(::Type{FullHEALPixField}) = FullHEALPixGrid
grid_type(::Type{FullHEALPixField{T}}) where T = FullHEALPixGrid
grid_type(::Type{FullHEALPixField{T, N}}) where {T, N} = FullHEALPixGrid

function Base.showarg(io::IO, F::Field{T, N, ArrayType, Grid}, toplevel) where {T, N, ArrayType, Grid<:FullHEALPixGrid{A}} where A <: AbstractArchitecture
    print(io, "FullHEALPixField{$T, $N}")
    toplevel && print(io, " as ", nonparametric_type(ArrayType))
    toplevel && print(io, " on ", F.grid.architecture)
end

# SIZE
nlat_odd(::Type{<:FullHEALPixGrid}) = true
get_npoints(::Type{<:FullHEALPixGrid}, nlat_half::Integer) = 4nlat_half * (2nlat_half-1)
get_nlat_half(::Type{<:FullHEALPixGrid}, npoints::Integer) = round(Int, 1/4 + sqrt(1/16 + npoints/8))
get_nlon(::Type{<:FullHEALPixGrid}, nlat_half::Integer) = 4nlat_half

## COORDINATES
get_latd(::Type{<:FullHEALPixGrid}, nlat_half::Integer) = get_latd(HEALPixGrid, nlat_half)
get_lond(::Type{<:FullHEALPixGrid}, nlat_half::Integer) = get_lond(FullGaussianGrid, nlat_half)

# QUADRATURE (use weights from reduced grids though!)
get_quadrature_weights(::Type{<:FullHEALPixGrid}, nlat_half::Integer) =
    equal_area_weights(HEALPixGrid, nlat_half)

Adapt.@adapt_structure FullHEALPixGrid