"""A `FullOctaHEALPixGrid` is like a `OctaHEALPixGrid` but with every latitude ring having the same number of longitude
points (a full grid). This grid is mostly defined for output to minimize the interpolation needed from an
`OctaHEALPixGrid` to a full grid required for output. A `FullOctaHEALPixGrid` has none of the equal-area properties of the `OctaHEALPixGrid`.
It only shares the latitudes with the `OctaHEALPixGrid` but uses the longitudes from the `FullGaussianGrid`
without offset, i.e. the first longitude point on every ring is at 0˚E.
$(TYPEDFIELDS)"""
struct FullOctaHEALPixGrid{A, V, W} <: AbstractFullGrid{A}
    nlat_half::Int      # number of latitudes on one hemisphere
    architecture::A     # information about device, CPU/GPU
    rings::V            # precomputed ring indices
    whichring::W        # precomputed ring index for each grid point ij
end

nonparametric_type(::Type{<:FullOctaHEALPixGrid}) = FullOctaHEALPixGrid

# FIELD
const FullOctaHEALPixField{T, N} = Field{T, N, ArrayType, Grid} where {ArrayType, Grid<:FullOctaHEALPixGrid}

# define grid_type (i) without T, N, (ii) with T, (iii) with T, N but not with <:?Field
# to not have precendence over grid_type(::Type{Field{...})
grid_type(::Type{FullOctaHEALPixField}) = FullOctaHEALPixGrid
grid_type(::Type{FullOctaHEALPixField{T}}) where T = FullOctaHEALPixGrid
grid_type(::Type{FullOctaHEALPixField{T, N}}) where {T, N} = FullOctaHEALPixGrid

function Base.showarg(io::IO, F::Field{T, N, ArrayType, Grid}, toplevel) where {T, N, ArrayType, Grid<:FullOctaHEALPixGrid{A}} where A <: AbstractArchitecture
    print(io, "FullOctaHEALPixField{$T, $N}")
    toplevel && print(io, " on ", nonparametric_type(ArrayType))
    toplevel && print(io, " on ", A)
end

# SIZE
nlat_odd(::Type{<:FullOctaHEALPixGrid}) = true
get_npoints(::Type{<:FullOctaHEALPixGrid}, nlat_half::Integer) = 4nlat_half * (2nlat_half-1)
get_nlat_half(::Type{<:FullOctaHEALPixGrid}, npoints::Integer) = round(Int, 1/4 + sqrt(1/16 + npoints/8))
get_nlon(::Type{<:FullOctaHEALPixGrid}, nlat_half::Integer) = 4nlat_half

## COORDINATES
get_latd(::Type{<:FullOctaHEALPixGrid}, nlat_half::Integer) = get_latd(OctaHEALPixGrid, nlat_half)
get_lond(::Type{<:FullOctaHEALPixGrid}, nlat_half::Integer) = get_lond(FullGaussianGrid, nlat_half)

# QUADRATURE (use weights from reduced grids though!)
get_quadrature_weights(::Type{<:FullOctaHEALPixGrid}, nlat_half::Integer) =
    equal_area_weights(OctaHEALPixGrid, nlat_half)

Adapt.@adapt_structure FullOctaHEALPixGrid
