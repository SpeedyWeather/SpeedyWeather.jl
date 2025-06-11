"""A `FullOctaHEALPixArray` is an array of full grids, subtyping `AbstractFullGridArray` that use
OctaHEALPix latitudes for each ring. This type primarily equists to interpolate data from
the (reduced) OctaHEALPixGrid onto a full grid for output.

First dimension of the underlying `N`-dimensional `data` represents the horizontal dimension,
in ring order (0 to 360˚E, then north to south), other dimensions are used for the vertical
and/or time or other dimensions. The resolution parameter of the horizontal grid is
`nlat_half` (number of latitude rings on one hemisphere, Equator included) and the ring indices
are precomputed in `rings`. Fields are
$(TYPEDFIELDS)"""
struct FullOctaHEALPixGrid{A, V, W} <: AbstractFullGrid{A}
    nlat_half::Int      # number of latitudes on one hemisphere
    architecture::A     # information about device, CPU/GPU
    rings::V            # precomputed ring indices
    whichring::W        # precomputed ring index for each grid point ij
end

nonparametric_type(::Type{<:FullOctaHEALPixGrid}) = FullOctaHEALPixGrid

# FIELD
const FullOctaHEALPixField{T, N} = Field{T, N, Architecture, Grid} where {Architecture, Grid<:FullOctaHEALPixGrid}

# define grid_type (i) without T, N, (ii) with T, (iii) with T, N but not with <:?Field
# to not have precendence over grid_type(::Type{Field{...})
grid_type(::Type{FullOctaHEALPixField}) = FullOctaHEALPixGrid
grid_type(::Type{FullOctaHEALPixField{T}}) where T = FullOctaHEALPixGrid
grid_type(::Type{FullOctaHEALPixField{T, N}}) where {T, N} = FullOctaHEALPixGrid

Base.show(io::IO, F::Type{<:FullOctaHEALPixField{T, N}}) where {T, N} =
    print(io, "FullOctaHEALPixField{$T, $N}")

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