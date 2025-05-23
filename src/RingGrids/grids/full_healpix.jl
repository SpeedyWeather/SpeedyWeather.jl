"""A `FullHEALPixArray` is an array of full grids, subtyping `AbstractFullGridArray`, that use
HEALPix latitudes for each ring. This type primarily equists to interpolate data from
the (reduced) HEALPixGrid onto a full grid for output.

First dimension of the underlying `N`-dimensional `data` represents the horizontal dimension,
in ring order (0 to 360˚E, then north to south), other dimensions are used for the vertical
and/or time or other dimensions. The resolution parameter of the horizontal grid is
`nlat_half` (number of latitude rings on one hemisphere, Equator included) and the ring indices
are precomputed in `rings`. Fields are
$(TYPEDFIELDS)"""
struct FullHEALPixGrid{A, V} <: AbstractFullGrid{A}
    nlat_half::Int      # number of latitudes on one hemisphere
    architecture::A     # information about device, CPU/GPU
    rings::V            # precomputed ring indices
end

nonparametric_type(::Type{<:FullHEALPixGrid}) = FullHEALPixGrid

# FIELD
const FullHEALPixField{T, N} = Field{T, N, Architecture, Grid} where {Architecture, Grid<:FullHEALPixGrid}

# define grid_type (i) without T, N, (ii) with T, (iii) with T, N but not with <:?Field
# to not have precendence over grid_type(::Type{Field{...})
grid_type(::Type{FullHEALPixField}) = FullHEALPixGrid
grid_type(::Type{FullHEALPixField{T}}) where T = FullHEALPixGrid
grid_type(::Type{FullHEALPixField{T, N}}) where {T, N} = FullHEALPixGrid

Base.show(io::IO, F::Type{<:FullHEALPixField{T, N}}) where {T, N} = print(io, "FullHEALPixField{$T, $N}")

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