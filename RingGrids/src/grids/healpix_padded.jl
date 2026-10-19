"""A `HEALPixPaddedGrid` shares the latitudes of the `HEALPixGrid` but refines the rings
with fewer than `npoints_min = 20` longitudes near the poles so that the HEALPix points stay
a subset of the padded grid points (they *nest*). Ring `j` has `k·n` longitudes with `n`
the number of HEALPix longitudes on that ring and `k = cld(20, n)` the smallest integer
refinement giving at least 20 points: 20, 24, 24, 32, 20, 24, 28, ... from the pole,
identical to HEALPix from ring 5 on. (20 points on rings 2-4 cannot nest as 20 is not a
multiple of 8, 12, 16.) The few extra points around the poles allow these rings to resolve
the zonal wavenumbers that fields carry there, see Courtier and Naughton, 1994,
DOI:10.1002/qj.49712051913.

The longitudes follow the HEALPix convention with or without an offset of half a grid spacing
from 0˚E. A refined ring has an offset if the HEALPix ring has one and `k` is odd, such that
HEALPix point `i` on ring `j` is padded point `k(i-1) + first` with `first` from
`healpix_nesting(Grid, nlat_half, j)`. Use `copyto!(::HEALPixField, ::HEALPixPaddedField)`
to subsample and `copyto!(::HEALPixPaddedField, ::HEALPixField)` to linearly interpolate
in longitude. Requires `nlat_half >= 16` so that only the 4 rings around each pole are
padded. The grid uses the equal-area ring weights of HEALPix, but is not equal-area
as the refined pixels are smaller. Fields are
$(TYPEDFIELDS)"""
struct HEALPixPaddedGrid{A, V, W, IntType} <: AbstractReducedGrid{A}
    nlat_half::IntType              # number of latitudes on one hemisphere
    architecture::A                 # information about device, CPU/GPU
    rings::V                        # precomputed ring indices
    whichring::W                    # precomputed ring index for each grid point ij
end

Architectures.nonparametric_type(::Type{<:HEALPixPaddedGrid}) = HEALPixPaddedGrid
full_grid_type(::Type{<:HEALPixPaddedGrid}) = FullHEALPixGrid

# FIELD
const HEALPixPaddedField{T, N} = Field{T, N, ArrayType, Grid} where {ArrayType, Grid <: HEALPixPaddedGrid}

# define grid_type (i) without T, N, (ii) with T, (iii) with T, N but not with <:?Field
# to not have precendence over grid_type(::Type{Field{...})
grid_type(::Type{HEALPixPaddedField}) = HEALPixPaddedGrid
grid_type(::Type{HEALPixPaddedField{T}}) where {T} = HEALPixPaddedGrid
grid_type(::Type{HEALPixPaddedField{T, N}}) where {T, N} = HEALPixPaddedGrid

function Base.showarg(io::IO, F::Field{T, N, ArrayType, Grid}, toplevel) where {T, N, ArrayType, Grid <: HEALPixPaddedGrid{A}} where {A <: AbstractArchitecture}
    print(io, "HEALPixPaddedField{$T, $N}")
    toplevel && print(io, " as ", nonparametric_type(ArrayType))
    return toplevel && print(io, " on ", F.grid.architecture)
end


"""An `OctaHEALPixPaddedGrid` is to the `OctaHEALPixGrid` what the `HEALPixPaddedGrid` is to
the `HEALPixGrid`: same latitudes but rings with fewer than 20 longitudes are refined by the
smallest integer factor `k = cld(20, n)` so that the `OctaHEALPixGrid` points nest, i.e. rings
have 20, 24, 24, 32, 20, 24, 28, ... points from the pole. `copyto!` between both grids
subsamples or linearly interpolates in longitude. Requires `nlat_half >= 8`. Fields are
$(TYPEDFIELDS)"""
struct OctaHEALPixPaddedGrid{A, V, W, IntType} <: AbstractReducedGrid{A}
    nlat_half::IntType              # number of latitudes on one hemisphere
    architecture::A                 # information about device, CPU/GPU
    rings::V                        # precomputed ring indices
    whichring::W                    # precomputed ring index for each grid point ij
end

Architectures.nonparametric_type(::Type{<:OctaHEALPixPaddedGrid}) = OctaHEALPixPaddedGrid
full_grid_type(::Type{<:OctaHEALPixPaddedGrid}) = FullOctaHEALPixGrid

# FIELD
const OctaHEALPixPaddedField{T, N} = Field{T, N, ArrayType, Grid} where {ArrayType, Grid <: OctaHEALPixPaddedGrid}

grid_type(::Type{OctaHEALPixPaddedField}) = OctaHEALPixPaddedGrid
grid_type(::Type{OctaHEALPixPaddedField{T}}) where {T} = OctaHEALPixPaddedGrid
grid_type(::Type{OctaHEALPixPaddedField{T, N}}) where {T, N} = OctaHEALPixPaddedGrid

function Base.showarg(io::IO, F::Field{T, N, ArrayType, Grid}, toplevel) where {T, N, ArrayType, Grid <: OctaHEALPixPaddedGrid{A}} where {A <: AbstractArchitecture}
    print(io, "OctaHEALPixPaddedField{$T, $N}")
    toplevel && print(io, " as ", nonparametric_type(ArrayType))
    return toplevel && print(io, " on ", F.grid.architecture)
end

# shared implementation for both padded grids, relative to their unpadded grid
const PaddedHEALPixGrids = Union{HEALPixPaddedGrid, OctaHEALPixPaddedGrid}

"""$(TYPEDSIGNATURES) The grid whose points nest in the padded grid `Grid`."""
unpadded_grid_type(::Type{<:HEALPixPaddedGrid}) = HEALPixGrid
unpadded_grid_type(::Type{<:OctaHEALPixPaddedGrid}) = OctaHEALPixGrid

## SIZE
nlat_odd(::Type{<:PaddedHEALPixGrids}) = true

# Rings with fewer than 20 longitudes are refined by the smallest integer factor k giving at least
# 20 points (k = cld(20, 4j)), these are always the 4 rings around each pole with 4, 8, 12, 16
# points, given nlat_half >= nlat_half_min, with HEALPix point i on padded point k(i-1) + first
const PADDED_NLONS = (20, 24, 24, 32)               # longitudes on rings 1-4 (and mirrored)
const PADDED_REFINEMENT = (5, 3, 2, 2)              # refinement factor k on rings 1-4
const PADDED_FIRST = (3, 2, 2, 2)                   # in-ring index of the first HEALPix point
const PADDED_EXTRA_BEFORE = (0, 16, 32, 44)         # extra points north of rings 1-4
const PADDED_EXTRA_RING = (16, 16, 12, 16)          # extra points on rings 1-4
const PADDED_EXTRA_HEMISPHERE = 60                  # extra points per hemisphere

"""$(TYPEDSIGNATURES) Minimum `nlat_half` of the padded grids such that only the 4 rings around
each pole are padded and the Equator has the most longitudes (at least 32)."""
nlat_half_min(::Type{<:HEALPixPaddedGrid}) = 16
nlat_half_min(::Type{<:OctaHEALPixPaddedGrid}) = 8

# number of longitudes on the unpadded grid, without the @assert of get_nlon_per_ring for kernels
@inline nlon_unpadded(::Type{<:HEALPixPaddedGrid}, nlat_half, j) = min(4j, 2nlat_half, 8nlat_half - 4j)
@inline nlon_unpadded(::Type{<:OctaHEALPixPaddedGrid}, nlat_half, j) = min(4j, 8nlat_half - 4j)

# ring index counted from the nearest pole, 1, 2, 3, 4 for the padded rings
@inline polar_ring(nlat_half, j) = min(j, 2nlat_half - j)

"""$(TYPEDSIGNATURES) Refinement factor `k` and in-ring index `first` of the first unpadded
point on ring `j` of the padded `Grid`, such that unpadded point `i` (in-ring index) is padded
point `k(i-1) + first`."""
@inline function healpix_nesting(::Type{<:PaddedHEALPixGrids}, nlat_half::Integer, j::Integer)
    m = polar_ring(nlat_half, j)
    return m <= 4 ? (PADDED_REFINEMENT[m], PADDED_FIRST[m]) : (1, 1)
end

"""$(TYPEDSIGNATURES) Number of points the padded `Grid` has in addition to its unpadded
grid north of ring `j`."""
@inline function extra_points_before(::Type{<:PaddedHEALPixGrids}, nlat_half::Integer, j::Integer)
    m = polar_ring(nlat_half, j)
    m > 4 && return PADDED_EXTRA_HEMISPHERE
    # northern polar rings, or southern ones counting backwards from the total extra points
    return j == m ? PADDED_EXTRA_BEFORE[m] : 2PADDED_EXTRA_HEMISPHERE - PADDED_EXTRA_BEFORE[m] - PADDED_EXTRA_RING[m]
end

function get_nlon_per_ring(Grid::Type{<:PaddedHEALPixGrids}, nlat_half::Integer, j::Integer)
    nlat = get_nlat(Grid, nlat_half)
    @assert 0 < j <= nlat "Ring $j is outside $(nonparametric_type(Grid)) with nlat_half=$nlat_half."
    m = polar_ring(nlat_half, j)
    return m <= 4 ? PADDED_NLONS[m] : nlon_unpadded(Grid, nlat_half, j)
end

get_npoints(Grid::Type{<:PaddedHEALPixGrids}, nlat_half::Integer) =
    get_npoints(unpadded_grid_type(Grid), nlat_half) + 2PADDED_EXTRA_HEMISPHERE

get_nlat_half(Grid::Type{<:PaddedHEALPixGrids}, npoints::Integer) =
    get_nlat_half(unpadded_grid_type(Grid), npoints - 2PADDED_EXTRA_HEMISPHERE)

## COORDINATES
get_latd(Grid::Type{<:PaddedHEALPixGrids}, nlat_half::Integer) = get_latd(unpadded_grid_type(Grid), nlat_half)

# offset of dlon/2 east of 0˚ if the unpadded ring has one and the refinement factor is odd
function hasoffset(Grid::Type{<:PaddedHEALPixGrids}, nlat_half::Integer, j::Integer)
    k, _ = healpix_nesting(Grid, nlat_half, j)
    return hasoffset(unpadded_grid_type(Grid), nlat_half, j) && isodd(k)
end

hasoffset(Grid::Type{<:PaddedHEALPixGrids}) = throw(
    ArgumentError(
        "$(nonparametric_type(Grid)) has rings with and without longitudinal offset, " *
            "use hasoffset(Grid, nlat_half, j) per ring instead."
    )
)

# ring-dependent offset, so kernels need it as a vector with one Bool per ring
function offset_maybe_vector(arch, grid::PaddedHEALPixGrids)
    (; nlat_half) = grid
    return on_architecture(arch, [hasoffset(nonparametric_type(grid), nlat_half, j) for j in 1:get_nlat(grid)])
end

function get_lond_per_ring(Grid::Type{<:PaddedHEALPixGrids}, nlat_half::Integer, j::Integer)
    nlon = get_nlon_per_ring(Grid, nlat_half, j)
    s = hasoffset(Grid, nlat_half, j) ? 1 : 2       # same convention as HEALPixGrid
    return [180 / (nlon ÷ 2) * (i - s / 2) for i in 1:nlon]
end

## QUADRATURE, same ring areas as HEALPix, the solid angle per point (generic
# get_solid_angles) then divides by the padded nlon
get_quadrature_weights(Grid::Type{<:PaddedHEALPixGrids}, nlat_half::Integer) =
    equal_area_weights(unpadded_grid_type(Grid), nlat_half)

## INDEXING
# closed form from the unpadded grid's ring indices shifted by the extra points north of ring j
@inline function each_index_in_ring(Grid::Type{<:PaddedHEALPixGrids}, j::Integer, nlat_half::Integer)
    index_1st = first(each_index_in_ring(unpadded_grid_type(Grid), j, nlat_half)) +
        extra_points_before(Grid, nlat_half, j)
    m = polar_ring(nlat_half, j)
    nlon = m <= 4 ? PADDED_NLONS[m] : nlon_unpadded(Grid, nlat_half, j)
    return index_1st:(index_1st + nlon - 1)
end

function each_index_in_ring!(
        rings::AbstractVector,
        Grid::Type{<:PaddedHEALPixGrids},
        nlat_half::Integer
    )
    nlat = length(rings)
    @boundscheck nlat == get_nlat(Grid, nlat_half) || throw(BoundsError)
    nlat_half >= nlat_half_min(Grid) || throw(
        DomainError(nlat_half, "$(nonparametric_type(Grid)) requires nlat_half >= $(nlat_half_min(Grid)).")
    )
    Grid <: HEALPixPaddedGrid && nside_healpix(nlat_half)     # throws for odd nlat_half

    index_end = 0
    @inbounds for j in 1:nlat
        index_1st = index_end + 1                   # 1st index is +1 from prev ring's last index
        index_end += get_nlon_per_ring(Grid, nlat_half, j)
        rings[j] = index_1st:index_end              # turn into UnitRange
    end
    return rings
end

## CONVERSION BETWEEN (Octa)HEALPixGrid AND (Octa)HEALPixPaddedGrid

function check_nesting(healpix::AbstractField, padded::AbstractField)
    @assert healpix.grid.nlat_half == padded.grid.nlat_half "Unpadded and padded grids differ in nlat_half: " *
        "$(healpix.grid.nlat_half) vs $(padded.grid.nlat_half)"
    @assert size(healpix)[2:end] == size(padded)[2:end] "Non-horizontal dimensions differ: " *
        "$(size(healpix)[2:end]) vs $(size(padded)[2:end])"
    return nothing
end

"""$(TYPEDSIGNATURES) Copy `padded` onto `healpix` by subsampling: every `(Octa)HEALPixGrid`
point is also a `(Octa)HEALPixPaddedGrid` point, so this is exact."""
Base.copyto!(healpix::HEALPixField, padded::HEALPixPaddedField) = subsample_padded!(healpix, padded)
Base.copyto!(healpix::OctaHEALPixField, padded::OctaHEALPixPaddedField) = subsample_padded!(healpix, padded)

function subsample_padded!(healpix::AbstractField, padded::AbstractField)
    check_nesting(healpix, padded)
    Grid = nonparametric_type(padded.grid)        # as Val for the kernel to specialize on it
    launch!(
        architecture(healpix), RingGridWorkOrder, size(healpix), subsample_padded_kernel!,
        healpix, padded, healpix.grid.whichring, Val(Grid), healpix.grid.nlat_half
    )
    return healpix
end

@kernel inbounds = true function subsample_padded_kernel!(healpix, padded, whichring, ::Val{Grid}, nlat_half) where {Grid}
    I = @index(Global, Cartesian)
    ij = I[1]                                   # grid point index on the unpadded grid
    k = CartesianIndex(Base.tail(Tuple(I)))     # all non-horizontal dimensions
    j = whichring[ij]
    ring_first = first(each_index_in_ring(unpadded_grid_type(Grid), j, nlat_half))
    refinement, first_nested = healpix_nesting(Grid, nlat_half, j)

    # same ring on the padded grid starts later by the extra points north of it
    ij_padded = ring_first + extra_points_before(Grid, nlat_half, j) + refinement * (ij - ring_first) + first_nested - 1
    healpix[ij, k] = padded[ij_padded, k]
end

"""$(TYPEDSIGNATURES) Copy `healpix` onto `padded`, nested points are copied, the points in
between are linearly interpolated in longitude along each ring (the latitudes are identical)."""
Base.copyto!(padded::HEALPixPaddedField, healpix::HEALPixField) = interpolate_padded!(padded, healpix)
Base.copyto!(padded::OctaHEALPixPaddedField, healpix::OctaHEALPixField) = interpolate_padded!(padded, healpix)

function interpolate_padded!(padded::AbstractField, healpix::AbstractField)
    check_nesting(healpix, padded)
    Grid = nonparametric_type(padded.grid)        # as Val for the kernel to specialize on it
    launch!(
        architecture(padded), RingGridWorkOrder, size(padded), interpolate_padded_kernel!,
        padded, healpix, padded.grid.whichring, Val(Grid), padded.grid.nlat_half
    )
    return padded
end

@kernel inbounds = true function interpolate_padded_kernel!(padded, healpix, whichring, ::Val{Grid}, nlat_half) where {Grid}
    I = @index(Global, Cartesian)
    ij = I[1]                                   # grid point index on the padded grid
    k = CartesianIndex(Base.tail(Tuple(I)))     # all non-horizontal dimensions
    j = whichring[ij]
    ring_first = first(each_index_in_ring(unpadded_grid_type(Grid), j, nlat_half))
    nlon = nlon_unpadded(Grid, nlat_half, j)
    refinement, first_nested = healpix_nesting(Grid, nlat_half, j)

    # padded in-ring index relative to the first nested point, in units of padded points
    d = ij - (ring_first + extra_points_before(Grid, nlat_half, j)) + 1 - first_nested
    q = fld(d, refinement)                      # 0-based unpadded point west of (or on) ij
    r = d - q * refinement                      # padded points east of that unpadded point
    i_west = ring_first + mod(q, nlon)
    i_east = ring_first + mod(q + 1, nlon)
    w = convert(eltype(padded), r) / refinement # linear interpolation weight
    padded[ij, k] = (1 - w) * healpix[i_west, k] + w * healpix[i_east, k]
end
