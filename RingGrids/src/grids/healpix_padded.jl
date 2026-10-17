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
in longitude. The grid uses the equal-area ring weights of HEALPix, but is not equal-area
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
subsamples or linearly interpolates in longitude. Fields are
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

"""$(TYPEDSIGNATURES) Minimum number of longitude points per ring of a `HEALPixPaddedGrid`."""
npoints_min(::Type{<:PaddedHEALPixGrids}) = 20

"""$(TYPEDSIGNATURES) Refinement factor `k` of ring `j` relative to the `HEALPixGrid`, i.e.
the padded ring has `k` times the HEALPix longitudes."""
function refinement_factor(Grid::Type{<:PaddedHEALPixGrids}, nlat_half::Integer, j::Integer)
    nlon_healpix = get_nlon_per_ring(unpadded_grid_type(Grid), nlat_half, j)
    return cld(npoints_min(Grid), nlon_healpix)
end

function get_nlon_per_ring(Grid::Type{<:PaddedHEALPixGrids}, nlat_half::Integer, j::Integer)
    return refinement_factor(Grid, nlat_half, j) * get_nlon_per_ring(unpadded_grid_type(Grid), nlat_half, j)
end

get_npoints(Grid::Type{<:PaddedHEALPixGrids}, nlat_half::Integer) =
    nlat_half == 0 ? 0 : sum(get_nlon_per_ring(Grid, nlat_half, j) for j in 1:get_nlat(Grid, nlat_half))

function get_nlat_half(Grid::Type{<:PaddedHEALPixGrids}, npoints::Integer)
    # no closed form, so search, npoints grows monotonically with nlat_half
    step = Grid <: HEALPixPaddedGrid ? 2 : 1          # HEALPix requires even nlat_half
    nlat_half = step
    while get_npoints(Grid, nlat_half) < npoints
        nlat_half += step
    end
    return nlat_half
end

# for small grids the refined polar rings can have more points than the Equator
get_nlon_max(Grid::Type{<:PaddedHEALPixGrids}, nlat_half::Integer) =
    maximum(get_nlon_per_ring(Grid, nlat_half, j) for j in 1:nlat_half)

## COORDINATES
get_latd(Grid::Type{<:PaddedHEALPixGrids}, nlat_half::Integer) = get_latd(unpadded_grid_type(Grid), nlat_half)

# offset of dlon/2 east of 0˚ if HEALPix has one on that ring and the refinement factor is odd
function hasoffset(Grid::Type{<:PaddedHEALPixGrids}, nlat_half::Integer, j::Integer)
    return hasoffset(unpadded_grid_type(Grid), nlat_half, j) && isodd(refinement_factor(Grid, nlat_half, j))
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

"""$(TYPEDSIGNATURES) Refinement factor `k` and in-ring index `first` of the first HEALPix
point on ring `j` of the `HEALPixPaddedGrid`, such that HEALPix point `i` (in-ring index)
is padded point `k(i-1) + first`."""
function healpix_nesting(Grid::Type{<:PaddedHEALPixGrids}, nlat_half::Integer, j::Integer)
    k = refinement_factor(Grid, nlat_half, j)
    # HEALPix point 1 is at dlon_healpix/2 = k*dlon/2 (offset) or 0˚ (no offset)
    first = hasoffset(unpadded_grid_type(Grid), nlat_half, j) ? (isodd(k) ? (k + 1) ÷ 2 : k ÷ 2 + 1) : 1
    return k, first
end

## QUADRATURE, same ring areas as HEALPix, the solid angle per point (generic
# get_solid_angles) then divides by the padded nlon
get_quadrature_weights(Grid::Type{<:PaddedHEALPixGrids}, nlat_half::Integer) =
    equal_area_weights(unpadded_grid_type(Grid), nlat_half)

## INDEXING
function each_index_in_ring!(
        rings::AbstractVector,
        Grid::Type{<:PaddedHEALPixGrids},
        nlat_half::Integer
    )
    nlat = length(rings)
    @boundscheck nlat == get_nlat(Grid, nlat_half) || throw(BoundsError)
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

# refinement factor k and first nested index per ring on architecture arch, for kernels
function healpix_nesting_on_architecture(arch, grid::PaddedHEALPixGrids)
    nesting = [healpix_nesting(nonparametric_type(grid), grid.nlat_half, j) for j in 1:get_nlat(grid)]
    return on_architecture(arch, first.(nesting)), on_architecture(arch, last.(nesting))
end

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
    arch = architecture(healpix)
    ring_first_healpix, _ = eachring_on_architecture(arch, healpix.grid)
    ring_first_padded, _ = eachring_on_architecture(arch, padded.grid)
    ks, firsts = healpix_nesting_on_architecture(arch, padded.grid)
    launch!(
        arch, RingGridWorkOrder, size(healpix), subsample_healpix_kernel!, healpix, padded,
        healpix.grid.whichring, ring_first_healpix, ring_first_padded, ks, firsts
    )
    return healpix
end

@kernel inbounds = true function subsample_healpix_kernel!(
        healpix, padded, whichring, ring_first_healpix, ring_first_padded, ks, firsts
    )
    I = @index(Global, Cartesian)
    ij = I[1]                                   # grid point index on HEALPix
    k = CartesianIndex(Base.tail(Tuple(I)))     # all non-horizontal dimensions
    j = whichring[ij]
    i = ij - ring_first_healpix[j]              # 0-based in-ring index on HEALPix
    healpix[ij, k] = padded[ring_first_padded[j] + ks[j] * i + firsts[j] - 1, k]
end

"""$(TYPEDSIGNATURES) Copy `healpix` onto `padded`, nested points are copied, the points in
between are linearly interpolated in longitude along each ring (the latitudes are identical)."""
Base.copyto!(padded::HEALPixPaddedField, healpix::HEALPixField) = interpolate_padded!(padded, healpix)
Base.copyto!(padded::OctaHEALPixPaddedField, healpix::OctaHEALPixField) = interpolate_padded!(padded, healpix)

function interpolate_padded!(padded::AbstractField, healpix::AbstractField)
    check_nesting(healpix, padded)
    arch = architecture(padded)
    ring_first_healpix, ring_length_healpix = eachring_on_architecture(arch, healpix.grid)
    ring_first_padded, _ = eachring_on_architecture(arch, padded.grid)
    ks, firsts = healpix_nesting_on_architecture(arch, padded.grid)
    launch!(
        arch, RingGridWorkOrder, size(padded), interpolate_healpix_padded_kernel!, padded, healpix,
        padded.grid.whichring, ring_first_healpix, ring_length_healpix, ring_first_padded, ks, firsts
    )
    return padded
end

@kernel inbounds = true function interpolate_healpix_padded_kernel!(
        padded, healpix, whichring, ring_first_healpix, ring_length_healpix, ring_first_padded, ks, firsts
    )
    I = @index(Global, Cartesian)
    ij = I[1]                                   # grid point index on HEALPixPadded
    k = CartesianIndex(Base.tail(Tuple(I)))     # all non-horizontal dimensions
    j = whichring[ij]
    nlon = ring_length_healpix[j]

    # padded in-ring index relative to the first nested point, in units of padded points
    d = ij - ring_first_padded[j] + 1 - firsts[j]
    q = fld(d, ks[j])                           # 0-based HEALPix point west of (or on) ij
    r = d - q * ks[j]                           # padded points east of that HEALPix point
    i_west = ring_first_healpix[j] + mod(q, nlon)
    i_east = ring_first_healpix[j] + mod(q + 1, nlon)
    w = convert(eltype(padded), r) / ks[j]      # linear interpolation weight
    padded[ij, k] = (1 - w) * healpix[i_west, k] + w * healpix[i_east, k]
end
