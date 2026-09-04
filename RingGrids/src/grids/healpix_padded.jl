"""A `HEALPixPaddedGrid` shares every latitude with the `HEALPixGrid` but carries more longitude
points on the polar-cap rings. HEALPix's cap ring `j` has `4j` longitudes, so the innermost rings
hold 4, 8, 12, ... points while the fields on them still carry spherical-harmonic orders well above
their Nyquist limit `(4j-1)÷2`. The content above that limit is simply lost when a product is
formed on the grid, and because the loss is one-signed for a squared quantity it acts as a
systematic polar sink that no quadrature weight can repair — see
`healpix_quadrature/polar_truncation.jl`, where the entire zonal-mean error of an analysed product
sits on the innermost 8 rings of each hemisphere.

This grid removes that error at its source by giving every cap ring `min(4j + npoints_padding, nlon_belt)`
longitudes, in the spirit of the `OctahedralGaussianGrid`'s 16 extra points per ring. The
equatorial belt, all latitudes, and the ring longitude offsets are HEALPix's. It is therefore
**not** an equal-area grid: cap pixels are smaller than belt pixels. What it preserves is the
HEALPix *ring* structure, so a field on it maps onto true HEALPix pixels by a per-ring Fourier
resampling, which is exact for band-limited fields.

`rings` are the precomputed ring indices and `whichring[ij]` the ring index of grid point `ij`,
as for every reduced grid. Fields are
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

grid_type(::Type{HEALPixPaddedField}) = HEALPixPaddedGrid
grid_type(::Type{HEALPixPaddedField{T}}) where {T} = HEALPixPaddedGrid
grid_type(::Type{HEALPixPaddedField{T, N}}) where {T, N} = HEALPixPaddedGrid

function Base.showarg(io::IO, F::Field{T, N, ArrayType, Grid}, toplevel) where {T, N, ArrayType, Grid <: HEALPixPaddedGrid{A}} where {A <: AbstractArchitecture}
    print(io, "HEALPixPaddedField{$T, $N}")
    toplevel && print(io, " as ", nonparametric_type(ArrayType))
    return toplevel && print(io, " on ", F.grid.architecture)
end

## SIZE

"""$(TYPEDSIGNATURES) Longitude points added to every polar-cap ring of a `HEALPixPaddedGrid`,
capped at the belt's `nlon`. 16 matches the `OctahedralGaussianGrid`, whose rings start at 20
rather than 4 points for exactly this reason. Must be a multiple of 4 so that a cap ring keeps an
even `nlon` (the real FFT needs it) and the HEALPix half-pixel offset stays meaningful."""
npoints_padding(::Type{<:HEALPixPaddedGrid}) = 16

nlat_odd(::Type{<:HEALPixPaddedGrid}) = true

function get_nlon_per_ring(Grid::Type{<:HEALPixPaddedGrid}, nlat_half::Integer, j::Integer)
    nlat = get_nlat(Grid, nlat_half)
    @assert 0 < j <= nlat "Ring $j is outside H$nlat_half grid."
    nlon_belt = 2nlat_half                                  # HEALPix belt, unchanged
    nlon_healpix = min(4j, nlon_belt, 8nlat_half - 4j)      # what HEALPixGrid would use
    # pad the cap rings only, never past the belt's nlon
    return min(nlon_healpix + npoints_padding(Grid), nlon_belt)
end

function get_npoints(Grid::Type{<:HEALPixPaddedGrid}, nlat_half::Integer)
    nlat_half == 0 && return 0
    return sum(get_nlon_per_ring(Grid, nlat_half, j) for j in 1:get_nlat(Grid, nlat_half))
end

function get_nlat_half(Grid::Type{<:HEALPixPaddedGrid}, npoints::Integer)
    # no closed form with the min(), so invert by search; npoints grows monotonically in nlat_half
    nlat_half = round(Int, sqrt(npoints / 3))
    nlat_half += mod(nlat_half, 2)                          # even only, as for HEALPixGrid
    while get_npoints(Grid, nlat_half) < npoints
        nlat_half += 2
    end
    while nlat_half > 2 && get_npoints(Grid, nlat_half - 2) >= npoints
        nlat_half -= 2
    end
    return nlat_half
end

## COORDINATES — HEALPix's, unchanged

get_latd(::Type{<:HEALPixPaddedGrid}, nlat_half::Integer) = get_latd(HEALPixGrid, nlat_half)

function get_lond_per_ring(Grid::Type{<:HEALPixPaddedGrid}, nlat_half::Integer, j::Integer)
    nside = nside_healpix(nlat_half)
    nlon = get_nlon_per_ring(Grid, nlat_half, j)
    # same offset convention as HEALPixGrid: s = 1 for the polar caps, alternating in the belt
    s = (j < nside) || (j >= 3nside) ? 1 : ((j - nside) % 2 + 1)
    return [180 / (nlon ÷ 2) * (i - s / 2) for i in 1:nlon]
end

function hasoffset(::Type{<:HEALPixPaddedGrid}, nlat_half::Integer, j::Integer)
    nside = nside_healpix(nlat_half)
    s = (j < nside) || (j >= 3nside) ? 1 : ((j - nside) % 2 + 1)
    return s == 1
end

hasoffset(::Type{<:HEALPixPaddedGrid}) = throw(
    ArgumentError(
        "HEALPixPaddedGrid has rings with and without longitudinal offset, " *
            "use hasoffset(Grid, nlat_half, j) per ring instead."
    )
)

function offset_maybe_vector(arch, grid::HEALPixPaddedGrid)
    (; nlat_half) = grid
    return on_architecture(arch, [hasoffset(HEALPixPaddedGrid, nlat_half, j) for j in 1:get_nlat(grid)])
end

## QUADRATURE WEIGHTS
#
# The ring *areas* are HEALPix's — the latitudes are unchanged, so each ring still subtends the
# same band of colatitude — but they are now shared among more points on the cap rings. So the
# ring weight is HEALPix's and only the per-point solid angle changes. Writing `w[j]` for the
# fraction of the sphere in ring j (summing to 2 over colatitude, as `equal_area_weights` does):
#
#     w[j] = 2 · nlon_healpix[j] / npoints_healpix
#
# which is what `equal_area_weights(HEALPixGrid, nlat_half)` returns.
get_quadrature_weights(Grid::Type{<:HEALPixPaddedGrid}, nlat_half::Integer) =
    equal_area_weights(HEALPixGrid, nlat_half)

# ...and the solid angle of one grid point is that ring weight divided by its own point count,
# not by HEALPix's. The generic `get_solid_angles` does exactly this, so it needs no method here;
# what it must NOT inherit is HEALPixGrid's constant-4π/npoints shortcut.

## INDEXING

function each_index_in_ring(
        Grid::Type{<:HEALPixPaddedGrid},
        j::Integer,                     # ring index north to south
        nlat_half::Integer
    )
    nlat = get_nlat(Grid, nlat_half)
    @boundscheck 0 < j <= nlat || throw(BoundsError)
    index_end = 0
    @inbounds for jj in 1:j
        index_end += get_nlon_per_ring(Grid, nlat_half, jj)
    end
    index_1st = index_end - get_nlon_per_ring(Grid, nlat_half, j) + 1
    return index_1st:index_end
end

function each_index_in_ring!(
        rings::AbstractVector,
        Grid::Type{<:HEALPixPaddedGrid},
        nlat_half::Integer
    )
    nlat = length(rings)
    @boundscheck nlat == get_nlat(Grid, nlat_half) || throw(BoundsError)
    @boundscheck iseven(nlat_half) ||
        throw(AssertionError("HEALPixPaddedGrid only defined for even nlat_half, nlat_half=$nlat_half provided."))

    index_end = 0
    @inbounds for j in 1:nlat
        index_1st = index_end + 1
        index_end += get_nlon_per_ring(Grid, nlat_half, j)
        rings[j] = index_1st:index_end
    end
    return rings
end
