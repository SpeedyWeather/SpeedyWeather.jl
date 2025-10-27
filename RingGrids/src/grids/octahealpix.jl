"""An `OctaHEALPixGrid` is an equal-area discretization of the sphere.
It has 4 faces, each `nlat_half x nlat_half` in size,
covering 90˚ in longitude, pole to pole. As part of the HEALPix family of grids,
the grid points are equal area. They start with 4 longitude points on the northern-most ring,
increase by 4 points per ring  towards the Equator with one ring on the Equator before reducing
the number of points again towards the south pole by 4 per ring. There is no equatorial belt for
OctaHEALPix grids. The southern hemisphere is symmetric to the northern, mirrored around the Equator.
OctaHEALPix grids have a ring on the Equator. For more details see
Górski et al. 2005, DOI:10.1086/427976, the OctaHEALPix grid belongs to the family of
HEALPix grids with Nθ = 1, Nφ = 4 but is not explicitly mentioned therein.

`rings` are the precomputed ring indices, for nlat_half = 3 (in contrast to HEALPix this can be odd)
it is `rings = [1:4, 5:12, 13:24, 25:32, 33:36]`. For efficient looping see `eachring` and `eachgrid`.
`whichring` is a precomputed vector of ring indices for each grid point ij, i.e. `whichring[ij]` gives
the ring index j of grid point ij. Fields are
$(TYPEDFIELDS)"""
struct OctaHEALPixGrid{A, V, W} <: AbstractReducedGrid{A}
    nlat_half::Int      # number of latitudes on one hemisphere
    architecture::A     # information about device, CPU/GPU
    rings::V            # precomputed ring indices
    whichring::W        # precomputed ring index for each grid point ij
end

nonparametric_type(::Type{<:OctaHEALPixGrid}) = OctaHEALPixGrid
full_grid_type(::Type{<:OctaHEALPixGrid}) = FullOctaHEALPixGrid

# FIELD
const OctaHEALPixField{T, N} = Field{T, N, ArrayType, Grid} where {ArrayType, Grid<:OctaHEALPixGrid}

# define grid_type (i) without T, N, (ii) with T, (iii) with T, N but not with <:?Field
# to not have precendence over grid_type(::Type{Field{...})
grid_type(::Type{OctaHEALPixField}) = OctaHEALPixGrid
grid_type(::Type{OctaHEALPixField{T}}) where T = OctaHEALPixGrid
grid_type(::Type{OctaHEALPixField{T, N}}) where {T, N} = OctaHEALPixGrid

function Base.showarg(io::IO, F::Field{T, N, ArrayType, Grid}, toplevel) where {T, N, ArrayType, Grid<:OctaHEALPixGrid{A}} where A <: AbstractArchitecture
    print(io, "OctaHEALPixField{$T, $N}")
    toplevel && print(io, " as ", nonparametric_type(ArrayType))
    toplevel && print(io, " on ", F.grid.architecture)
end

## SIZE
nlat_odd(::Type{<:OctaHEALPixGrid}) = true
get_npoints(::Type{<:OctaHEALPixGrid}, nlat_half::Integer) = 4*nlat_half^2
get_nlat_half(::Type{<:OctaHEALPixGrid}, npoints::Integer) = round(Int, sqrt(npoints/4))

# number of longitude 
function get_nlon_per_ring(Grid::Type{<:OctaHEALPixGrid}, nlat_half::Integer, j::Integer)
    nlat = get_nlat(Grid, nlat_half)
    @assert 0 < j <= nlat "Ring $j is outside P$nlat_half grid."
    # j = j > nlat_half ? nlat - j + 1 : j      # flip north south due to symmetry
    return min(4j, 8nlat_half-4j)
end

matrix_size(::Type{<:OctaHEALPixGrid}, nlat_half::Integer) = (2nlat_half, 2nlat_half)

## COORDINATES
function get_latd(::Type{<:OctaHEALPixGrid}, nlat_half::Integer)
    nlat = get_nlat(OctaHEALPixGrid, nlat_half)
    latd = zeros(nlat)

    # Górski et al. 2005 eq 4 but without the 1/3 and Nside=nlat_half
    for j in 1:nlat_half        latd[j] = 90 - acosd(1-(j/nlat_half)^2) end # north + Equator
    for j in nlat_half+1:nlat   latd[j] = -latd[nlat-j+1]               end # southern hemisphere

    return latd
end

function get_lond_per_ring(Grid::Type{<:OctaHEALPixGrid}, nlat_half::Integer, j::Integer)
    nlon = get_nlon_per_ring(Grid, nlat_half, j)
    # equidistant longitudes with equal offsets from 0˚ and 360˚,
    # e.g. 45, 135, 225, 315 for nlon=4
    return collect(180/nlon:360/nlon:360)
end

## INDEXING
function each_index_in_ring(::Type{<:OctaHEALPixGrid},     # function for OctaHEALPix grids
                            j::Integer,                     # ring index north to south
                            nlat_half::Integer)             # resolution param

    @boundscheck 0 < j < 2nlat_half || throw(BoundsError)   # ring index valid?
    if j <= nlat_half                                       # northern hemisphere incl Equator
        index_1st = 2j*(j-1) + 1                            # first in-ring index i
        index_end = 2j*(j+1)                                # last in-ring index i
    else                                                    # southern hemisphere 
        n = 4nlat_half^2                                    # total number of points
        j = 2nlat_half - j                                  # count ring index from south pole
        index_1st = n - 2j*(j+1) + 1                        # count backwards
        index_end = n - 2j*(j-1)
    end
    return index_1st:index_end                              # range of i's in ring
end

function each_index_in_ring!(   rings,
                                Grid::Type{<:OctaHEALPixGrid},
                                nlat_half::Integer) # resolution param

    nlat = length(rings)
    @boundscheck nlat == get_nlat(Grid, nlat_half) || throw(BoundsError)

    index_end = 0
    @inbounds for j in 1:nlat_half                  # North incl Eq only
        index_1st = index_end + 1                   # 1st index is +1 from prev ring's last index
        index_end += 4j                             # add number of grid points per ring
        rings[j] = index_1st:index_end              # turn into UnitRange
    end
    @inbounds for (j, j_rev) in zip(nlat_half+1:nlat,       # South only
                                    nlat-nlat_half:-1:1)    # reverse index

        index_1st = index_end + 1                   # 1st index is +1 from prev ring's last index
        index_end += 4j_rev                         # add number of grid points per ring
        rings[j] = index_1st:index_end              # turn into UnitRange
    end
end

Adapt.@adapt_structure OctaHEALPixGrid

# REORDERING: Ring to Nested or Matrix and vice versa

abstract type AbstractHEALPixOrder end
struct RingOrder <: AbstractHEALPixOrder end
struct NestedOrder <: AbstractHEALPixOrder end
struct MatrixOrder <: AbstractHEALPixOrder end

reorder(::RingOrder,   ij, grid::OctaHEALPixGrid) = nest2ring(ij, grid)
reorder(::NestedOrder, ij, grid::OctaHEALPixGrid) = ring2nest(ij, grid)
reorder(::MatrixOrder, ij, grid::OctaHEALPixGrid) = ring2xy(ij, grid)

reorder(order, field::OctaHEALPixField) = reorder!(similar(field), order, field)

function reorder!(
    out::OctaHEALPixField,
    order,
    field::OctaHEALPixField,
)
    @boundscheck out.grid == field.grid || throw(BoundsError("Reordering requires identical grids, got $(out.grid) and $(field.grid)."))
    @assert ispow2(get_nlat_half(field.grid)) "Reordering only supported for nlat_half power of 2, got $(get_nlat_half(field.grid))."

    @inbounds for k in eachlayer(field)
        for ij in eachgridpoint(field)
            out_indices = reorder(order, ij, field.grid)
            out[out_indices, k] = field[ij, k]
        end    
    end
    return out
end

ring_order(  field::OctaHEALPixField) = reorder(RingOrder(),   field)
nested_order(field::OctaHEALPixField) = reorder(NestedOrder(), field)
matrix_order(field::OctaHEALPixField) = reorder(MatrixOrder(), field)

ring_order!(  out, field::OctaHEALPixField) = reorder!(out, RingOrder(),   field)
nested_order!(out, field::OctaHEALPixField) = reorder!(out, NestedOrder(), field)
matrix_order!(out, field::OctaHEALPixField) = reorder!(out, MatrixOrder(), field)

Matrix(field::OctaHEALPixField) = reshape(field.data, matrix_size(field)...)

# quadrant of ij in ring order, TODO needed?
function quadrant_ring(ij::Integer, grid::OctaHEALPixGrid)
    j = RingGrids.whichring(grid)[ij]  # ring index j of ij
    ring = eachring(grid)[j]           # ij indices of ring j
    nlon = length(ring)                # number of grid points in ring
    i = ij - ring[1]                   # 0-based index in ring
    q = mod(4i ÷ nlon, 4)              # quadrant q, either 0, 1, 2, 3
    iq = i - q*(nlon÷4)                # 0-based index i relative to quadrant
    return q+1, iq+1                   # convert to 1-based
end

"""$TYPEDSIGNATURES
Convert ring index ij to matrix indices row, column (r, c) and quadrant q. All 1-based.
r=1, c=1, is at the north pole, r increases south-eastwards, c increases south-westwards.
Quadrants are numbered 1 to 4 starting at the prime meridian and increasing eastwards."""
function ring2rcq(ij::Integer, grid::OctaHEALPixGrid)
    nside = get_nlat_half(grid)     # resolution param, nside=nlat_half for OctaHEALPix
    j = whichring(grid)[ij]         # ring index j of ij
    ring = eachring(grid)[j]        # ij indices of ring j
    nlon = length(ring)             # number of grid points in ring
    i = ij - ring[1]                # 0-based index in ring
    q = mod(4i ÷ nlon, 4)           # quadrant q, either 0, 1, 2, 3
    iq = i - q*(nlon÷4)             # 0-based index i but relative to quadrant
    q += 1; iq += 1                 # convert to 1-based
    r = min(j, nside) - iq + 1      # row in matrix m (1-based)
    c = iq + max(0, j-nside)        # column in matrix m (1-based)
    return r, c, q
end 

"""$TYPEDSIGNATURES
Convert matrix indices row, column (r, c) and quadrant q to ring index ij. All 1-based.
r=1, c=1, is at the north pole, r increases south-eastwards, c increases south-westwards.
Quadrants are numbered 1 to 4 starting at the prime meridian and increasing eastwards."""
function rcq2ring(r, c, q, grid::OctaHEALPixGrid)
    j = r + c - 1                   # 1-based ring index
    nside = get_nlat_half(grid)     # resolution param, nside=nlat_half for OctaHEALPix
    ring = eachring(grid)[j]        # ij indices of latitude ring j
    iq = min(nside - r, c - 1) + 1  # 1-based in-ring index i, relative to quadrant
    nlon = length(ring)             # number of longitude points in ring
    i = iq  + (q - 1)*(nlon ÷ 4)    # in-ring index i (1-based)
    ij = ring[i]                    # convert to running index ij
    return ij
end 

# unpack nlat_half from grid
rcq2nest(r, c, q, grid::OctaHEALPixGrid) = rcq2nest(r, c, q, get_nlat_half(grid))


"""$TYPEDSIGNATURES
Convert matrix indices row, column (r, c) and quadrant q to nested index ij. All 1-based.
r=1, c=1, is at the north pole, r increases south-eastwards, c increases south-westwards.
Quadrants are numbered 1 to 4 starting at the prime meridian and increasing eastwards."""
function rcq2nest(r, c, q, nside::Integer)
    # @assert ispow2(nside)                   # nside must be power of 2
    # TODO UInt32 restricts nside to 2^16, ~150m resolution, remove for higher resolution but slower
    br = interleave_with_zeros(UInt32(r-1))         # bit representation of row (0-based) occupies odd bits
    bc = interleave_with_zeros(UInt32(c-1)) << 1    # bit representation of column (0-based) occupies even bits
    bq = (q-1) << 2trailing_zeros(nside)    # the 2 quadrant bits occupy highest bits
    b = bq | bc | br                        # combine bits
    return b + 1                            # convert to 1-based nested index
end

# unpack nlat_half from grid
nest2rcq(ij, grid::OctaHEALPixGrid) = nest2rcq(ij, get_nlat_half(grid))

"""$TYPEDSIGNATURES
Convert nested index ij to matrix indices row, column (r, c) and quadrant q. All 1-based.
r=1, c=1, is at the north pole, r increases south-eastwards, c increases south-westwards.
Quadrants are numbered 1 to 4 starting at the prime meridian and increasing eastwards."""
function nest2rcq(ij, nside)
    # @assert ispow2(nside)           # nside must be power of 2
    ij -= 1                         # convert to 0-based ij
    shift = 2trailing_zeros(nside)  # number of bits occupied for row and column
    q = ij >> shift                 # 0-based quadrant is encoded in first 2 bits
    bcr = ij - (q << shift)         # bits for column, row, quadrant bits removed
    bcr %= UInt32                   # TODO this restricts to nside = 2^16, ~150m
                                    # remove for higher resolution but also much slower
    r = deinterleave(bcr) + 1       # 1-based row index by deinterleaving the odd bits
    c = deinterleave(bcr >> 1) + 1  # 1-based column index by deinterleaving the even bits
    q += 1                          # 1-based quadrant
    return r, c, q
end

"""$TYPEDSIGNATURES
Interleave an integer's bits with zeros, e.g.

    00000000 00000000 00000000 00000111

becomes
    
    00000000 00000000 00000000 00010101

The highest half of the bits are zeros will be discarded."""
function interleave_with_zeros(ui::Integer)
    r = zero(ui)
    nbits = 8*sizeof(ui)
    for s in 0:nbits-1  # TODO is there a more efficient way?
        r |= (ui & (one(ui) << s)) << s
    end
    return r
end

"""$TYPEDSIGNATURES
Deinterleave an integer's bits, e.g.    

    00000000 00000000 00000000 00010101

becomes

    00000000 00000000 00000000 00000111."""
function deinterleave(ui::Integer)
    r = zero(ui)
    nbits = 8*sizeof(ui)
    for s in 0:nbits-1  # TODO is there a more efficient way?
        r |= (ui & (one(ui) << 2s)) >> s
    end
    return r
end

"""$TYPEDSIGNATURES
Convert ring index ij to nested index ij of grid. All 1-based."""
ring2nest(ij::Integer, grid::OctaHEALPixGrid) = rcq2nest(ring2rcq(ij, grid)..., grid)

"""$TYPEDSIGNATURES
Convert nested index ij to ring index ij of grid. All 1-based."""
nest2ring(ij::Integer, grid::OctaHEALPixGrid) = rcq2ring(nest2rcq(ij, grid)..., grid)

rcq2xy(r, c, q, grid::OctaHEALPixGrid; kwargs...) = rcq2xy(r, c, q, get_nlat_half(grid); kwargs...)

function rcq2xy(r, c, q, nside;
    quadrant_rotation = (0, 1, 2, 3),                    # = 0˚, 90˚, 180˚, 270˚ anti-clockwise
    matrix_quadrant = ((2, 2), (1, 2), (1, 1), (2, 1)),  # north polar-centric view
)
    # rotate indices in quadrant
    x, y = rotate_matrix_indices(r, c, nside, quadrant_rotation[q])

    # shift grid quadrant to matrix quadrant
    sr, sc = matrix_quadrant[q]
    x += (sr-1)*nside           # shift row into matrix quadrant
    y += (sc-1)*nside           # shift column into matrix quadrant
    xy = (y-1)*2nside + x       # to single running index xy that can be reshaped into matrix
    return xy
end

"""$TYPEDSIGNATURES
Convert ring index ij to matrix index xy of grid. All 1-based.
xy is a running index in a 2D matrix of size (2*nlat_half, 2*nlat_half),
with a polar-centric view on the north pole in the middle of that matrix,
the South Pole divided into 4 in the corners. Like a stereographic projection."""
ring2xy(grid::OctaHEALPixGrid, ij::Integer; kwargs...) = rcq2xy(ring2rcq(ij, grid)..., grid; kwargs...)