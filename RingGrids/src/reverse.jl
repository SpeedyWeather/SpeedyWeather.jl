"""$(TYPEDSIGNATURES)
Reverse `field` in-place along dimension `sym`, called via `reverse!(field, dims=:lat)`.
`:lat` (or `:latitude`) mirrors the field about the equator,
`:lon` (or `:longitude`) reverses the direction of every ring."""
Base._reverse!(field::AbstractField, sym::Symbol) = Base._reverse!(field, Val(sym))

"""$(TYPEDSIGNATURES)
Mirror `field` about the equator in-place by swapping every northern ring
with its southern counterpart, as called via `reverse!(field, dims=:lat)`."""
function Base._reverse!(field::AbstractField, ::Val{:lat})
    arch = architecture(field)
    ring_first, ring_length = eachring_on_architecture(arch, field.grid)
    nrings = get_nlat(field)

    # one thread per northern ring j (and layer k), swapping with its mirrored southern ring;
    # for odd nrings the equator ring maps onto itself and is skipped
    worksize = (nrings ÷ 2, size(field)[2:end]...)
    launch!(arch, ArrayWorkOrder, worksize, reverse_lat_kernel!, field, nrings, ring_first, ring_length)
    return field
end

@kernel inbounds = true function reverse_lat_kernel!(field, nrings, ring_first, ring_length)
    I = @index(Global, Cartesian)
    j = I[1]                                    # northern ring index
    k = CartesianIndex(Base.tail(Tuple(I)))     # all non-horizontal dimensions
    j_south = nrings - j + 1                    # mirrored southern ring, same length by symmetry
    ij_north = ring_first[j]
    ij_south = ring_first[j_south]

    for i in 1:ring_length[j]
        field[ij_north, k], field[ij_south, k] = field[ij_south, k], field[ij_north, k]
        ij_north += 1
        ij_south += 1
    end
end

"""$(TYPEDSIGNATURES)
Reverse `field` in longitude in-place (mirror at the 0˚ meridian), as called via
`reverse!(field, dims=:lon)`. On rings with a longitudinal offset (first point dlon/2
east of 0˚, see `hasoffset`) all points within the ring are reversed; on rings whose
first point lies on 0˚ that point stays and only the remaining points are reversed.
Either way the ring is mirrored exactly at 0˚, consistent with the spectral
`reverse!(::LowerTriangularArray, dims=:lon)`."""
function Base._reverse!(field::AbstractField, ::Val{:lon})
    arch = architecture(field)
    ring_first, ring_length = eachring_on_architecture(arch, field.grid)

    # true/false if the dlon/2 offset from 0˚ is identical on every ring,
    # otherwise a Bool vector with one offset per ring (HEALPixGrid)
    offset = offset_maybe_vector(arch, field.grid)

    # one thread per ring j (and layer k) as the in-place reversal is sequential within a ring
    worksize = (get_nlat(field), size(field)[2:end]...)
    launch!(arch, ArrayWorkOrder, worksize, reverse_lon_kernel!, field, ring_first, ring_length, offset)
    return field
end

# offset rings mirror at 0˚ by reversing all points; rings starting at 0˚ keep their
# first point (it lies on the mirror axis) and reverse the others, hence + (1 - offset)
@kernel inbounds = true function reverse_lon_kernel!(field, ring_first, ring_length, offset::Bool)
    I = @index(Global, Cartesian)
    j = I[1]                                    # ring index
    k = CartesianIndex(Base.tail(Tuple(I)))     # all non-horizontal dimensions
    i0 = ring_first[j] + (1 - offset)
    reverse_ring!(field, k, i0, ring_first[j] + ring_length[j] - 1)
end

# as above but for grids where the offset is ring-dependent (HEALPixGrid)
@kernel inbounds = true function reverse_lon_kernel!(field, ring_first, ring_length, offset::AbstractVector)
    I = @index(Global, Cartesian)
    j = I[1]                                    # ring index
    k = CartesianIndex(Base.tail(Tuple(I)))     # all non-horizontal dimensions
    i0 = ring_first[j] + (1 - offset[j])
    reverse_ring!(field, k, i0, ring_first[j] + ring_length[j] - 1)
end

# also allow long names longitude, latitude
Base._reverse!(field::AbstractField, ::Val{:latitude}) = Base._reverse!(field, Val{:lat}())
Base._reverse!(field::AbstractField, ::Val{:longitude}) = Base._reverse!(field, Val{:lon}())
