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
Reverse the direction of every ring of `field` in-place (west-east mirror),
as called via `reverse!(field, dims=:lon)`."""
function Base._reverse!(field::AbstractField, ::Val{:lon})
    arch = architecture(field)
    ring_first, ring_length = eachring_on_architecture(arch, field.grid)

    # one thread per ring j (and layer k) as the in-place reversal is sequential within a ring
    worksize = (get_nlat(field), size(field)[2:end]...)
    launch!(arch, ArrayWorkOrder, worksize, reverse_lon_kernel!, field, ring_first, ring_length)
    return field
end

@kernel inbounds = true function reverse_lon_kernel!(field, ring_first, ring_length)
    I = @index(Global, Cartesian)
    j = I[1]                                    # ring index
    k = CartesianIndex(Base.tail(Tuple(I)))     # all non-horizontal dimensions
    i0 = ring_first[j]
    reverse_ring!(field, k, i0, i0 + ring_length[j] - 1)
end

# also allow long names longitude, latitude
Base._reverse!(field::AbstractField, ::Val{:latitude}) = Base._reverse!(field, Val{:lat}())
Base._reverse!(field::AbstractField, ::Val{:longitude}) = Base._reverse!(field, Val{:lon}())
