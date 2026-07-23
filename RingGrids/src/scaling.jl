# functions to scale the latitude of any grid, in-place
scale_coslat!(field::AbstractField) = _scale_coslat!(field, power = 1)
scale_coslat²!(field::AbstractField) = _scale_coslat!(field, power = 2)
scale_coslat⁻¹!(field::AbstractField) = _scale_coslat!(field, power = -1)
scale_coslat⁻²!(field::AbstractField) = _scale_coslat!(field, power = -2)

# and via a deepcopy
scale_coslat(field::AbstractField) = scale_coslat!(deepcopy(field))
scale_coslat²(field::AbstractField) = scale_coslat²!(deepcopy(field))
scale_coslat⁻¹(field::AbstractField) = scale_coslat⁻¹!(deepcopy(field))
scale_coslat⁻²(field::AbstractField) = scale_coslat⁻²!(deepcopy(field))

# powers the cosine of latitude; `get_lat` returns a CPU vector by design so far,
# so move it onto the field's architecture before launching the kernel.
function _scale_coslat!(field::AbstractField; power = 1)
    coslat_pow = cos.(get_lat(field)) .^ power
    return _scale_lat!(field, on_architecture(architecture(field), coslat_pow))
end

"""
$(TYPEDSIGNATURES)
Generic latitude scaling applied to `field` in-place with latitude-like vector `v`."""
function _scale_lat!(field::AbstractField, v::AbstractVector)
    @boundscheck get_nlat(field) == length(v) || throw(DimensionMismatch(field, v))

    arch = architecture(field)
    launch!(arch, RingGridWorkOrder, size(field), scale_lat_kernel!, field, v, whichring(field.grid))
    return field
end

@kernel inbounds = true function scale_lat_kernel!(field, v, whichring)
    I = @index(Global, Cartesian)
    ij = I[1]                                   # grid point index
    k = CartesianIndex(Base.tail(Tuple(I)))     # all non-horizontal dimensions
    j = whichring[ij]   # get ring index for grid point ij
    field[ij, k] *= v[j]
end
