# functions to scale the latitude of any grid, in-place
scale_coslat!(  field::AbstractField) = _scale_coslat!(field, power=1)
scale_coslat²!( field::AbstractField) = _scale_coslat!(field, power=2)
scale_coslat⁻¹!(field::AbstractField) = _scale_coslat!(field, power=-1)
scale_coslat⁻²!(field::AbstractField) = _scale_coslat!(field, power=-2)

# and via a deepcopy
scale_coslat(  field::AbstractField) = scale_coslat!(  deepcopy(field))
scale_coslat²( field::AbstractField) = scale_coslat²!( deepcopy(field))
scale_coslat⁻¹(field::AbstractField) = scale_coslat⁻¹!(deepcopy(field))
scale_coslat⁻²(field::AbstractField) = scale_coslat⁻²!(deepcopy(field))

# powers the cosine of latitude
_scale_coslat!(field::AbstractField; power=1) = _scale_lat!(field, cos.(get_lat(field)) .^ power)

"""
$(TYPEDSIGNATURES)
Generic latitude scaling applied to `field` in-place with latitude-like vector `v`."""
function _scale_lat!(field::AbstractField, v::AbstractVector)
    @boundscheck get_nlat(field) == length(v) || throw(DimensionMismatch(field, v))
    
    rings = eachring(field)
    @inbounds for k in eachlayer(field)
        for (j, ring) in enumerate(rings)
            vj = convert(eltype(field), v[j])
            for ij in ring
                field[ij, k] *= vj
            end
        end
    end

    return field
end