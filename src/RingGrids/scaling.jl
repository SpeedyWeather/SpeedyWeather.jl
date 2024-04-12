# alias functions to scale the latitude of any gridded map A
scale_coslat!(  A::AbstractGridArray) = _scale_coslat!(A, power=1)
scale_coslat²!( A::AbstractGridArray) = _scale_coslat!(A, power=2)
scale_coslat⁻¹!(A::AbstractGridArray) = _scale_coslat!(A, power=-1)
scale_coslat⁻²!(A::AbstractGridArray) = _scale_coslat!(A, power=-2)

function _scale_coslat!(A::Grid; power=1) where {Grid<:AbstractGridArray}
    coslat = sin.(get_colat(Grid, A.nlat_half))    # sin(colat) = cos(lat)
    coslat .^= power
    return _scale_lat!(A, coslat)
end

"""
$(TYPEDSIGNATURES)
Generic latitude scaling applied to `A` in-place with latitude-like vector `v`."""
function _scale_lat!(A::AbstractGridArray{NF}, v::AbstractVector) where NF
    @boundscheck get_nlat(A) == length(v) || throw(BoundsError)
    @inbounds for (j, ring) in enumerate(eachring(A))
        vj = convert(NF, v[j])
        for ij in ring
            A[ij] *= vj
        end
    end
    return A
end 