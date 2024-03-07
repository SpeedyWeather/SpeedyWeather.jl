# alias functions to scale the latitude of any gridded map A
scale_coslat!(  A::AbstractGrid) = _scale_coslat!(A, power=1)
scale_coslat²!( A::AbstractGrid) = _scale_coslat!(A, power=2)
scale_coslat⁻¹!(A::AbstractGrid) = _scale_coslat!(A, power=-1)
scale_coslat⁻²!(A::AbstractGrid) = _scale_coslat!(A, power=-2)

function _scale_coslat!(A::Grid; power=1) where {Grid<:AbstractGrid}
    coslat = sin.(get_colat(Grid, A.nlat_half))    # sin(colat) = cos(lat)
    coslat .^= power
    return _scale_lat!(A, coslat)
end

"""
$(TYPEDSIGNATURES)
Generic latitude scaling applied to `A` in-place with latitude-like vector `v`."""
function _scale_lat!(A::AbstractGrid{NF}, v::AbstractVector) where NF
    @boundscheck get_nlat(A) == length(v) || throw(BoundsError)
    
    rings = eachring(A)
    
    @inbounds for (j, ring) in enumerate(rings)
        vj = convert(NF, v[j])
        for ij in ring
            A[ij] *= vj
        end
    end

    return A
end 