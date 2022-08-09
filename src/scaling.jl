# alias functions to scale the latitude of lat-lon map `A`
scale_coslat!(  A::AbstractMatrix,G::Geometry) = _scale_lat!(A,G.coslat)
scale_coslat²!( A::AbstractMatrix,G::Geometry) = _scale_lat!(A,G.coslat²)
scale_coslat⁻¹!(A::AbstractMatrix,G::Geometry) = _scale_lat!(A,G.coslat⁻¹)
scale_coslat⁻²!(A::AbstractMatrix,G::Geometry) = _scale_lat!(A,G.coslat⁻²)

"""
    _scale_lat!(A::AbstractMatrix{NF},v::AbstractVector) where {NF<:AbstractFloat}

Generic latitude scaling applied to `A` in-place with latitude-like vector `v`."""
function _scale_lat!(A::AbstractMatrix{NF},v::AbstractVector) where {NF<:AbstractFloat}
    nlon,nlat = size(A)
    @boundscheck nlat == length(v) || throw(BoundsError)

    @inbounds for j in 1:nlat
        vj = convert(NF,v[j])
        for i in 1:nlon
            A[i,j] *= vj
        end
    end
end 

"""
    scale!( progn::PrognosticVariables{NF},
            var::Symbol,
            s::Number) where NF

Scale the variable `var` inside `progn` with scalar `s`.
"""
function scale!(progn::PrognosticVariables{NF},
                var::Symbol,
                s::Number) where NF

    s_NF = convert(Complex{NF},s)
    
    if var == :pres     # surface pressure is not stored in layers
        for leapfrog_step in progn.pres
            leapfrog_step *= s_NF
        end
    else
        for layer in progn.layers
            for leapfrog_step in layer.leapfrog
                @eval $leapfrog_step.$var .*= $s_NF
            end
        end
    end
end