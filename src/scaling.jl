# alias functions to scale the latitude of any gridded map A
scale_coslat!(  A::AbstractGrid,G::Geometry) = _scale_lat!(A,G.coslat)
scale_coslat²!( A::AbstractGrid,G::Geometry) = _scale_lat!(A,G.coslat²)
scale_coslat⁻¹!(A::AbstractGrid,G::Geometry) = _scale_lat!(A,G.coslat⁻¹)
scale_coslat⁻²!(A::AbstractGrid,G::Geometry) = _scale_lat!(A,G.coslat⁻²)

# matrix versions used for output
scale_coslat!(  A::AbstractMatrix,G::Geometry) = A.*G.coslat'
scale_coslat²!( A::AbstractMatrix,G::Geometry) = A.*G.coslat²'
scale_coslat⁻¹!(A::AbstractMatrix,G::Geometry) = A.*G.coslat⁻¹'
scale_coslat⁻²!(A::AbstractMatrix,G::Geometry) = A.*G.coslat⁻²'

"""
    _scale_lat!(A::AbstractGrid,v::AbstractVector)

Generic latitude scaling applied to `A` in-place with latitude-like vector `v`."""
function _scale_lat!(A::AbstractGrid{NF},v::AbstractVector) where {NF<:AbstractFloat}
    @boundscheck length(get_nlat(A)) == length(v) || throw(BoundsError)
    
    rings = eachring(A)
    
    @inbounds for (j,ring) in enumerate(rings)
        vj = convert(NF,v[j])
        for ij in ring
            A[ij] *= vj
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