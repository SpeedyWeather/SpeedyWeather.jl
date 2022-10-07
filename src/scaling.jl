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
    @boundscheck get_nlat(A) == length(v) || throw(BoundsError)
    
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
    
    if var == :pres     # surface pressure is not stored in layers
        for pres in progn.pres.leapfrog
            scale!(pres,s)              # pres*s but in-place
        end
    else
        for layer in progn.layers
            for leapfrog_step in layer.leapfrog
                variable = getfield(leapfrog_step,var)
                scale!(variable,s)      # var*s but in-place
            end
        end
    end
end

"""
    scale!(progn::PrognosticVariables{NF}) where NF 

Scales all variables within `progn` with the scale set when initializing it.
"""
function scale!(progn::PrognosticVariables{NF}) where NF 
    for (var,scale) in progn.scale 
        scale!(progn, var, scale)
    end 
end 

"""
    unscale!(progn::PrognosticVariables{NF}) where NF 

Scales  all variables within `progn` with the inverse of scale set when initializing it.
"""
function unscale!(progn::PrognosticVariables{NF}) where NF 
    for (var,scale) in progn.scale 
        scale!(progn, var, inv(scale))
    end
end  

"""
    unscale(progn::PrognosticVariables{NF}) where NF 

Scales all variables within `progn` with the inverse of scale set when initializing it.
"""
function unscale(progn::PrognosticVariables{NF}) where NF 
    unscaled_progn = copy(progn)
    for (var,scale) in progn.scale 
        scale!(unscaled_progn, var, inv(scale))
    end
    return unscaled_progn
end  