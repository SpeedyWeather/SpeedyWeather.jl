# alias functions to scale the latitude of any gridded map A
scale_coslat!(  A::AbstractGrid) = _scale_coslat!(A,power=1)
scale_coslat²!( A::AbstractGrid) = _scale_coslat!(A,power=2)
scale_coslat⁻¹!(A::AbstractGrid) = _scale_coslat!(A,power=-1)
scale_coslat⁻²!(A::AbstractGrid) = _scale_coslat!(A,power=-2)

function _scale_coslat!(A::Grid;power=1) where {Grid<:AbstractGrid}
    coslat = sin.(RingGrids.get_colat(Grid,A.nlat_half))    # sin(colat) = cos(lat)
    coslat .^= power
    return _scale_lat!(A,coslat)
end

"""
$(TYPEDSIGNATURES)
Generic latitude scaling applied to `A` in-place with latitude-like vector `v`."""
function _scale_lat!(A::AbstractGrid{NF},v::AbstractVector) where NF
    @boundscheck get_nlat(A) == length(v) || throw(BoundsError)
    
    rings = eachring(A)
    
    @inbounds for (j,ring) in enumerate(rings)
        vj = convert(NF,v[j])
        for ij in ring
            A[ij] *= vj
        end
    end

    return A
end 

"""
$(TYPEDSIGNATURES)
Scale the variable `var` inside `progn` with scalar `scale`.
"""
function scale!(progn::PrognosticVariables{NF},
                var::Symbol,
                scale::Real) where NF
    if var == :pres
        for pres in progn.pres.timesteps
            pres .*= scale
        end
    else
        for layer in progn.layers
            for step in layer.timesteps
                variable = getfield(step,var)
                variable .*= scale
            end
        end
    end
end

"""
$(TYPEDSIGNATURES)
Scales the prognostic variables vorticity and divergence with
the Earth's radius which is used in the dynamical core."""
function scale!(progn::PrognosticVariables,
                scale::Real)
    scale!(progn,:vor,scale)
    scale!(progn,:div,scale)
    progn.scale[] = scale   # store scaling information
end

"""
$(TYPEDSIGNATURES)
Undo the radius-scaling of vorticity and divergence from scale!(progn,scale::Real)."""
function unscale!(progn::PrognosticVariables)
    scale = progn.scale[]
    scale!(progn,:vor,inv(scale))
    scale!(progn,:div,inv(scale))
    progn.scale[] = 1       # set scale back to 1=unscaled
end

"""
$(TYPEDSIGNATURES)
Undo the radius-scaling for any variable. Method used for netcdf output."""
function unscale!(  variable::AbstractArray,
                    scale::Real)
    variable ./= scale
end