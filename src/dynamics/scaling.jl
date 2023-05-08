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
        for pres in progn.pres.timesteps
            pres .*= s  # pres*s but in-place
        end
    else
        for layer in progn.layers
            for step in layer.timesteps
                variable = getfield(step,var)
                variable .*= s  # var*s but in-place
            end
        end
    end
end

"""
    scale!( progn::PrognosticVariables,
            model::ModelSetup)

Scales the prognostic variables vorticity and divergence with
the Earth's radius which is used in the dynamical core."""
function scale!(progn::PrognosticVariables,
                model::ModelSetup)

    (; radius ) = model.geometry
    scale!(progn,:vor,radius)
    scale!(progn,:div,radius)
end

"""
    unscale!(   progn::PrognosticVariables,
                model::ModelSetup)

Undo the radius-scaling of vorticity and divergence from scale!(progn,model)."""
function unscale!(  progn::PrognosticVariables,
                    model::ModelSetup)

    (; radius ) = model.geometry
    scale!(progn,:vor,inv(radius))
    scale!(progn,:div,inv(radius))
end

"""
    unscale!(   variable::AbstractArray,
                model::ModelSetup)
    
Undo the radius-scaling for any variable. Method used for netcdf output."""
function unscale!(  variable::AbstractArray,
                    model::ModelSetup)
    variable ./= model.geometry.radius
end