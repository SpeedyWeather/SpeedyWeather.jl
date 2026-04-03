"""$(TYPEDSIGNATURES)
Scales the prognostic variables vorticity and divergence with
the Earth's radius which is used in the dynamical core."""
@propagate_inbounds function scale_prognostic!(vars::Variables, scale::Real)
    progn = vars.prognostic             # for convenience
    new_scale = scale / progn.scale[]   # undo previous scale and new scale in one go
    haskey(progn, :vor) && (progn.vor .*= new_scale)
    haskey(progn, :div) && (progn.div .*= new_scale)
    # no need to actually scale the diagnostic variables as they will be
    # overwritten by the transform of the prognostic variables anyway
    progn.scale[] = scale               # store scaling information
    return vars
end

scale!(var::AbstractField, scale::Real) = (var .*= scale)

"""$(TYPEDSIGNATURES)
Scale the tendencies of u, v, temp, humid with scalar `scale`.
Intended use to scale the tendencies of the parameterizations
by the radius for the dynamical core."""
@propagate_inbounds function scale_tendencies!(vars::NamedTuple, scale::Real)
    haskey(vars, :u) && (vars.u .*= scale)
    haskey(vars, :v) && (vars.v .*= scale)
    haskey(vars, :temp) && (vars.temp .*= scale)
    haskey(vars, :humid) && (vars.humid .*= scale)
    return nothing
end

@propagate_inbounds function scale_tendencies!(ij, vars::NamedTuple, scale::Real)
    haskey(vars, :u) && (
        for k in eachlayer(vars.u)
            vars.u[ij, k] *= scale
        end
    )
    haskey(vars, :v) && (
        for k in eachlayer(vars.v)
            vars.v[ij, k] *= scale
        end
    )
    haskey(vars, :temp) && (
        for k in eachlayer(vars.temp)
            vars.temp[ij, k] *= scale
        end
    )
    haskey(vars, :humid) && (
        for k in eachlayer(vars.humid)
            vars.humid[ij, k] *= scale
        end
    )
    return nothing
end

"""$(TYPEDSIGNATURES)
Undo the radius-scaling of vorticity and divergence from `scale_prognostic!(vars, scale::Real)`."""
function unscale!(vars::Variables)
    progn = vars.prognostic             # for convenience
    inv_scale = inv(progn.scale[])
    haskey(progn, :vor) && (progn.vor .*= inv_scale)
    haskey(progn, :div) && (progn.div .*= inv_scale)

    # and the corresponding grid variables if they exist
    haskey(vars.grid, :vor) && (vars.grid.vor .*= inv_scale)
    haskey(vars.grid, :div) && (vars.grid.div .*= inv_scale)

    # TODO unscale the tendencies too?

    progn.scale[] = 1                   # set scale back to 1=unscaled
    return vars
end

"""
$(TYPEDSIGNATURES)
Scale the variable `var` with scalar `scale`.
"""
@propagate_inbounds function scale!(
        variable::Union{LowerTriangularArray, Field},
        scale::Real
    )
    return variable.data .*= scale
end

"""
$(TYPEDSIGNATURES)
Undo the scaling of the variable `var` with scalar `scale`.
"""
@propagate_inbounds function unscale!(
        variable::Union{LowerTriangularArray, Field},
        scale::Real
    )
    return variable.data ./= scale
end
