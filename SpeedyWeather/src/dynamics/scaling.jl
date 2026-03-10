"""$(TYPEDSIGNATURES)
Scales the prognostic variables vorticity and divergence with
the Earth's radius which is used in the dynamical core."""
@propagate_inbounds function scale!(vars::Variables, scale::Real)
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
Scale the tendencies inside `diagn` with scalar `scale`.
Intended use to scale the tendencies of the parameterizations
by the radius for the dynamical core."""
@propagate_inbounds function scale!(ij, vars::NamedTuple, scale::Real)
    for var in vars
        for k in eachlayer(var)
            var[ij, k] *= scale
        end
    end
    return nothing
end

"""$(TYPEDSIGNATURES)
Undo the radius-scaling of vorticity and divergence from `scale!(vars, scale::Real)`."""
function unscale!(vars::Variables)
    progn = vars.prognostic             # for convenience
    inv_scale = inv(progn.scale[])
    haskey(progn, :vor) && (progn.vor .*= inv_scale)
    haskey(progn, :div) && (progn.div .*= inv_scale)
    progn.scale[] = 1                   # set scale back to 1=unscaled
    return vars
end
