"""$(TYPEDSIGNATURES)
Scales the prognostic variables vorticity and divergence with
the Earth's radius which is used in the dynamical core."""
@propagate_inbounds function scale_prognostic!(vars::Variables, scale::Real)
    progn = vars.prognostic             # for convenience
    new_scale = scale / progn.scale[]   # undo previous scale and new scale in one go
    haskey(progn, :vorticity) && (progn.vorticity .*= new_scale)
    haskey(progn, :divergence) && (progn.divergence .*= new_scale)
    # no need to actually scale the diagnostic variables as they will be
    # overwritten by the transform of the prognostic variables anyway
    progn.scale[] = scale               # store scaling information
    return vars
end

function scale_tendencies!(vars::Variables, model::AbstractModel)
    scale = vars.prognostic.scale[]
    (; tendencies) = vars
    TS = model.time_stepping

    # spectral
    for varname in tendency_names(vars)
        var = get_tendency_step(getfield(tendencies, varname), TS, DummyParameterization())
        scale!(var, scale)
    end

    # grid
    for varname in tendency_and_uv_names(vars)
        var = get_tendency_step(getfield(tendencies.grid, varname), TS, DummyParameterization())
        scale!(var, scale)
    end

    # tracers
    for varname in tracer_tendency_names(vars)
        var = get_tendency_step(getfield(tendencies.grid_tracers, varname), TS, DummyParameterization())
        scale!(var, scale)
    end
    return nothing
end

function unscale_tendencies!(vars::Variables)
    scale = vars.prognostic.scale[]
    (; tendencies) = vars

    # spectral
    for varname in tendency_names(vars)
        unscale!(getfield(tendencies, varname), scale)
    end

    # grid
    for varname in tendency_and_uv_names(vars)
        unscale!(getfield(tendencies.grid, varname), scale)
    end

    # tracers
    for varname in tracer_tendency_names(vars)
        unscale!(getfield(tendencies.grid_tracers, varname), scale)
    end
    return nothing
end


"""$(TYPEDSIGNATURES)
Undo the radius-scaling of vorticity and divergence from `scale_prognostic!(vars, scale::Real)`."""
function unscale!(vars::Variables)
    progn = vars.prognostic             # for convenience
    inv_scale = inv(progn.scale[])
    haskey(progn, :vorticity) && (progn.vorticity .*= inv_scale)
    haskey(progn, :divergence) && (progn.divergence .*= inv_scale)

    # and the corresponding grid variables if they exist
    haskey(vars.grid, :vorticity) && (vars.grid.vorticity .*= inv_scale)
    haskey(vars.grid, :divergence) && (vars.grid.divergence .*= inv_scale)

    # also unscale tendencies
    unscale_tendencies!(vars)

    progn.scale[] = 1                   # set scale back to 1=unscaled
    return vars
end

"""$(TYPEDSIGNATURES)
Scale the variable `var` with scalar `scale`."""
@inline function scale!(
        variable::Union{LowerTriangularArray, AbstractField},
        scale::Real
    )
    variable.data .*= scale
    return variable
end

"""$(TYPEDSIGNATURES)
Undo the scaling of the variable `var` with scalar `scale`."""
@inline function unscale!(
        variable::Union{LowerTriangularArray, AbstractField},
        scale::Real
    )
    return scale!(variable, inv(scale))
end
