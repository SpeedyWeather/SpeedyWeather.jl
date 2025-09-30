"""
$(TYPEDSIGNATURES)
Scale the variable `var` inside `progn` with scalar `scale`.
"""
function scale!(
    progn::PrognosticVariables,
    var::Symbol,
    scale::Real,
)
    var = getfield(progn, var)
    var .*= scale
end

"""
$(TYPEDSIGNATURES)
Scale the variable `var` inside `diagn` with scalar `scale`.
"""
function scale!(
    diagn::DiagnosticVariables,
    var::Symbol,
    scale::Real,
)
    variable = getfield(diagn.grid, var)
    variable .*= scale
end

"""
$(TYPEDSIGNATURES)
Scales the prognostic variables vorticity and divergence with
the Earth's radius which is used in the dynamical core."""
function scale!(progn::PrognosticVariables,
                diagn::DiagnosticVariables,
                scale::Real)
                
    new_scale = scale/progn.scale[]     # undo previous scale and new scale in one go
    scale!(progn, :vor, new_scale)
    scale!(progn, :div, new_scale)
    progn.scale[] = scale               # store scaling information
    diagn.scale[] = scale               # store scaling information
    # no need to actually scale the diagnostic variables as they will be
    # overwritten by the transform of the prognostic variables anyway
end

"""
$(TYPEDSIGNATURES)
Undo the radius-scaling of vorticity and divergence from scale!(progn, scale::Real)."""
function unscale!(progn::PrognosticVariables)
    inv_scale = inv(progn.scale[])
    scale!(progn, :vor, inv_scale)
    scale!(progn, :div, inv_scale)
    progn.scale[] = 1                   # set scale back to 1=unscaled
end

"""
$(TYPEDSIGNATURES)
Undo the radius-scaling of vorticity and divergence from scale!(diagn, scale::Real)."""
function unscale!(diagn::DiagnosticVariables)
    inv_scale = inv(diagn.scale[])
    scale!(diagn, :vor_grid, inv_scale)
    scale!(diagn, :div_grid, inv_scale)
    diagn.scale[] = 1                   # set scale back to 1=unscaled
end

"""
$(TYPEDSIGNATURES)
Undo the radius-scaling for any variable. Method used for netcdf output."""
function unscale!(  variable::AbstractArray,
                    scale::Real)
    variable ./= scale
end