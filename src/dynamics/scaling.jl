"""
$(TYPEDSIGNATURES)
Scale the variable `var` inside `progn` with scalar `scale`.
"""
function scale!(
    progn::PrognosticVariables,
    var::Symbol,
    scale::Real,
)
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
    new_scale = scale/progn.scale[]     # undo previous scale and new scale in one go
    scale!(progn,:vor,new_scale)
    scale!(progn,:div,new_scale)
    progn.scale[] = scale               # store scaling information
end

"""
$(TYPEDSIGNATURES)
Undo the radius-scaling of vorticity and divergence from scale!(progn,scale::Real)."""
function unscale!(progn::PrognosticVariables)
    inv_scale = inv(progn.scale[])
    scale!(progn,:vor,inv_scale)
    scale!(progn,:div,inv_scale)
    progn.scale[] = 1                   # set scale back to 1=unscaled
end

"""
$(TYPEDSIGNATURES)
Undo the radius-scaling for any variable. Method used for netcdf output."""
function unscale!(  variable::AbstractArray,
                    scale::Real)
    variable ./= scale
end