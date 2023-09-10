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