"""
$(TYPEDSIGNATURES)
A restart file `restart.jld2` with the prognostic variables is written
to the output folder (or current path) that can be used to restart the model.
`restart.jld2` will then be used as initial conditions. The prognostic variables
are bitrounded for compression and the 2nd leapfrog time step is discarded.
Variables in restart file are unscaled."""
function write_restart_file!(
    output::AbstractOutput,
    progn::PrognosticVariables,
)
    
    # exit immediately if no output or no restart file desired
    output.active && output.write_restart || return nothing  
    
    # move 2nd leapfrog to 1st to compress restart file
    get_step(progn.vor, 1)   .= get_step(progn.vor, 2)
    get_step(progn.div, 1)   .= get_step(progn.div, 2)
    get_step(progn.temp, 1)  .= get_step(progn.temp, 2)
    get_step(progn.humid, 1) .= get_step(progn.humid, 2)
    get_step(progn.pres, 1)  .= get_step(progn.pres, 2)

    # bitround 1st leapfrog step to output precision
    if eltype(progn) <: Base.IEEEFloat  # currently not defined for other formats...
        round!(get_step(progn.vor, 1),    7)      # hardcode some defaults for now
        round!(get_step(progn.div, 1),    7)
        round!(get_step(progn.temp, 1),  12)
        round!(get_step(progn.humid, 1), 10)
        round!(get_step(progn.pres, 1),  14)
    end

    # remove 2nd leapfrog step by filling with zeros
    fill!(get_step(progn.vor, 2),   0)
    fill!(get_step(progn.div, 2),   0)
    fill!(get_step(progn.temp, 2),  0)
    fill!(get_step(progn.humid, 2), 0)
    fill!(get_step(progn.pres, 2),  0)

    jldopen(joinpath(output.run_path, "restart.jld2"), "w"; compress=true) do f
        f["prognostic_variables"] = progn
        f["version"] = output.pkg_version
        f["description"] = "Restart file created for SpeedyWeather.jl"
    end
end