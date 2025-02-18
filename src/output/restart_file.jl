"""
$(TYPEDSIGNATURES)
A restart file `restart.jld2` with the prognostic variables is written
to the output folder (or current path) that can be used to restart the model.
`restart.jld2` will then be used as initial conditions. The prognostic variables
are bitrounded for compression and the 2nd leapfrog time step is discarded.
Variables in restart file are unscaled."""
function write_restart_file!(
    output::AbstractOutput,
    progn::PrognosticVariables{T},
) where T
    
    # exit immediately if no output or no restart file desired
    output.active && output.write_restart || return nothing  
    
    # move 2nd leapfrog to 1st to compress restart file
    copyto!(progn.vor[1],   progn.vor[2])
    copyto!(progn.div[1],   progn.div[2])
    copyto!(progn.temp[1],  progn.temp[2])
    copyto!(progn.humid[1], progn.humid[2])
    copyto!(progn.pres[1],  progn.pres[2])

    # bitround 1st leapfrog step to output precision
    if T <: Base.IEEEFloat  # currently not defined for other formats...
        round!(progn.vor[1],    7)  # hardcode some defaults for now
        round!(progn.div[1],    7)
        round!(progn.temp[1],  12)
        round!(progn.humid[1], 10)
        round!(progn.pres[1],  14)
    end

    # remove 2nd leapfrog step by filling with zeros
    fill!(progn.vor[2],   0)
    fill!(progn.div[2],   0)
    fill!(progn.temp[2],  0)
    fill!(progn.humid[2], 0)
    fill!(progn.pres[2],  0)

    jldopen(joinpath(output.run_path, "restart.jld2"), "w"; compress=true) do f
        f["prognostic_variables"] = progn
        f["version"] = output.pkg_version
        f["description"] = "Restart file created for SpeedyWeather.jl"
    end
end