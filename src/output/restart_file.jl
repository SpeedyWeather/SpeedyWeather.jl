export RestartFile

"""RestartFile callback. Writes a restart file as .jld2.
Options are $(TYPEDFIELDS)"""
@kwdef mutable struct RestartFile <: AbstractCallback
    "[OPTION] File name for restart file, should end with .jld2"
    filename::String = "restart.jld2"

    "[OPTION] Path for restart file, uses model.output.run_path if not specified"
    path::String = ""

    "[OPTION] Apply lossless compression in JLD2?"
    compress::Bool = true

    "[OPTION] Only write with model.output.active = true?"
    write_only_with_output::Bool = true

    "[DERIVED] package version to track possible incompatibilities"
    pkg_version::VersionNumber = isnothing(pkgversion(SpeedyWeather)) ? v"0.0.0" : pkgversion(SpeedyWeather)
end

initialize!(::RestartFile, args...) = nothing
callback!(::RestartFile, args...) = nothing

"""
$(TYPEDSIGNATURES)
A restart file `restart.jld2` with the prognostic variables is written
to the output folder (or current path) that can be used to restart the model.
`restart.jld2` will then be used as initial conditions.
Variables in restart file are not scaled."""
function finalize!(
    restart::RestartFile,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    # escape in case of no output
    restart.write_only_with_output && (model.output.active || return nothing)

    (; compress, filename) = restart

    # use output run path if not specified
    path = restart.path == "" ? model.output.run_path : restart.path
    mkpath(path)

    jldopen(joinpath(path, filename), "w"; compress) do f
        f["prognostic_variables"] = progn
        f["version"] = restart.pkg_version
        f["description"] = "Restart file created by SpeedyWeather.jl"
    end

    model.output.active || @info "Restart file written to $(joinpath(path, filename)) although output=false"
end