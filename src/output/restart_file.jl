export RestartFile

"""RestartFile callback. Writes a restart file as .jld2.
Options are $(TYPEDFIELDS)"""
@kwdef struct RestartFile <: AbstractCallback
    "File name for restart file, should end with .jld2"
    filename::String = "restart.jld2"

    "Path for restart file, uses model.output.run_path if not specified"
    path::String = ""

    "Apply lossless compression in JLD2?"
    compress::Bool = true
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
    (; compress, filename) = restart

    # use output run path if not specified
    path = restart.path == "" ? model.output.run_path : restart.path

    jldopen(joinpath(path, filename), "w"; compress) do f
        f["prognostic_variables"] = progn
        f["version"] = model.output.pkg_version
        f["description"] = "Restart file created by SpeedyWeather.jl"
    end
end